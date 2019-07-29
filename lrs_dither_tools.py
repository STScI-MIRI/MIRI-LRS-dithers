import numpy as np
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
from os import access, R_OK
import glob
from astropy.table import Table, Column
import astropy.io.ascii as ascii
import datetime
import sys
import pdb


import pysiaf
import miricoord
import miricoord.miricoord.imager.mirim_tools as mt

plt.close('all')


class LRSPattern(object):
	
	'''
	This code defines a class for an LRS dither pattern. It is initialised from a table of positions (this should be an astropy table), or from an ascii file that uses the standard pattern template (see LRS_pattern_template.txt). 
	
	Initialisation:
	---------------
	- file: a filename. If the file passed is not in the standard template, the initialization will fail as the metadata will likely not be found.
	- t: an Astropy table with position index, x and y location. if the class is initialised with a table, metadata must be provided with the call to LRSPattern. 
	
	NOTE: a input file will have priority over an input table. if both are provided, the table will be IGNORED.
	
	
	Attributes:
	-----------
	- npts: number of pointings in the pattern (integer)
	- mode: the LRS mode - 'slit' or 'slitless', or both. If both are written then 2 patterns will be generated. (string)
	- frame: the coordinate reference frame of the coordinates (string). Options:
			* 'tel' for telescope coordinates (aka v2v3), in arcsec
			* 'det-abs' for *absolute* detector coordinates, in pixels
			* 'det-rel' for *relative* detector coordinates, in pixels. to translate between relative and absolute coordinates, the code will look at the 'ref' attribute and add the required absolute coordinates.
			* 'idl' for ideal coordinates, in arcsec. NOTE: slit and slitless mode have distinct ideal coordinate systems
	- car: the commissioning activity request(s) this pattern is used for. this should be specified in the format 'CAR-xxx'
	- cal: the calibration activity request(s) this pattern is used for. this should be specified in the format 'CAL-xxx'
	- ref: the reference point for this pattern. can be a single point or a list. choose from 'c', 'nod1', 'nod2'.
	- notes: any comments in a string format. PLEASE USE to add relevant information to the pattern!
	- nref: number of reference points
	
	
	Methods:
	--------
	TO DO! THESE ARE IDEAS!
	- convert coordinates: convert between frames
	- plot: visualize the pattern
	- regenerate: regenerate in case of new SIAF update.
	- save: save to ascii file
	
	'''
	
	def __init__(self, table=None, file=None, mode=None, frame=None, car=None, cal=None, name=None, ref=None, notes=None):
		
		# check that either a table or a file are provided
		assert (table is not None) or (file is not None), "You must provide either a table or a pattern file to initialize an LRSPattern instance"
		
		if (file is not None):
			
			# if a filename is provided, the other keywords will be ignored and we'll attempt to initialise the pattern from the file
			
			t = ascii.read(file)
			# t is an astropy table. the comment lines are placed in the table metadata. populate the class attributes from the metadata
			header = ascii.read(t.meta['comments'], delimiter=':', format='no_header', names=['key', 'val'])
			self.patt = t
			self.mode = header['val'][header['key']=='Mode']
			self.npts = len(t)
			self.car = header['val'][header['key']=='CARs']
			self.cal = header['val'][header['key']=='CALs']
			self.notes = header['val'][header['key']=='Comments']
			self.ref = header['val'][header['key']=='Reference']
			self.frame = header['val'][header['key']=='Frame']
			self.name = file.split('.')[0]
			
			
			# first check if there are multiple reference positions:
			ref_tmp = (self.ref[0]).split(',')
			# strip out whitespaces in case there are any:
			refs = [r.strip() for r in ref_tmp]
			self.nref = len(refs)
			
			# if there is only 1 reference position, do nothing. if there are N>1 positions, copy the x and y columns N-1 times to make space for the extra pointings. NOTE: if there ar emultiple reference positions, the coordinate frame shouls be relative not absolute - otherwise that wouldn't make sense.
			
			if (self.nref > 1):
				assert ('rel' in self.frame[0]), "Cannot have multiple references for the pattern with an absolute coordinate frame"
				i = 0
				while (i < self.nref-1):
					self.patt.add_column(t['x'].copy(), name='x{0}'.format(i+1))
					self.patt.add_column(t['y'].copy(), name='y{0}'.format(i+1))	
					i += 1		
				
				
			
		else:
			
			# if no filename is provided, initialise from the keywords
			
			self.patt = t
			self.npts = len(t)
			self.mode = mode
			self.car = car
			self.cal = cal
			self.notes = notes
			self.frame = frame
			self.reps = reps
			self.name = name
			self.ref = ref
			
			# first check if there are multiple reference positions:
			ref_tmp = (self.ref[0]).split(',')
			# strip out whitespaces in case there are any:
			refs = [r.strip() for r in ref_tmp]
			self.nref = len(refs)
			
			if (self.nref > 1):
				assert ('rel' in self.frame[0]), "Cannot have multiple references for the pattern with an absolute coordinate frame"
				i = 0
				while (i < self.nref-1):
					self.patt.add_column(t['x'].copy(), name='x{0}'.format(i+1))
					self.patt.add_column(t['y'].copy(), name='y{0}'.format(i+1))	
					i += 1	
					
					
		#self.run_checks()
		
		
	def run_checks(self):
		
		'''In this function we will run some sanity checks on the coordinates and reference frames provided, so that we don't waste time or create unrealistic patterns.
		
		Parameters:
		-----------
		the pattern, an LRSPattern instance
		
		Output:
		-------
		None
		'''
		
		xcols = [xx for xx in self.patt.colnames if 'x' in xx]
		ycols = [yy for yy in self.patt.colnames if 'y' in yy]
		
		if (self.frame=='det-abs'):
			
			for c in self.patt.colnames[1:]:
				assert (np.all(self.patt[c]) > 0.), "Coordinate frame/values inconsistency: Coordinates must be positive in an absolute frame"
				
			# for x columns the values have to lie below 1032, for y columns below 1024. strictly speaking it's below 1032.5, but let's assume that if we're trying to point half a pixel from the edge, something's gone wrong.
			for c in self.patt[xcols]:
				assert (np.all(self.patt[c]) <= 1032.), "Coordinate frame/values inconsistency: Pixel coordinates out of bounds"
			for c in self.patt[ycols]:
				assert (np.all(self.patt[c]) <= 1024.), "Coordinate frame/values inconsistency: Pixel coordinates out of bounds"
		
		
		elif (self.frame=='det-rel'):
			 #with the relative detector coordinates values can be negative but they still shouldn't exceed the array size.
			
			for c in self.patt[xcols]:
				assert (np.all(self.patt[c]) <= 1032.), "Coordinate frame/values inconsistency: Pixel coordinates out of bounds"
			for c in self.patt[ycols]:
				assert (np.all(self.patt[c]) <= 1024.), "Coordinate frame/values inconsistency: Pixel coordinates out of bounds"
			#pdb.set_trace()

		# for frame 'tel', check that the coordinates are within the MIRI Imager detector
		elif (self.frame=='tel'):
			for c in self.patt[xcols]:
				assert (np.all(self.patt[c]) >= -486.) and (np.all(self.patt[c]) <= -381.), "Coordinate frame/values inconsistency: Telescope coordinates out of range of the MIRI Imager detector"
			for c in self.patt[ycols]:
				assert (np.all(self.patt[c]) >= -436.) and (np.all(self.patt[c]) <= -314.), "Coordinate frame/values inconsistency: Telescope coordinates out of range of the MIRI Imager detector"
			
			
		
		
		
		
			
	
	def to_absolute(self):
		
		'''
		This function will change the coordinates from relative to absolute. This is ONLY applicable for patterns defined in the 'det-rel' frame; it will take the pattern, identify the reference coordinates for the given mode, and calculate the absolute detector pixels.
		
		Parameters:
		-----------
		None
		
		Returns:
		--------
		- the input pattern with updated positions, and the reference frame set to 'det-abs' in the metadata
		
		
		'''
		
		# first check that the input pattern is defined in relative coordinates
		assert ('rel' in self.frame[0]), "The input pattern is not in relative coordinates!"
		
		print('Converting pattern to absolute coordinates.....')
		# if the input coordinates are relative, load the coordinates so we can do the translation
		coords = lrs_gencoords(mode=self.mode, frame=self.frame[0], verbose=False)
		
		# convert the reference attribute to a list
		ref_tmp = (self.ref[0]).split(',')
		# strip out whitespaces in case there are any:
		refs = [r.strip() for r in ref_tmp]
		
		# the (number of columns - 1) needs to be == (2 * the number of entries in the reference list)
		n_coordcols = len(self.patt.colnames) - 1
		assert (n_coordcols == 2*self.nref), "number of columns in the pattern doesn't match the number of reference points listed"
		

		for i, r in enumerate(refs):
			if (i==0):
				self.patt['x'] += coords[r]['x']
				self.patt['y'] += coords[r]['y']
			else:
				colx = 'x{0}'.format(i)
				coly = 'y{0}'.format(i)
				self.patt[colx] += coords[r]['x']
				self.patt[coly] += coords[r]['y']
		
		# now update the frame attribute from relative to absolute:
		self.frame = '{0}-abs'.format(self.frame[0][:3])
		
		
			
	
	def plot(self, out=None):
		
		'''
		This function will visualize the pattern, in the coordinate system provided. If the 'ref' attribute contains more than 1 entry, then the pattern will be duplicated and translated to all locations.
		
		Parameters:
		-----------
		- out: 	if a filename is provided, the plot will be saved to this file. Default = None. (string)
		
		Future work:
		------------
		Specify the coordinate frame and perform the conversion before plotting.
		
		
		'''
		
		# generate the key coordinates
		coords = lrs_gencoords(mode=self.mode, frame=self.frame[:3])
		
		# figure out if there are multiple refrence positions
		ref_tmp = (self.ref[0]).split(',')
		# strip out whitespaces in case there are any:
		refs = [r.strip() for r in ref_tmp]
		nref = len(refs)
		print(refs)
		pdb.set_trace()
		
		
		
		fig, ax = plt.subplots(figsize=[12,4])
		
		if (self.frame == 'det-abs') or (self.frame == 'det-rel'):
			units = 'pixels'
		else:
			units = 'arcsec'
		
		# check whether the pattern is for slit or slitless
		if (self.mode == 'slit'):
			coords = lrs_gencoords(mode='slit', frame=pltframe)
			
			# if the coordinates are in relative pixel coordinates, translate the pattern 
			if (self.frame == 'det-rel'):
				self.to_absolute()
				
			cornersx = np.array((coords['ll']['x'], coords['ul']['x'], coords['ur']['x'], coords['lr']['x']))
			cornersy = np.array((coords['ll']['y'], coords['ul']['y'], coords['ur']['y'], coords['lr']['y']))
			corners = np.stack((cornersx, cornersy), axis=-1)
			# create the matplotlib patch for plotting
			slit_rect = Polygon(corners, closed=True, lw=2., color='g', fill=False, label='slit edge')
			ax.scatter(cornersx, cornersy, marker='+', color='g')
			ax.scatter(coords['c']['x'], coords['c']['y'], marker='o', color='r', label='slit centre')
			for n in ['nod1', 'nod2']:
				ax.scatter(coords[n]['x'], coords[n]['y'], marker='o', edgecolor='r', facecolor='white')
			ax.add_patch(slit_rect)
			
		elif (self.mode == 'slitless'):
			coords = lrs_gencoords(mode='slitless', frame=pltframe)
			
			# as above, if the points are given in relative pixel values, translate the pattern
			if (self.frame == 'det-rel'):
				#self.patt['x'] += coords['c']['x']
				#self.patt['y'] += coords['c']['y']
				self.to_absolute()
		
		if (len(self.ref) == 1):
			
			# if there's only 1 reference point, we don't need to expand the pattern to account for different reference points
			# if there's only 1 reference point we assume it's the centre/nominal pointing location
			
			#indx = np.arange(len(self.patt))
			pts = ax.scatter(self.patt['x'], self.patt['y'], marker='x', c=self.patt['Pointing'], cmap='plasma', label='pointings ({})'.format(self.npts))
			cbar = fig.colorbar(pts, ax=ax)
			cbar.set_label('pointing index')
			ax.grid(color='k', linestyle='-', linewidth=0.5, alpha=0.5)
			ax.set_xlabel(units)
			ax.set_ylabel(units)
			ax.legend(loc='upper right')
				 
		
		ax.set_title(self.name)
		

		if (out is not None):
			fig.savefig(out)
    
		fig.show()
    
		return fig
		
		


#----------------------------------------------------------------------------	
def lrs_gencoords(mode='slit', frame='tel', plot=False, verbose=False):
	
	'''
	Function that returns a dictionary with the key coordinates for the MIRI LRS. 
	
	Parameters:
	-----------
	- mode: 'slit' (default) or 'slitless'. this will determine which set of coordinates is returned
	- frame: the reference frame in which the coordinates are specified. Options: 'det-abs' (detector [pixels], absolute), 'det-rel' (detector [pixels], relative) 'idl' (Ideal [arcsec]), 'tel' (v2v3 [arcsec])
	- verbose: prints some output to screen (boolean, default False)
	
	
	Returns:
	-------
	- coord: a dictionary with the relevant coordinates. For SLIT, this will contain the slit centre and slit corner locations in the requested frame and units (see options above). For SLITLESS, this contains only the nominal pointing position.
	
	'''
	
	# checks that the parameters provided are valid
	assert (mode in ['slit', 'slitless']), "Mode not recognised. Options are: 'slit', 'slitless'. "
	assert (frame in ['tel', 'idl', 'det-abs', 'det-rel']), "Coordinate frame not recognised. Options are: 'det-rel' (detector relative), 'det-abs' (detector absolute), 'idl' (ideal), 'tel' (telescope/v2v3)."
	
	# the read_aperture function (below) converts 'slit' or 'slitless' to the appropriate SIAF aperture name.
	ap = read_aperture(mode=mode)
	
	# the slit corner coordinates are hard coded unfortunately (in detector units):
	if (mode == 'slit'):
		# populate the dictionary with the slit corner and centre coordinates, in PIXELS
		coord_dict_det = {'ll': {'x': 304.77, 'y': 298.38},
			    	'ul': {'x': 304.77, 'y': 303.03}, 
			       	'ur': {'x': 347.49, 'y': 303.03},
			       	'lr': {'x': 347.49, 'y': 298.38},
					'c': {'x':326.13, 'y': 300.70}}
		
		# calculate and add the coordinates of the nods as well
		coord_dict_det = generate_nods(coord_dict_det, verbose=verbose)
		
		
	else:
		# For slitless, we start with only the centre coordinate
		coord_dict_det = {'c': {'x': 38.5, 'y': 829.0}}
	
	if (frame == 'tel'):
		coord_dict = coord_dict_det.copy()
		for loc in coord_dict_det.keys():
			x,y = mt.xytov2v3(coord_dict_det[loc]['x']-1.,coord_dict_det[loc]['y']-1., 'F770W')
			if verbose:
				print('{0}, {1}'.format(x, y))
			coord_dict[loc] = {'x': x, 'y': y}
	
	elif (frame == 'idl'):
		# to go to ideal coordinatines, we need to go to v2v3 first. so this is an additional step.
		coord_dict = coord_dict_det.copy()
		for loc in coord_dict_det.keys():
			xtel,ytel = mt.xytov2v3(coord_dict_det[loc]['x']-1.,coord_dict_det[loc]['y']-1., 'F770W')
			x,y = mt.v2v3toIdeal(xtel, ytel, ap.AperName)
			coord_dict[loc] = {'x': x, 'y': y}
			
	
	else:
		coord_dict = coord_dict_det.copy()
	
	if plot:
		p = plot_pattern(slit=coord_dict, patt=None)
	
	
	return coord_dict

#----------------------------------------------------------------------------
def read_aperture(mode='slit', verbose=True):
	
	'''
	Function that loads and returns the right SIAF aperture for LRS as specified in the 'mode' parameter: 'slit' or 'slitless'

	'''
	if verbose:
		print("Pysiaf uses PRD version {}".format(pysiaf.JWST_PRD_VERSION))
	
	# check that the mode is a valid option
	assert (mode in ['slit', 'slitless']), "Mode not supported. Please use 'slit' or 'slitless'."
	
	# Read in the SIAF file using PySiaf
	instrument = 'MIRI'
	siaf = pysiaf.Siaf(instrument)
	
	# Load in the aperture corresponding to the specified mode:
	if (mode == 'slit'):
		ap = siaf['MIRIM_SLIT']
		
	elif (mode == 'slitless'):
		ap = siaf['MIRIM_SLITLESSPRISM']
		
	else:
		raise IOError('Mode not supported!')
		
	return ap
#----------------------------------------------------------------------------	
def plot_pattern(slit=None, patt=None):
	
	'''
	Function to plot the LRS slit and/or a dither pattern
	
	Parameters
	----------
	slit: a dictionary with slit corner coordinates, in any frame. [optional]
	patt: a Table object of size (N, 2) [optional]
	
	slit and patt can't both be None at the same time. 
	* If slit == None and the pattern is defined, the function will assume
	slitless operation. 
	* If patt == None and slit is defined, the function will plot the outline of the slit with no pattern
	* If patt and slit are defined, the function will assuem slit operation, and plot both the slit and the pattern
	
	
	
	Output
	------
	As described above
	
	'''
	
	fig, ax = plt.subplots(figsize=[12,4])
    
	if (slit is None) & (patt is None):
		raise IOError("Slit and Pattern can't both be None!")
	
	# pull out the corner coordinates from the slit dictionary and format for plotting
	if slit is not None:
		print("Slit coordinates are specified. Assuming slit operation.")
		cornersx = np.array((slit['ll']['x'][0], slit['ul']['x'][0], slit['ur']['x'][0], slit['lr']['x'][0]))
		cornersy = np.array((slit['ll']['y'][0], slit['ul']['y'][0], slit['ur']['y'][0], slit['lr']['y'][0]))
		corners = np.stack((cornersx, cornersy), axis=-1)
    
		# create the matplotlib patch for plotting
		slit_rect = Polygon(corners, closed=True, lw=2., color='g', fill=False, label='slit edge')
		
		ax.scatter(cornersx, cornersy, marker='+', color='g')
		ax.scatter(slit['c']['x'][0], slit['c']['y'][0], marker='o', color='r', facecolor=None, label='slit centre')
		ax.add_patch(slit_rect)
    
	if patt is not None:
		# create an index vector for colormapping the order of the points
		indx = np.arange(len(patt))
		pts = ax.scatter(pattern['xidl'], pattern['yidl'], marker='x', c=indx, cmap='plasma', label='pointings')
		cbar = fig.colorbar(pts, ax=ax)
		cbar.set_label('pointing index')
	else: 
		print("No pattern provided, plotting slit only.")
	
	
	ax.grid(color='k', linestyle='-', linewidth=0.5, alpha=0.5)
	ax.set_xlabel('arcsec')
	ax.set_ylabel('arcsec')
	ax.legend()
		

#	if save:
#	   plt.savefig(save) 
    
	fig.show()
    
	return fig
	
#----------------------------------------------------------------------------	

def generate_nods(slit, verbose=False):
	
	'''
	Function that takes a dictionary with slit coordinates as input and returns the locations of the nods, at 30 and 70% along the slit.
	
	Parameters:
	-----------
	- slit (dict): a dictionary with slit corner coordinates
	
	Output:
	-------
	- outslit (dict): a new dictionary; a copy of the input dictionary 'slit' with 2 additional coordinate sets appended
	'''
	
	
	outslit = slit.copy()
	
	# calculate the coordinates of the midpoints along the slit's edges:

	midl_x = np.mean([slit['ul']['x'], slit['ll']['x']])
	midl_y = np.mean([slit['ul']['y'], slit['ll']['y']])
	midl = [midl_x, midl_y]
	
	midr_x = np.mean([slit['ur']['x'], slit['lr']['x']])
	midr_y = np.mean([slit['ur']['y'], slit['lr']['y']])
	midr = [midr_x, midr_y]

	# now create a grid of 11 points between these points (each represents 10% along the slit) and pick the 4th and 8th in (x,y):
	xpts = np.linspace(midl[0], midr[0], num=11)
	ypts = np.linspace(midl[1], midr[1], num=11)
	slit_grid = np.stack((np.array(xpts), np.array(ypts)), axis=-1)
	
	outslit.update({'nod1': {'x': slit_grid[3,0], 'y': slit_grid[3,1]}, 'nod2': {'x': slit_grid[7,0], 'y': slit_grid[7 ,1]}})
	if verbose:
		print(outslit)
	return outslit
	
	