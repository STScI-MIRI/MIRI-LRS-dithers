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
from lrs_dither_tools import lrs_gencoords
import miricoord
import miricoord.miricoord.imager.mirim_tools as mt

plt.close('all')


class LRSPattern(object):
	
	'''
	This code defines a class for an LRS dither pattern. It is initialised from a table of positions (this should be an astropy table), or from an ascii file that uses the standard pattern template (see LRS_pattern_template.txt). 
	
	Initialisation:
	---------------
	- file: a filename. If the file passed is not in the standard template, the initialization will fail as the metadata will likely not be found.
	- pattern: an Astropy table with position index, x and y location. if the class is initialised with a table, metadata must be provided with the call to LRSPattern. 
	
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
	
	def __init__(self, pattern=None, file=None, mode=None, frame=None, car=None, cal=None, name=None, ref=None, notes=None):
		
		# check that either a table or a file are provided
		assert (pattern is not None) or (file is not None), "You must provide either a table or a pattern file to initialize an LRSPattern instance"
		
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
			
			self.patt = pattern
			self.npts = len(pattern)
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
		
		
	
			

	def to_absolute(self):
		
		'''
		This function will change the coordinates from relative to absolute. This is ONLY applicable for patterns defined in the 'det-rel' frame; it will take the pattern, identify the reference coordinates for the given mode, and calculate the absolute detector pixels.
		
		Parameters:
		-----------
		none
		
		Returns:
		--------
		- the input pattern with updated positions, and the reference frame set to 'det-abs' in the metadata. if the out keyword was set, the original pattern will be the same and a new pattern is created
		
		
		
		'''
		
		# first check that the input pattern is defined in relative coordinates
		assert ('rel' in self.frame[0]), "The input pattern is not in relative coordinates!"
		
		print('Converting pattern to absolute coordinates.....')
		# if the input coordinates are relative, load the coordinates so we can do the translation
		coords = lrs_gencoords(mode=self.mode, frame=self.frame[0][:3], verbose=False)
		
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
		self.frame = ['{0}-abs'.format(self.frame[0][:3])]
			
		return
		
	
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
		#coords = lrs_gencoords(mode=self.mode, frame=self.frame[0][:3])
		
		# figure out if there are multiple refrence positions
		ref_tmp = (self.ref[0]).split(',')
		# strip out whitespaces in case there are any:
		refs = [r.strip() for r in ref_tmp]
		nref = len(refs)
		print(refs)
		#pdb.set_trace()
		
		
		
		fig, ax = plt.subplots(figsize=[12,4])
		
		if (self.frame[0] == 'det-abs') or (self.frame[0] == 'det-rel'):
			units = 'pixels'
			pltframe = 'det'
		else:
			units = 'arcsec'
			pltframe = self.frame[0]
		
		# check whether the pattern is for slit or slitless
		if (self.mode == 'slit'):
			coords = lrs_gencoords(mode='slit', frame=pltframe)
			
			# if the coordinates are in relative pixel coordinates, translate the pattern 
			if (self.frame[0] == 'det-rel'):
				print('NOTE: converting the pattern to absolute detector coordinates')
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
				self.to_absolute()
		
		if (nref == 1):
			
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
		
		else:
			# if there are multiple reference points, need to plot from multiple sets of columns
			for i, r in enumerate(refs):
				if (i==0):
					pts = ax.scatter(self.patt['x'], self.patt['y'], marker='x', c=self.patt['Pointing'], cmap='plasma', label='pointings ({0}, {1})'.format(self.npts, r))
				else:
					colx = 'x{0}'.format(i)
					coly = 'y{0}'.format(i)
					pts = ax.scatter(self.patt[colx], self.patt[coly], marker='x', c=self.patt['Pointing'], cmap='plasma', label='pointings ({0}, {1})'.format(self.npts, r))
					
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
		
		
	def run_checks(self):
		
		'''In this function we will run some sanity checks on the coordinates and reference frames provided, so that we don't waste time or create unrealistic patterns. DOESN'T CURRENTLY WORK.
		
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
			
			
		return



#----------------------------------------------------------------------------	