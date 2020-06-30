import numpy as np
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import os
from os import access, R_OK
import glob
from astropy.table import Table, Column, vstack
import astropy.io.ascii as ascii
import datetime
import sys
import pdb


import pysiaf
from lrs_dither_tools import lrs_gencoords
import miricoord
import miricoord.imager.mirim_tools as mt

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
			self.frame = header['val'][header['key']=='Frame'][0]
			if ('/' in file):
				self.name = (file.split('/')[-1]).split('.')[0]
			else:
				self.name = file.split('.')[0]
			
			
			# first check if there are multiple reference positions - NOTE this is only relevant for patterns in frame det-rel
			ref_tmp = (self.ref[0]).split(',')
			# strip out whitespaces in case there are any:
			refs = [r.strip() for r in ref_tmp]
			self.nref = len(refs)
			
			# if there is only 1 reference position, do nothing. if there are N>1 positions, copy the x and y columns N-1 times to make space for the extra pointings. NOTE: if there ar emultiple reference positions, the coordinate frame shouls be relative not absolute - otherwise that wouldn't make sense.
			
			if (self.nref > 1):
				assert ('rel' in self.frame), "Cannot have multiple references for the pattern with an absolute coordinate frame"
				i = 1
				while (i <= self.nref-1):
					self.patt.add_column(t['x'].copy(), name='x{0}'.format(i+1))
					self.patt.add_column(t['y'].copy(), name='y{0}'.format(i+1))	
					i += 1		
					
					# for each additional reference point, create a copy of the table that we will then stack with the original vertically (vstack)
					#new_patt = self.patt[:self.npts].copy()
					#new_patt['Pointing'] += self.patt['Pointing'][:self.npts]
					#self.patt = vstack([self.patt, new_patt])
					#i += 1
				#self.patt['Pointing'] = np.arange(1, len(self.patt)+1)
					
					#self.patt.add_column(t['x'].copy(), name='x{0}'.format(i+1))
					#self.patt.add_column(t['y'].copy(), name='y{0}'.format(i+1))	
					#i += 1		
				
				
			
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
				assert ('rel' in self.frame), "Cannot have multiple references for the pattern with an absolute coordinate frame"
				i = 1
				while (i <= self.nref-1):
					#new_patt = self.patt[:self.npts].copy()
					#new_patt['Pointing'] += self.patt['Pointing'][:self.npts]
					self.patt.add_column(t['x'].copy(), name='x{0}'.format(i+1))
					self.patt.add_column(t['y'].copy(), name='y{0}'.format(i+1))	
					#self.patt = vstack([self.patt, new_patt])
					i += 1	
				#self.patt['Pointing'] = np.arange(1, len(self.patt)+1)
					
					
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
		assert ('rel' in self.frame), "The input pattern is not in relative coordinates!"
		
		print('Converting pattern to absolute (1-INDEXED) coordinates.....')
		# if the input coordinates are relative, load the coordinates so we can do the translation
		coords = lrs_gencoords(mode=self.mode, frame=self.frame[:3], verbose=False)
		
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
				colx = 'x{0}'.format(i+1)
				coly = 'y{0}'.format(i+1)
				self.patt[colx] += coords[r]['x']
				self.patt[coly] += coords[r]['y']
		
		# now update the frame attribute from relative to absolute:
		self.frame = '{0}-abs'.format(self.frame[:3])
			
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
		#coords = lrs_gencoords(mode=self.mode, frame=self.frame[:3])
		
		# figure out if there are multiple refrence positions
		ref_tmp = (self.ref[0]).split(',')
		# strip out whitespaces in case there are any:
		refs = [r.strip() for r in ref_tmp]
		nref = len(refs)
		#pdb.set_trace()
		
		
		
		fig, ax = plt.subplots(figsize=[12,4])
		
		if (self.frame == 'det-abs') or (self.frame == 'det-rel'):
			units = 'pixels'
			pltframe = 'det'
		else:
			units = 'arcsec'
			pltframe = self.frame
		
		#print(pltframe)
		#pdb.set_trace()
		
		# check whether the pattern is for slit or slitless
		if (self.mode == 'slit'):
			coords = lrs_gencoords(mode='slit', frame=pltframe)
			print(coords)
			
			# if the coordinates are in relative pixel coordinates, translate the pattern 
			if (self.frame == 'det-rel'):
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
					colx = 'x{0}'.format(i+1)
					coly = 'y{0}'.format(i+1)
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
	
	def to_coordinates(self, new_frame=None):
		
		'''
		This function will convert the pattern coordinates to a new frame, specified in the frame keyword.
		
		Parameters:
		-----------
		- frame:	the output coordinate frame. can be 'tel', 'idl', or 'det'. the frame 'det' will be interpreted as 'det-abs'
		
		Output:
		--------
		- the LRS pattern, with pointings and frame attribute updated
		
		
		'''
		
		assert (new_frame != self.frame), "Output frame is the same as the current frame!"
		assert (new_frame in ['det', 'tel', 'idl']), "Frame not recognised. Please specify 'det', 'tel' or 'idl'"
		
		
		colpt = self.patt['Pointing'].copy()
		ncols = len(self.patt.colnames)
		
		if (new_frame == 'tel'):
			print('Converting coordinates to Telescope (v2v3) frame')
			
			# NOTE: convert the detector coordinates to ZERO-INDEXED coordinates before conversion. 
			#pdb.set_trace()
			
			# CASE 1.1: FROM DETECTOR FRAME TO TELESCOPE FRAME
			if ('det' in self.frame):
				
				if ('rel' in self.frame):
					self.to_absolute()
				
				if (self.nref == 1):
					
					new_patt = mt.xytov2v3(self.patt['x']-1., self.patt['y']-1., 'F770W')
					colx = Column(new_patt[0], name=self.patt.colnames[1])
					coly = Column(new_patt[1], name=self.patt.colnames[2])
					new_patt_table = Table([colpt, colx, coly])
					
					# update the Pattern attributes:
					self.patt = new_patt_table
					self.frame = new_frame
				
				else:
					new_patt_table = Table([colpt])
					for i in range(self.nref):
						col_keys = [self.patt.keys()[(i*2)+1], self.patt.keys()[(i*2)+2]]
						new_patt = mt.xytov2v3(self.patt[col_keys[0]]-1., self.patt[col_keys[1]]-1., 'F770W')
						colx = Column(new_patt[0], name=col_keys[0])
						coly = Column(new_patt[1], name=col_keys[1])
						new_patt_table.add_columns([colx, coly])
					
					self.patt = new_patt_table
					self.frame = new_frame
						
			# CASE 1.2: FROM IDEAL FRAME TO TELESCOPE FRAME
			elif ('idl' in self.frame):
				
				if (self.nref == 1):
					
					if (self.mode[0] == 'slit'):
						new_patt = mt.Idealtov2v3(self.patt['x'], self.patt['y'], 'MIRIM_SLIT')
					elif (self.mode[0] == 'slitless'):
						new_patt = mt.Idealtov2v3(self.patt['x'], self.patt['y'], 'MIRIM_SLITLESSPRISM')
					colx = Column(new_patt[0], name=self.patt.colnames[1])
					coly = Column(new_patt[1], name=self.patt.colnames[2])
					new_patt_table = Table([colpt, colx, coly])
					
					# update the Pattern attributes:
					self.patt = new_patt_table
					self.frame = new_frame
				
				else: 
					
					new_patt_table = Table([colpt])
					for i in range(self.nref):
						col_keys = [self.patt.keys()[(i*2)+1], self.patt.keys()[(i*2)+2]]
						if (self.mode[0] == 'slit'):
							new_patt = mt.Idealtov2v3(self.patt[col_keys[0]], self.patt[col_keys[1]], 'MIRIM_SLIT')
						elif (self.mode[0] == 'slitless'):
							new_patt = mt.Idealtov2v3(self.patt[col_keys[0]], self.patt[col_keys[1]], 'MIRIM_SLITLESSPRISM')
						colx = Column(new_patt[0], name=col_keys[0])
						coly = Column(new_patt[1], name=col_keys[1])
						new_patt_table.add_columns([colx, coly])
					self.patt = new_patt_table
					self.frame = new_frame
					
			
			# CASE 1.3: FROM TELESCOPE FRAME TO TELESCOPE FRAME: DO NOTHING!		
			elif ('tel' in self.frame):
				
				print('Pattern already in frame {}'.format(new_frame))
				
				
					
		# CASE 2: TO IDEAL FRAME			
		elif (new_frame == 'idl'):
			print('Converting coordinates to Ideal frame')
			
			#CASE 2.1: FROM DETECTOR FRAME TO IDEAL FRAME - THIS CONVERSION FIRST GOES TO V2V3 THEN TO IDEAL.
			if ('det' in self.frame):
				
				if ('rel' in self.frame):
					self.to_absolute()
				
				if (self.nref == 1):
					
					tel_patt = mt.xytov2v3(self.patt['x']-1., self.patt['y']-1., 'F770W')
					if (self.mode[0] == 'slit'):
						new_patt = mt.v2v3toIdeal(tel_patt[0], tel_patt[1], 'MIRIM_SLIT')
					elif (self.mode[0] == 'slitless'):
						new_patt = mt.v2v3toIdeal(tel_patt[0], tel_patt[1], 'MIRIM_SLITLESSPRISM')
					colx = Column(new_patt[0], name=self.patt.colnames[1])
					coly = Column(new_patt[1], name=self.patt.colnames[2])
					new_patt_table = Table([colpt, colx, coly])
					
					# update the Pattern attributes:
					self.patt = new_patt_table
					self.frame = new_frame
				
				else:
					
					#pdb.set_trace()
					
					new_patt_table = Table([colpt])
					for i in range(self.nref):
						col_keys = [self.patt.keys()[(i*2)+1], self.patt.keys()[(i*2)+2]]
						tel_patt = mt.xytov2v3(self.patt[col_keys[0]]-1., self.patt[col_keys[1]]-1., 'F770W')
						if (self.mode[0] == 'slit'):
							new_patt = mt.v2v3toIdeal(tel_patt[0], tel_patt[1], 'MIRIM_SLIT')
						elif (self.mode[0] == 'slitless'):
							new_patt = mt.v2v3toIdeal(tel_patt[0], tel_patt[1], 'MIRIM_SLITLESSPRISM')
						colx = Column(new_patt[0], name=col_keys[0])
						coly = Column(new_patt[1], name=col_keys[1])
						new_patt_table.add_columns([colx, coly])
					# Update the pattern attributes
					self.patt = new_patt_table
					self.frame = new_frame
					
					
					
					
					

			# CASE 2.2: FROM TELESCOPE FRAME TO IDEAL FRAME
			elif ('tel' in self.frame):
				
				if (self.nref == 1):
					
					if (self.mode[0] == 'slit'):
						new_patt = mt.v2v3toIdeal(self.patt['x'], self.patt['y'], 'MIRIM_SLIT')
					elif (self.mode[0] == 'slitless'):
						new_patt = mt.v2v3toIdeal(self.patt['x'], self.patt['y'], 'MIRIM_SLITLESSPRISM') 
					colx = Column(new_patt[0], name=self.patt.colnames[1])
					coly = Column(new_patt[1], name=self.patt.colnames[2])
					new_patt_table = Table([colpt, colx, coly])
					# update the Pattern attributes:
					self.patt = new_patt_table
					self.frame = new_frame
				
				else:
										
					new_patt_table = Table([colpt])
					for i in range(self.nref):
						col_keys = [self.patt.keys()[(i*2)+1], self.patt.keys()[(i*2)+2]]
						if (self.mode[0] == 'slit'):
							new_patt = mt.v2v3toIdeal(self.patt[col_keys[0]], self.patt[col_keys[1]], 'MIRIM_SLIT')
						elif (self.mode[0] == 'slitless'):
							new_patt = mt.v2v3toIdeal(self.patt[col_keys[0]], self.patt[col_keys[1]], 'MIRIM_SLITLESSPRISM')
						colx = Column(new_patt[0], name=col_keys[0])
						coly = Column(new_patt[1], name=col_keys[1])
						new_patt_table.add_columns([colx, coly])
					# Update the pattern attributes
					self.patt = new_patt_table
					self.frame = new_frame
					
			# CASE 2.3: PATTERN ALREADY IN IDEAL FRAME. DO NOTHING!
			elif ('idl' in self.frame):
				print('Pattern already in frame {}'.format(new_frame))
			
			
		#CASE 3: TO DETECTOR FRAME
		elif (new_frame == 'det'):
			print('Converting coordinates to Detector frame')
			# NOTE the tel -> detector conversion returns ZERO-INDEXED coordinates, so need to add 1 to put into SIAF frame
			
			
			# CASE 3.1: ALREADY IN DETECTOR FRAME, BUT RELATIVE COORDINATES. THIS IS THE SAME AS RUNNING THE TO_ABSOLUTE() METHOD
			if ('rel' in self.frame):
				self.to_absolute()
			
			# CASE 3.2: ALREADY IN ABSOLUTE DETECTOR FRAME. DO NOTHING.
			elif ('abs' in self.frame):
				
				print('Pattern already in frame {}'.format(new_frame))
				
				
			# CASE 3.3: FROM TELESCOPE FRAME TO ABSOLUTE DETECTOR FRAME	
			elif ('tel' in self.frame):
					
				if (self.nref == 1):
						
					new_patt = mt.v2v3toxy(self.patt['x'], self.patt['y'], 'F770W')
					colx = Column(new_patt[0]+1., name=self.patt.colnames[1])
					coly = Column(new_patt[1]+1., name=self.patt.colnames[2])
					new_patt_table = Table([colpt, colx, coly])
				

				else:
					new_patt_table = Table([colpt])
					for i in range(self.nref):
						col_keys = [self.patt.keys()[(i*2)+1], self.patt.keys()[(i*2)+2]]
						new_patt = mt.v2v3toxy(self.patt[col_keys[0]], self.patt[col_keys[1]], 'F770W')
						colx = Column(new_patt[0]+1., name=col_keys[0])
						coly = Column(new_patt[1]+1., name=col_keys[1])
						new_patt_table.add_columns([colx, coly])
				
				# Update the pattern attributes
				self.patt = new_patt_table
				self.frame = new_frame
					
			
			# CASE 3.4: FROM IDEAL FRAME TO ABSOLUTE DETECTOR FRAME. 
			elif ('idl' in self.frame):
				if (self.nref == 1):	
					if (self.mode[0] == 'slit'):
						tel_patt = mt.Idealtov2v3(self.patt['x'], self.patt['y'], 'MIRIM_SLIT')
					elif (self.mode[0] == 'slitless'):
						tel_patt = mt.Idealtov2v3(self.patt['x'], self.patt['y'], 'MIRIM_SLITLESSPRISM')
					new_patt = mt.v2v3toxy(tel_patt[0], tel_patt[1], 'F770W')
					colx = Column(new_patt[0]+1., name=self.patt.colnames[1])
					coly = Column(new_patt[1]+1., name=self.patt.colnames[2])
					new_patt_table = Table([colpt, colx, coly])
				
				else:
					new_patt_table = Table([colpt])
					for i in range(self.nref):
						col_keys = [self.patt.keys()[(i*2)+1], self.patt.keys()[(i*2)+2]]
						if (self.mode[0] == 'slit'):
							tel_patt = mt.Idealtov2v3(self.patt[col_keys[0]], self.patt[col_keys[1]], 'MIRIM_SLIT')
						elif (self.mode[0] == 'slitless'):
							tel_patt = mt.Idealtov2v3(self.patt[col_keys[0]], self.patt[col_keys[1]], 'MIRIM_SLITLESSPRISM')
						new_patt = mt.v2v3toxy(tel_patt[0], tel_patt[1], 'F770W')
						colx = Column(new_patt[0]+1., name=col_keys[0])
						coly = Column(new_patt[1]+1., name=col_keys[1])
						new_patt_table.add_columns([colx, coly])
				
				
				
				# update the Pattern attributes:
				self.patt = new_patt_table
				self.frame = new_frame
					
				
					
					
		return
		
		

#----------------------------------------------------------------------------	

	def write(self, frame=None, outfile=None):
		
		'''
		Function that will write out a single dither pattern to file, in any supported coordinate frame. 
		
		Parameters:
		-----------
		- frame:	the output reference frame (supported: 'det-abs', 'idl', 'tel')
		- outfile:	output filename (optional)
		
		Output:
		-------
		outfile:	output filename. if no name is provided, a filename will be created from the mode, format and creation date.



		'''
		
		today = datetime.date.today().isoformat()
		
		# create an iterable list for the reference positions. works fine for a single reference location as well:
		refs_tmp = (self.ref[0]).split(',')
		# this will just remove any whitespace from the reference locations:
		refs = [r.strip() for r in refs_tmp]
		
		# define a filename either based on the outfile parameter or based on the pattern name, today's date, frame and reference location. i.e. if there's > 1 ref, multiple files will be created
		if outfile:
			out = ['output-patterns/{0}_{1}.txt'.format(outfile, rr) for rr in refs]

		else: 
			out = ['output-patterns/{0}_{1}_{2}_{3}.txt'.format(self.name, rr, frame, today) for rr in refs]
			
		
		# convert to the right coordinate frame if needed:
		if (self.frame != frame):
			self.to_coordinates(new_frame=frame)
			
			
			
		
		# if there's only one reference position, take the pattern, convert coordinates and write out to file:
		if (len(refs) == 1):
			new_patt = Table([self.patt.columns[0], self.patt.columns[1], self.patt.columns[2]], names=['Pointing', 'X', 'Y'], dtype=['i4', 'f8', 'f8'])
			ascii.write(new_patt, output=out[0], format='basic', delimiter='\t', formats={'X': '%.3f', 'Y': '%.3f'}, overwrite=True)
			
			#new_patt.write(out[0], format='ascii.no_header', formats={'X': '%:.3f', 'Y': '%:.3f'})
		
		
		# if there's multiple reference positions, we have multiple output files. first let's check that the number of refrence points corresponds to the number of output file, and that the number of ref points matches the number of columns in the pattern
		else:
			assert len(refs) == len(out), "Number of reference points doesn't match output files"
			assert len(refs) == (len(self.patt.columns)-1)/2, "Number of reference points doesn't match pattern shape"
			
			# if checks are okay, loop through output files and write out the appropriate columns from the pattern
			
			for i, oo in enumerate(out):
				# create a new output table based on the pattern. we always wabt column[0] (Pointing number), and then 2 columns for x and y
				ref_patt = Table([self.patt.columns[0], self.patt.columns[(i*2)+1], self.patt.columns[(i*2)+2]], names=['Pointing', 'X', 'Y'], dtype=['i4', 'f8', 'f8'])
				ascii.write(ref_patt, output=oo, format='basic', delimiter='\t', formats={'X': '%.3f', 'Y': '%.3f'}, overwrite=True)
				
				
			
		
		return
		
######################################################################################

def write_prd(input_dir='', output_dir=''):
	
	'''This function will take all the patterns in a given directory and write them to a file suitable for PRD ingestion. 
	
	Parameters
	----------
	- input_dir: 	a directory with pattern files
	- output_dir:	where to place the output file
	
	
	
	'''	
	
	today = datetime.date.today().isoformat()
	
	if (output_dir==''):
		
		output_dir = 'output-patterns/'
	
	if ('slitless' in input_dir):
		out_tmp = 'MiriLrsDithers_slitless_{}.txt'.format(today)
	else:
		out_tmp = 'MiriLrsDithers_slit_{}.txt'.format(today)
		
	outfile = os.path.join(output_dir, out_tmp)
	
	pfiles = glob.glob(input_dir+'*.txt')

	for i, pf in enumerate(pfiles):
		pp = LRSPattern(file=pf)
		pp.to_coordinates(new_frame='idl')
		# extract the name of the pattern. if there's only 1 reference, then this is the name from the metadata with the _ replaced by whitespace
		if (pp.nref == 1):
			pname = pp.name.replace('_', ' ')
			
			# for the first file we want to write the output file, the rest is appended
			if (i==0):
				print(pname,file=open(outfile,"w+"))
				for j in range(pp.npts):
					print('{0:<3}{1:>8.6f}         {2:>8.6f}'.format(pp.patt['Pointing'][j], pp.patt['x'][j], pp.patt['y'][j]), file=open(outfile, "a"))
		
			else:
				print(pname,file=open(outfile,"a"))
				for j in range(pp.npts):
					print('{0:<3}{1:>8.6f}         {2:>8.6f}'.format(pp.patt['Pointing'][j], pp.patt['x'][j], pp.patt['y'][j]), file=open(outfile, "a"))

			print(' ', file=open(outfile, "a"))
			
		else:
			# create an iterable list for the reference positions. works fine for a single reference location as well:
			refs_tmp = (pp.ref[0]).split(',')
			# this will just remove any whitespace from the reference locations:
			refs = [r.strip() for r in refs_tmp]
			pname = [(pp.name+'_'+rr).replace('_', ' ') for rr in refs]
			
			for r in range(pp.nref):
				col_keys = [pp.patt.keys()[(r*2)+1], pp.patt.keys()[(r*2)+2]]
				
				if (i==0):
					print(pname[r],file=open(outfile,"w+"))
					for j in range(pp.npts):
						print('{0:<3}{1:>8.6f}         {2:>8.6f}'.format(pp.patt['Pointing'][j], pp.patt[col_keys[0]][j], pp.patt[col_keys[1]][j]), file=open(outfile, "a"))
				else:
					print(pname[r],file=open(outfile,"a"))
					for j in range(pp.npts):
						print('{0:<3}{1:>8.6f}         {2:>8.6f}'.format(pp.patt['Pointing'][j], pp.patt[col_keys[0]][j], pp.patt[col_keys[1]][j]), file=open(outfile, "a"))
				print(' ', file=open(outfile, "a"))
				
			
			# if there are multiple references, then split them out and print individually
			
		
		
	
		
	
	
	
	