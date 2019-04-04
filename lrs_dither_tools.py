import numpy as np
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
from os import access, R_OK
import glob
from astropy.table import Table, Column
import astropy.io.ascii as ascii
import datetime
import sys


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
			* 'det' for detector coordinates, in pixels
			* 'idl' for ideal coordinates, in arcsec. NOTE: slit and slitless mode have distinct ideal coordinate systems
	- car: the commissioning activity request(s) this pattern is used for. this should be specified in the format 'CAR-xxx'
	- cal: the calibration activity request(s) this pattern is used for. this should be specified in the format 'CAL-xxx'
	- ref: the reference point for this pattern. can be a single point or a list. choose from 'c', 'nod1', 'nod2'.
	- notes: any comments in a string format. PLEASE USE to add relevant information to the pattern!
	
	
	Methods:
	--------
	TO DO! THESE ARE IDEAS!
	- convert coordinates: convert between frames
	- plot: visualize the pattern
	- regenerate: regenerate in case of new SIAF update.
	- save: save to ascii file
	
	'''
	
	def __init__(self, table=None, file=None, mode=None, frame=None, car=None, cal=None, ref=None, notes=None):
		
		# check that either a table or a file are provided
		assert (table is not None) or (file is not None), "You must provide either a table or a pattern file to initialize an LRSPattern instance"
		
		if (file is not None):
			
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
			
					
				
		else:
			self.patt = t
			self.npts = len(t)
			self.mode = mode
			self.car = car
			self.cal = cal
			self.notes = notes
			self.frame = frame
			self.reps = reps
		


#----------------------------------------------------------------------------	
def lrs_gencoords(mode='slit', frame='tel', plot=False):
	
	'''
	Function that returns a dictionary with the key coordinates for the MIRI LRS. 
	
	Parameters:
	-----------
	- mode: 'slit' (default) or 'slitless'. this will determine which set of coordinates is returned
	- frame: the reference frame in which the coordinates are specified. Options: 'det' (detector [pixels]), 'idl' (Ideal [arcsec]), 'tel' (v2v3 [arcsec])
	
	
	Returns:
	-------
	- coord: a dictionary with the relevant coordinates. For SLIT, this will contain the slit centre and slit corner locations in the requested frame and units (see options above). For SLITLESS, this contains only the nominal pointing position.
	
	'''
	
	# checks that the parameters provided are valid
	assert (mode in ['slit', 'slitless']), "Mode not recognised. Options are: 'slit', 'slitless'. "
	assert (frame in ['tel', 'idl', 'det']), "Coordinate frame not recognised. Options are: 'det' (detector), 'idl' (ideal), 'tel' (telescope/v2v3)."
	
	# the read_aperture function (below) converts 'slit' or 'slitless' to the appropriate SIAF aperture name.
	ap = read_aperture(mode=mode)
	
	# the slit corner coordinates are hard coded unfortunately (in telescope frame):
	if (mode == 'slit'):
		# populate the dictionary with the slit corner coordinates
		coord_dict_tel = {'ll': {'x': -411.99, 'y': -401.14},
	    	'ul': {'x': -411.95, 'y': -400.62}, 
	       	'ur': {'x': -416.68, 'y': -400.24},
	       	'lr': {'x': -416.72, 'y': -400.76}}
	
	else:
		# we don't need teh slit corner coordinates for slitless mode
		coord_dict_tel = {}
	
	coord_dict_tel.update({'c': {'x': ap.V2Ref, 'y': ap.V3Ref}})
	
	if (frame == 'idl'):
		# create a copy of the telescope frame dictionary, and loop through the locations to convert the coordinates to ideal frame using the miricoord tools
		coord_dict = coord_dict_tel.copy()
		for loc in coord_dict_tel.keys():
			x,y = mt.v2v3toIdeal(coord_dict_tel[loc]['x'],coord_dict_tel[loc]['y'], ap.AperName)
			coord_dict[loc] = {'x': x, 'y': y}
			
	elif (frame == 'det'):
		# create a copy of the telescope frame dictionary, and loop through the locations to convert the coordinates to detector pixels using the miricoord tools
		#TO DO!!
		coord_dict = coord_dict_tel.copy()
		print("Sorry detector coordinates are not yet supported, returning v2v3 instead")
	
	else:
		coord_dict = coord_dict_tel.copy()
	
	
	print(coord_dict)
	
	
	if plot:
		p = plot_pattern(slit=coord_dict, patt=None)
	
	
	return coord_dict

#----------------------------------------------------------------------------
def read_aperture(mode='slit'):
	
	'''
	Function that loads and returns the right SIAF aperture for LRS as specified in the 'mode' parameter: 'slit' or 'slitless'

	'''
	
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
		cornersx = np.array((slit['ll']['x'], slit['ul']['x'], slit['ur']['x'], slit['lr']['x']))
		cornersy = np.array((slit['ll']['y'], slit['ul']['y'], slit['ur']['y'], slit['lr']['y']))
		corners = np.stack((cornersx, cornersy), axis=-1)
		print(corners)
    
		# create the matplotlib patch for plotting
		slit_rect = Polygon(corners, closed=True, lw=2., color='g', fill=False, label='slit edge')
		
		ax.scatter(cornersx, cornersy, marker='+', color='g')
		ax.scatter(slit['c']['x'], slit['c']['y'], marker='o', color='r', facecolor=None, label='slit centre')
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
