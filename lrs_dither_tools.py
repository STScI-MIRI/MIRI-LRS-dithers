import numpy as np
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import os
import glob
from astropy.table import Table, Column
import datetime
import sys


import pysiaf
import miricoord

plt.close('all')


class LRSPattern(object):
	
	def __init__(self, t, mode=None, frame=None, car=None, cal=None, notes=None):
		self.npts = len(t)
		self.mode = mode
		self.CAR = car
		self.CAL = cal
		self.notes = notes
		self.frame = frame
		


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
	
	ap = read_aperture(mode=mode)
	
	# the slit corner coordinates are hard coded unfortunately (in telescope frame):
	if (mode == 'slit'):
			
		coord_dict = {'ll': {'x': -411.99, 'y': -401.14},
	    	'ul': {'x': -411.95, 'y': -400.62}, 
	       	'ur': {'x': -416.68, 'y': -400.24},
	       	'lr': {'x': -416.72, 'y': -400.76}}
	
	else:
		
		coord_dict = {}
	
	coord_dict.update({'c': {'x': ap.V2Ref, 'y': ap.V3Ref}})
	
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
