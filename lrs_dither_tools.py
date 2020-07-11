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
import miricoord.imager.mirim_tools as mt


def lrs_gencoords(mode='slit', frame='tel', plot=False, verbose=False):
	
	'''
	Function that returns a dictionary with the key coordinates for the MIRI LRS. 
	
	Parameters:
	-----------
	- mode: 'slit' (default) or 'slitless'. this will determine which set of coordinates is returned
	- frame: the reference frame in which the coordinates are specified. Options: 'det' (detector [pixels]), 'idl' (Ideal [arcsec]), 'tel' (v2v3 [arcsec])
	- verbose: prints some output to screen (boolean, default False)
	
	
	Returns:
	-------
	- coord: a dictionary with the relevant coordinates. For SLIT, this will contain the slit centre and slit corner locations in the requested frame and units (see options above). For SLITLESS, this contains only the nominal pointing position.
	
	'''
	
	# checks that the parameters provided are valid
	#print(frame)
	#pdb.set_trace()
	assert (mode in ['slit', 'slitless']), "Mode not recognised. Options are: 'slit', 'slitless'. "
	assert (frame in ['tel', 'idl', 'det']), "Coordinate frame not recognised. Options are: 'det' (detector), 'idl' (ideal), 'tel' (telescope/v2v3)."
	
	# the read_aperture function (below) converts 'slit' or 'slitless' to the appropriate SIAF aperture name.
	ap = read_aperture(mode=mode)
	
	# get the slit centre coordinates from the SIAF
	# if mode is slit then have to do this in v2v3 first, if slitless can get it directly in detector coordinates. remember to add 1!
	
	# the slit corner coordinates are hard coded unfortunately (in detector units):
	if (mode == 'slit'):
		
		xc, yc = ap.reference_point(to_frame='tel')
		xxc, yyc = mt.v2v3toxy(xc, yc, 'F770W')
		# populate the dictionary with the slit corner and centre coordinates, in PIXELS and 1-INDEXED.
		coord_dict_det = {'ll': {'x': 304.77, 'y': 298.38},
			    	'ul': {'x': 304.77, 'y': 303.03}, 
			       	'ur': {'x': 347.49, 'y': 303.03},
			       	'lr': {'x': 347.49, 'y': 298.38},
					'c': {'x':xxc[0]+1., 'y': yyc[0]+1.}}
		
		# calculate and add the coordinates of the nods as well
		coord_dict_det = generate_nods(coord_dict_det, verbose=verbose)
		
		print(coord_dict_det)
		
	else:
		# For slitless, we start with only the centre coordinate. DONT HAVE TO ADD 1 IF PULLING DIRECTLY FROM THE SIAF IN XY COORDS!!
		xc, yc = ap.reference_point(to_frame='det')
		coord_dict_det = {'c': {'x': xc, 'y': yc}}
	
	if (frame == 'tel'):
		coord_dict = coord_dict_det.copy()
		for loc in coord_dict_det.keys():
			x,y = mt.xytov2v3(coord_dict_det[loc]['x']-1.,coord_dict_det[loc]['y']-1., 'F770W')
			if verbose:
				print('{0}, {1}'.format(x, y))
			coord_dict[loc] = {'x': x[0], 'y': y[0]}
	
	elif (frame == 'idl'):
		# to go to ideal coordinatines, we need to go to v2v3 first. so this is an additional step.
		coord_dict = coord_dict_det.copy()
		for loc in coord_dict_det.keys():
			xtel,ytel = mt.xytov2v3(coord_dict_det[loc]['x']-1.,coord_dict_det[loc]['y']-1., 'F770W')
			x,y = mt.v2v3toIdeal(xtel, ytel, ap.AperName)
			coord_dict[loc] = {'x': x[0], 'y': y[0]}
			
	
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
	# (force the y coordinate to be the same as that of the slit centre in xy coordinates!)
	midl_x = np.mean([slit['ul']['x'], slit['ll']['x']])
	midl = [midl_x, slit['c']['y']]
	#midl_y = np.mean([slit['ul']['y'], slit['ll']['y']])
	#midl = [midl_x, midl_y]
	
	midr_x = np.mean([slit['ur']['x'], slit['lr']['x']])
	midr = [midl_x, slit['c']['y']]
	#midr_y = np.mean([slit['ur']['y'], slit['lr']['y']])
	#midr = [midr_x, midr_y]

	# now create a grid of 11 points between these points (each represents 10% along the slit) and pick the 4th and 8th in (x,y):
	xpts = np.linspace(midl[0], midr[0], num=11)
	ypts = np.linspace(midl[1], midr[1], num=11)
	slit_grid = np.stack((np.array(xpts), np.array(ypts)), axis=-1)
	
	outslit.update({'nod1': {'x': slit_grid[3,0], 'y': slit_grid[3,1]}, 'nod2': {'x': slit_grid[7,0], 'y': slit_grid[7 ,1]}})
	if verbose:
		print(outslit)
	return outslit

#----------------------------------------------------------------------------	
def generate_dither_file(mode=None, format=None, outfile=None, verbose=False):
    
    '''
    Function that will write out dither pattern files for selected purposes. Supported are: APT, MIRISim. 

    
    Parameters:
    -----------
    - mode:     'slit' or 'slitless'
    - format:   the tool the file will be read by. this will determine the output reference frame.
    - verbose:  Boolean, for extra output. default='False'
    - outfile:  output filename
    
    Output:
    -------
    outfile:    output filename. if no name is provided, a filename will be created from the mode, format and creation date.
    
    
    '''

    # take the mode, and read in all patterns in that pattern directory
    
    # for APT, convert all to Ideal coordinates. for MIRISim, absolute detector coordinates.
    
    # remember to split hem out if there are multiple reference locations
    
    # create a big file and print all patterns to file
    
    return outfile

#----------------------------------------------------------------------------	
