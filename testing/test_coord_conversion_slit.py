from LRSPattern import LRSPattern
from lrs_dither_tools import lrs_gencoords
import pysiaf

# THIS SCRIPT WILL TEST THE COORDINATE CONVERSIONS FOR SLIT 


f = '../test-patterns/test_slit_det.txt'

xx = lrs_gencoords(mode='slit', frame='det')

# First convert to DETECTOR ABSOLUTE and check that this matches the output from lrs_gencoords function (which pulls from the siaf)

print('******************************************************')
print('* TESTING: DETECTOR RELATIVE TO ABSOLUTE COORDINATES *')
print('******************************************************')
print('First method: the to_absolute() method')
pp = LRSPattern(file=f)
print('Input pattern in frame {}'.format(pp.frame))
print(pp.patt)
print(' ')

pp.to_absolute()
print('Output pattern in frame {}'.format(pp.frame))
print(pp.patt)
print(' ')

pat_ref = pp.patt.copy()


print('Second method: the convert_coords() method')
pp = LRSPattern(file=f)
print('Input pattern in frame {}'.format(pp.frame))
print(pp.patt)
print(' ')

pp.to_coordinates(new_frame='det')
print('Output pattern in frame {}'.format(pp.frame))
print(pp.patt)
print(' ')

# compare with the siaf coords
coord_comp = lrs_gencoords(mode='slit', frame='det')
print(coord_comp['c'])
print(' ')



# Second test: detector relative to TELESCOPE. compare against the siaf entries.
print('******************************************************')
print('* TESTING: DETECTOR RELATIVE TO TELESCOPE (V2V3)     *')
print('******************************************************')
pp = LRSPattern(file=f)
print('Input pattern in frame {}\n\n'.format(pp.frame))

pp.to_coordinates(new_frame='tel')
print('Output pattern in frame {}\n\n'.format(pp.frame))
print(pp.patt)
print(' ')

siaf = pysiaf.Siaf('MIRI')
ap = siaf['MIRIM_SLIT']
print('Comparison with SIAF value for slit centre:')
ref_tel = ap.reference_point(to_frame='tel')
print(ref_tel)
print ('DIFFERENCE IN SLIT CENTRE POS = {0:.5f} in X, {1:.5f} in Y'.format(pp.patt['x'][4]-ref_tel[0], pp.patt['y'][4]-ref_tel[1]))
print('\n\n')

print('Converting back to detector')
pp.to_coordinates(new_frame='det')

for i in range(len(pp.patt['x'])):
	print('Position {0}: delta_x = {1:.5f} px     delta_y = {2:.5f} px'.format(pp.patt['Pointing'][i], pp.patt['x'][i]-pat_ref['x'][i], pp.patt['y'][i]-pat_ref['y'][i]))




print('******************************************************')
print('* TESTING: DETECTOR RELATIVE TO IDEAL                *')
print('******************************************************')
pp = LRSPattern(file=f)

print('Input pattern in frame {}\n\n'.format(pp.frame))

pp.to_coordinates(new_frame='idl')
print('Output pattern in frame {}\n\n'.format(pp.frame))
print(pp.patt)
print(' ')

print('Comparison with SIAF value for slit centre:')
ref_idl = ap.reference_point(to_frame='idl')
print(ref_idl)

print ('DIFFERENCE IN SLIT CENTRE POS = {0:.5f} in X, {1:.5f} in Y'.format(pp.patt['x'][4]-ref_idl[0], pp.patt['y'][4]-ref_idl[1]))

print('Converting back to detector')
pp.to_coordinates(new_frame='det')

for i in range(len(pp.patt['x'])):
	print('Position {0}: delta_x = {1:.5f} px     delta_y = {2:.5f} px'.format(pp.patt['Pointing'][i], pp.patt['x'][i]-pat_ref['x'][i], pp.patt['y'][i]-pat_ref['y'][i]))






