from LRSPattern import LRSPattern
from lrs_dither_tools import lrs_gencoords
import miricoord
import miricoord.imager.mirim_tools as mt
import pdb
import pysiaf

# THIS SCRIPT WILL TEST THE COORDINATE CONVERSIONS FOR A SLIT PATTERN WITH MULTIPLE REFERENCES


f = '../test-patterns/test_slit_det_multiref.txt'

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
print('SLIT CENTRE: {}'.format(coord_comp['c']))
print('NOD 1: {}'.format(coord_comp['nod1']))
print('NOD 2: {}'.format(coord_comp['nod2']))

print(' ')

#pdb.set_trace()

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
ref_c_tel = ap.reference_point(to_frame='tel')
ref_nod1_tel = mt.xytov2v3(xx['nod1']['x']-1., xx['nod1']['y']-1., 'F770W')
ref_nod2_tel = mt.xytov2v3(xx['nod2']['x']-1., xx['nod2']['y']-1., 'F770W')

print('DIFFERENCE IN SLIT CENTRE POS = {0:.5f} in X, {1:.5f} in Y'.format(pp.patt['x'][4]-ref_c_tel[0], pp.patt['y'][4]-ref_c_tel[1]))
print('DIFFERENCE IN NOD 1 POS = {0:.5f} in X, {1:.5f} in Y'.format(pp.patt['x2'][4]-ref_nod1_tel[0][0], pp.patt['y2'][4]-ref_nod1_tel[1][0]))
print('DIFFERENCE IN NOD 1 POS = {0:.5f} in X, {1:.5f} in Y'.format(pp.patt['x3'][4]-ref_nod2_tel[0][0], pp.patt['y3'][4]-ref_nod2_tel[1][0]))
print('\n\n')

print('Converting back to detector')
pp.to_coordinates(new_frame='det')


print('Ref pos 1 check:')
for i in range(len(pp.patt['x'])):
	print('Position {0}: delta_x = {1:.5f} px     delta_y = {2:.5f} px'.format(pp.patt['Pointing'][i], pp.patt['x'][i]-pat_ref['x'][i], pp.patt['y'][i]-pat_ref['y'][i]))

print('Ref pos 2 check:')
for i in range(len(pp.patt['x'])):
	print('Position {0}: delta_x = {1:.5f} px     delta_y = {2:.5f} px'.format(pp.patt['Pointing'][i], pp.patt['x2'][i]-pat_ref['x2'][i], pp.patt['y2'][i]-pat_ref['y2'][i]))
	
print('Ref pos 3 check:')
for i in range(len(pp.patt['x'])):
	print('Position {0}: delta_x = {1:.5f} px     delta_y = {2:.5f} px'.format(pp.patt['Pointing'][i], pp.patt['x3'][i]-pat_ref['x3'][i], pp.patt['y3'][i]-pat_ref['y3'][i]))



#pdb.set_trace()

print('******************************************************')
print('* TESTING: DETECTOR RELATIVE TO IDEAL                *')
print('******************************************************')
pp = LRSPattern(file=f)

print('Input pattern in frame {}\n\n'.format(pp.frame))

pp.to_coordinates(new_frame='idl')
print('Output pattern in frame {}\n\n'.format(pp.frame))
print(pp.patt)
print(' ')
pdb.set_trace()
print('Comparison with SIAF value for slit centre:')
ref_c_idl = ap.reference_point(to_frame='idl')
ref_nod1_tel = mt.xytov2v3(xx['nod1']['x']-1., xx['nod1']['y']-1., 'F770W')
ref_nod2_tel = mt.xytov2v3(xx['nod2']['x']-1., xx['nod2']['y']-1., 'F770W')

ref_nod1_idl = mt.v2v3toIdeal(ref_nod1_tel[0][0], ref_nod1_tel[1][0], 'MIRIM_SLIT')
ref_nod2_idl = mt.v2v3toIdeal(ref_nod2_tel[0][0], ref_nod2_tel[1][0], 'MIRIM_SLIT')

print ('DIFFERENCE IN SLIT CENTRE POS = {0:.5f} in X, {1:.5f} in Y'.format(pp.patt['x'][4]-ref_c_idl[0], pp.patt['y'][4]-ref_c_idl[1]))
print('DIFFERENCE IN NOD 1 POS = {0:.5f} in X, {1:.5f} in Y'.format(pp.patt['x2'][4]-ref_nod1_idl[0], pp.patt['y2'][4]-ref_nod1_idl[1]))
print('DIFFERENCE IN NOD 1 POS = {0:.5f} in X, {1:.5f} in Y'.format(pp.patt['x3'][4]-ref_nod2_idl[0], pp.patt['y3'][4]-ref_nod2_idl[1]))
print('\n\n')

print('Converting back to detector')
pp.to_coordinates(new_frame='det')


print('Ref pos 1 check:')
for i in range(len(pp.patt['x'])):
	print('Position {0}: delta_x = {1:.5f} px     delta_y = {2:.5f} px'.format(pp.patt['Pointing'][i], pp.patt['x'][i]-pat_ref['x'][i], pp.patt['y'][i]-pat_ref['y'][i]))

print('Ref pos 2 check:')
for i in range(len(pp.patt['x'])):
	print('Position {0}: delta_x = {1:.5f} px     delta_y = {2:.5f} px'.format(pp.patt['Pointing'][i], pp.patt['x2'][i]-pat_ref['x2'][i], pp.patt['y2'][i]-pat_ref['y2'][i]))
	
print('Ref pos 3 check:')
for i in range(len(pp.patt['x'])):
	print('Position {0}: delta_x = {1:.5f} px     delta_y = {2:.5f} px'.format(pp.patt['Pointing'][i], pp.patt['x3'][i]-pat_ref['x3'][i], pp.patt['y3'][i]-pat_ref['y3'][i]))





