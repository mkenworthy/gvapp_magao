from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt

amp   = fits.getdata('mapp1_aperture.fits')
phase = fits.getdata('mapp1_phase.fits')

# use amplitude to define the outer edges of the pupil
# and then crop to this size

nx, ny = amp.shape

(x,y) = np.mgrid[0:nx,0:ny]

xgood = x[(amp>0)]
ygood = y[(amp>0)]

amp_crop =   amp[xgood.min():xgood.max()+1,ygood.min():ygood.max()+1]
pha_crop = phase[xgood.min():xgood.max()+1,ygood.min():ygood.max()+1]

file_amp = 'mapp1_amp_crop.fits'
file_pha = 'mapp1_pha_crop.fits'

fits.writeto(file_amp, amp_crop, overwrite=True)
fits.writeto(file_pha, pha_crop, overwrite=True)

print('Written out cropped gvAPP pattern to {} and {}'.format(file_amp, file_pha)) 



