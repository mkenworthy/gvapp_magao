from hcipy import *
import numpy as np
from matplotlib import pyplot as plt
from rotate import *

from astropy.io import fits

amp = read_fits('mapp1_amp_crop.fits')
phase = read_fits('mapp1_pha_crop.fits')

clio_dims = (512, 1024) # dimensions of Clio2 detector

# add grating to the APP
split = 35. # separation between the two coronagraphic PSFs in diffraction widths
pupil_grid_for_grating = make_pupil_grid(amp.shape, 1)
phase = phase + (pupil_grid_for_grating.y*2*np.pi*(split/2.)).reshape(phase.shape)

# make the three complex pupils for the three PSFs
aper1 = amp * np.exp(-1j * phase)
aper2 = amp * np.exp(1j * phase)
aper3 = amp + 0j

D = 6.3 # effective telescope diameter [m]

pupil_grid = make_pupil_grid(aper1.shape, D)

npix = 512 # initial output for PSF before rotation and translation into Clio frame (should be 2048)
pscale = 0.015846 # pixel scale arcsec/pix Morsinski 2015 ApJ Appendix C.4

wavelength_0 = 3.94E-6     # reference wavelength

pixels_per_lambda = (206265. * wavelength_0 / D) / pscale

number_of_airy_rings = npix / pixels_per_lambda / 2# should give npix as the FOV

# This variable is used to determine the scale of a pixel in the simulation (effectively it is focal_length/D the fnumber of course...)
pixel_size = 18e-6 # microns
telescope_focal_length = 206265. / (pscale / pixel_size)

# Make the Telescope

aper11 = Field(aper1.ravel(), pupil_grid)
aper22 = Field(aper2.ravel(), pupil_grid)
aper33 = Field(aper3.ravel(), pupil_grid)

# This detector grid is in the units of telescope_focal_length
detector_grid = make_focal_grid(pupil_grid, wavelength=wavelength_0, q=pixels_per_lambda, num_airy=number_of_airy_rings, focal_length=telescope_focal_length)
prop = FraunhoferPropagator(pupil_grid, detector_grid, wavelength_0=wavelength_0, focal_length=telescope_focal_length)

def make_gvapp_psf(wavelength, leakage=0.03):
        # make the three PSFS at a single wavelength
        wf1 = Wavefront(aper11, wavelength)
        wf1_foc = prop.forward(wf1)
        wf2 = Wavefront(aper22, wavelength)
        wf2_foc = prop.forward(wf2)
        wf3 = Wavefront(aper33, wavelength)
        wf3_foc = prop.forward(wf3)

        # Plot it with the on_sky_grid coordinates
        psf_tot = wf1_foc.power + wf2_foc.power + leakage * wf3_foc.power

        psf_image = psf_tot/psf_tot.max()

        return psf_image


dy = 26
dx = -257
clio_rot = -26.
# Propagate different wavelengths

# lambda = 3.94um central, 90nm
wavelengths = np.linspace(3.94E-6, 5E-6, 1)

for wavelength in wavelengths:

        psf_image = make_gvapp_psf(wavelength)
        
        psf_orig = psf_image.reshape(npix,npix)
        
#        imshow_field(np.log10(psf_image), on_sky_grid, vmin=-5)

        write_fits(psf_image, 'unrotated.fits')

        t2 = cen_rot(psf_orig, clio_rot, clio_dims, offset1=(0.5*npix, 0.5 * npix),offset2=(dy,dx))

        sirius_flux = 20788.
        leakage_flux = 0.0726
        
        write_fits(t2*(sirius_flux/leakage_flux), 'test2.fits')

        sirius = read_fits('Sirius2_split_005200.fits')

        write_fits(sirius-t2*(sirius_flux/leakage_flux), 'test3.fits')

        
