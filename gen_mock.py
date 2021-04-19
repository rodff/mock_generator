from astropy.io import fits
import numpy as np
import sys
import os

point = str(sys.argv[1])
bnd = str(sys.argv[2])
ID = str(sys.argv[3])

drt = 'point_{0}_{1}/gal_{2}'.format(point,bnd,ID)

# SAVING RANDOM STATE
np.random.seed()    # Set random state
arq = open('{0}/seed_noise.txt'.format(drt),'w')
arq.write('{0}'.format(np.random.get_state())) 
arq.close()

# Generate image for the model
os.system('makeimage {0}/true_{1}_{2}.dat --refimage {0}/cut_{1}_{2}.fits[1] --psf psf_gaussian_{3}_{2}.fits -o {0}/true_{1}_{2}.fits'.format(drt,ID,bnd,point))

# Load model data
mod = fits.getdata('{0}/true_{1}_{2}.fits'.format(drt,ID,bnd))

# Load stamp data
cut = fits.getdata('{0}/cut_{1}_{2}.fits'.format(drt,ID,bnd))
# Clip for negative values (if they happen to be)
mod = np.clip(mod,0.0,np.max(mod))
# Apply noise to model
n_mod = np.random.poisson(mod)

# Add noised model to stamp
mock = cut+n_mod

# Save
hdr = fits.getheader('{0}/cut_{1}_{2}.fits'.format(drt,ID,bnd),ext=0)
d_hdr = fits.getheader('{0}/cut_{1}_{2}.fits'.format(drt,ID,bnd),ext=1)

pars = dict(np.loadtxt('{0}/true_{1}_{2}.dat'.format(drt,ID,bnd),dtype=str))
d_hdr['re'] = float(pars['r_e'])
d_hdr['Ie'] = float(pars['I_e'])
d_hdr['n'] = float(pars['n'])
d_hdr['ell'] = float(pars['ell'])
d_hdr['PA'] = float(pars['PA'])

primary_hdu = fits.PrimaryHDU(header=hdr)
hdu = fits.ImageHDU(data=mock,header=d_hdr)     # Check the do_not_scale_image_data thing
hdul = fits.HDUList([primary_hdu,hdu])
hdul.writeto('{0}/mock_{1}_{2}.fits'.format(drt,ID,bnd),overwrite=True)
