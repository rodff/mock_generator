from astropy.io import fits
import numpy as np
import configparser
import sys
import os

cf_nm = str(sys.argv[1])
ID = str(sys.argv[2])

# CONFIGURATION
cfg = configparser.ConfigParser()
cfg.read(cf_nm)
D = dict(cfg['default'])
D_o = dict(cfg['optional'])   # Optional configurations
D_d = dict(cfg['data'])   # Optional configurations
outPath = str(D['outpath'])    # Output directory to save mock galaxies
drt = outPath+'/gal{0}'.format(ID)

# SAVING RANDOM STATE
np.random.seed()    # Set random state
arq = open('{0}/seed_noise_gal{1}.txt'.format(drt,ID),'w')
arq.write('{0}'.format(np.random.get_state())) 
arq.close()

# Generate image for the model
os.system('makeimage {0}/true_gal{1}.dat --refimage {0}/cut_gal{1}.fits[1] --psf {2} -o {0}/true_gal{1}.fits'.format(drt,ID,D['psf']))

# Load model data
mod = fits.getdata('{0}/true_gal{1}.fits'.format(drt,ID))

# Load stamp data
cut = fits.getdata('{0}/cut_gal{1}.fits'.format(drt,ID))
# Clip for negative values (if they happen to be)
mod = np.clip(mod,0.0,np.max(mod))
# Apply noise to model
n_mod = np.random.poisson(mod)

# Add noised model to stamp
mock = cut+n_mod

# Save
hdr = fits.getheader('{0}/cut_gal{1}.fits'.format(drt,ID),ext=0)
d_hdr = fits.getheader('{0}/cut_gal{1}.fits'.format(drt,ID),ext=1)

pars = dict(np.loadtxt('{0}/true_gal{1}.dat'.format(drt,ID),dtype=str))
d_hdr['re'] = float(pars['r_e'])
d_hdr['Ie'] = float(pars['I_e'])
d_hdr['n'] = float(pars['n'])
d_hdr['ell'] = float(pars['ell'])
d_hdr['PA'] = float(pars['PA'])

primary_hdu = fits.PrimaryHDU(header=hdr)
hdu = fits.ImageHDU(data=mock,header=d_hdr)     # Check the do_not_scale_image_data thing
hdul = fits.HDUList([primary_hdu,hdu])
hdul.writeto('{0}/mock_gal{1}.fits'.format(drt,ID),overwrite=True)
