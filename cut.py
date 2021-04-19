from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np
import configparser
import sys

cf_nm = str(sys.argv[1])
ID = str(sys.argv[2])

# CONFIGURATION
cfg = configparser.ConfigParser()
cfg.read(cf_nm)
D = dict(cfg['default'])
D_o = dict(cfg['optional'])   # Optional configurations

sz = int(D['size'])     # Size of the cutout in pixels
outPath = str(D['outpath'])    # Output directory to save mock galaxies
big_fn = str(D['image'])    # Image used as background for the mock galaxies
drt = outPath+'/gal{0}'.format(ID)

# Limits of predefined square region to create mocks (in degrees)
ra_min,ra_max = float(D['ra_min']),float(D['ra_max'])
dec_min,dec_max = float(D['dec_min']),float(D['dec_max'])

# Generate random ra, dec inside predefined region 
ra = np.random.uniform(ra_min,ra_max)
dec = np.random.uniform(dec_min,dec_max)
# Center of the cut in celestial coordinates
pos = SkyCoord(ra,dec, unit='deg', frame='fk5')

# Load big image
hdul = fits.open(big_fn)
if(str(D['multi_ext'])=='yes'):
    # Locate extension with corresponding location
    for i in range(1,len(hdul)):
        Twcs = WCS(hdul[i].header,hdul)
        if(Twcs.footprint_contains(pos)):
            Owcs = Twcs
            data = hdul[i].data
            Tn = i
            print('Object in extension {0}'.format(Tn))
            break
        else:
            pass
else:
    Twcs = WCS(hdul[0].header,hdul)
    data = hdul[0].data

x,y = Owcs.all_world2pix(ra,dec,1)

# While to force that the cut is entirely composed of valid pixels
bool_cut = False
while(bool_cut==False):
    try:
        # Cut image
        cut = Cutout2D(data, (x,y), sz, mode='strict')
        bool_cut=True
    except:
        # Generate random ra, dec inside predefined region on the pointing
        ra = np.random.uniform(ra_min,ra_max)
        dec = np.random.uniform(dec_min,dec_max)
        # Center of the cut in celestial coordinates
        pos = SkyCoord(ra,dec, unit='deg', frame='fk5')

        if(str(D['multi_ext'])=='yes'):
            # Locate extension with corresponding location
            for i in range(1,len(hdul)):
                Twcs = WCS(hdul[i].header,hdul)
                if(Twcs.footprint_contains(pos)):
                    Owcs = Twcs
                    data = hdul[i].data
                    Tn = i
                    print('Object in extension {0}'.format(Tn))
                    break
                else:
                    pass
        else:
            Twcs = WCS(hdul[0].header,hdul)
            data = hdul[0].data

        x,y = Owcs.all_world2pix(ra,dec,1)

# If a weight image is specified
if('wht_image' in D_o):
    w_big_fn = str(D_o['wht_image'])
    w_hdul = fits.open(w_big_fn)
    w_data = w_hdul[Tn].data
    w_cut = Cutout2D(w_data, (x,y), sz, mode='strict')
    fits.writeto('{0}/wht_gal{1}.fits'.format(drt,ID),w_cut.data,overwrite=True)
else:
    pass

# Just copying the header of the big image into the cutout
hdr = fits.getheader(big_fn,ext=0)
primary_hdu = fits.PrimaryHDU(header=hdr)
hdu = fits.ImageHDU(data=cut.data)
hdu.header['RA'] = ra
hdu.header['DEC'] = dec
hdul = fits.HDUList([primary_hdu,hdu])
hdul.writeto('{0}/cut_gal{1}.fits'.format(drt,ID),overwrite=True)
print('Saved at {0}/cut_gal{1}.fits'.format(drt,ID))