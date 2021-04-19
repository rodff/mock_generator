from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np
import sys

point = str(sys.argv[1])
bnd = str(sys.argv[2])
ID = str(sys.argv[3])

# INPUT
drt = 'point_{0}_{1}/gal_{2}'.format(point,bnd,ID)
big_fn = '../../point_{0}_{1}.fits'.format(point,bnd)
sz = 400
weight = True

# Limits of predefined region on the pointing (in degrees)
lim = np.loadtxt('point_{0}_{1}/coord_limits.dat'.format(point,bnd))
ra_min,ra_max = np.min(lim[:,0]),np.max(lim[:,0])
dec_min,dec_max = np.min(lim[:,1]),np.max(lim[:,1])

# Generate random ra, dec inside predefined region on the pointing
ra = np.random.uniform(ra_min,ra_max)
dec = np.random.uniform(dec_min,dec_max)

# Center of the cut in celestial coordinates
pos = SkyCoord(ra,dec, unit='deg', frame='fk5')

# Load big image
hdul = fits.open(big_fn)

# Locate tile
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

x,y = Owcs.all_world2pix(ra,dec,1)

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

        # Locate tile
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

        x,y = Owcs.all_world2pix(ra,dec,1)


if(weight==True):
    w_big_fn = '../../osw_point_{0}_{1}.fits'.format(point,bnd)
    w_hdul = fits.open(w_big_fn)
    w_data = w_hdul[Tn].data
    w_cut = Cutout2D(w_data, (x,y), sz, mode='strict')
    fits.writeto('{0}/wht_{1}_{2}.fits'.format(drt,ID,bnd),w_cut.data,overwrite=True)
else:
    pass

hdr = fits.getheader(big_fn,ext=0)
#d_hdr = hdul[Tn].header        # Problems using this header and matching WCS
primary_hdu = fits.PrimaryHDU(header=hdr)
hdu = fits.ImageHDU(data=cut.data)
hdu.header['RA'] = ra
hdu.header['DEC'] = dec
hdul = fits.HDUList([primary_hdu,hdu])
hdul.writeto('{0}/cut_{1}_{2}.fits'.format(drt,ID,bnd),overwrite=True)