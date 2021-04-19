import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
from astropy.io import fits

hdul = fits.open('carlsten_2020_data.fits')
dt = hdul[1].data
hdul.close()

px_s = 0.27     # DECAM pixel scale
pc_s = 46.94    # Parsec scale assuming distance modulus 29.93

re = dt['re_pc']
re = re[np.logical_not(np.isnan(re))]   # Remove nan
re = re[re<=2500] # Remove all values greater than 2.5 kpc
# Converting radii from pc to pixels considering distance at NGC 3115 and DECAM pixel scale
re = re/px_s/pc_s

plt.hist(re)
plt.show()

ecdf = ECDF(re)
np.save('Re_carlsten_2020.npy',ecdf.x)
np.save('CDF_Re_carlsten_2020.npy',ecdf.y)

M = dt['Mg']
M = M[np.logical_not(np.isnan(M))]   # Remove nan
M = M[M>=-16]

ecdf = ECDF(M,side='left')

plt.hist(M)
plt.show()

np.save('M_g_carlsten_2020.npy',ecdf.x)
np.save('CDF_M_g_carlsten_2020.npy',ecdf.y)
