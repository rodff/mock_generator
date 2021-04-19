import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
from astropy.io import fits

#defining fonts
import matplotlib
from matplotlib import rc
rc('font',**{'family':'times','size'   :25})
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rc('axes', labelsize=40)
matplotlib.rc('xtick', labelsize=30)
matplotlib.rc('ytick', labelsize=30)
matplotlib.rc('lines', lw=2.0,color='k')
matplotlib.rc('axes',lw=2.0)
matplotlib.rc('legend', fontsize=22)
plt.rcParams['figure.figsize'] = (12,10)


dmod = 29.93    # Distance modulus for NGC 3115

hdul = fits.open('carlsten_2020_data.fits')
dt = hdul[1].data
hdul.close()

px_s = 0.27     # DECAM pixel scale
pc_s = 46.94    # Parsec scale assuming distance modulus 29.93

re = dt['re_pc']
M = dt['Mg']

# To overcome ValueError: Big-endian buffer not supported on little-endian compiler
# Solution from https://tedboy.github.io/pandas/gotchas/gotchas10.html
re = re.byteswap().newbyteorder() # force native byteorder
M = M.byteswap().newbyteorder() # force native byteorder

print(M.shape,re.shape)

# Remove nan
idx = np.where(np.logical_not(np.isnan(re)))[0]
re,M = re[idx],M[idx]
# Remove nan
idx = np.where(np.logical_not(np.isnan(M)))[0]
re,M = re[idx],M[idx]
# Remove all values greater than 2.5 kpc
idx = np.where(re<=2500)[0]
re,M = re[idx],M[idx]
# Remove galaxies brighter than -16 mag
idx = np.where(M>=-16)[0]
re,M = re[idx],M[idx]

print(M.shape,re.shape)

# Converting radii from pc to pixels considering distance at NGC 3115 and DECAM pixel scale
#re = re/px_s/pc_s
#m = M+dmod

# Create linear regression object
regr = linear_model.LinearRegression()

# Train the model using the training sets
regr.fit(M.reshape(-1,1),np.log10(re.reshape(-1,1)))

Y = regr.predict(M.reshape(-1,1))
Y = Y.reshape(-1,)

deltas = Y-np.log10(re)
delta = np.max(np.abs(deltas))
#plt.hist(deltas,10)
#plt.show()

print('Score: \n', regr.score(M.reshape(-1,1),np.log10(re.reshape(-1,1))))
print('Mean Square Error: \n', mean_squared_error(M.reshape(-1,1),np.log10(re.reshape(-1,1))))
print('Max. delta: \n', delta)
print('Coefficients: \n', regr.coef_)
print('Intercepts: \n', regr.intercept_)

# Load results of the paper and plot alongside
df = pd.read_csv('../results_g.csv')
M_lsbs = df['m_g'].to_numpy()-dmod
Re_lsbs = df['Re_g_pc'].to_numpy()
#df = pd.read_csv('../uncertainties_g.csv')
#M_err = df['m_g_err'].to_numpy()
#Re_err = df['Re_g_pc_err'].to_numpy()

plt.scatter(M,np.log10(re),s=80,c='green',label='Carlsten et al. (2020)')
plt.scatter(M_lsbs,np.log10(Re_lsbs),s=175,marker='*',c='red',label='This work')
plt.plot(M,Y,c='black',lw=1.5)
#plt.plot(M,Y+delta,c='grey',lw=1)
#plt.plot(M,Y-delta,c='grey',lw=1)
plt.ylabel(r'log $\rm R_{e,g}$ [pc]')
plt.xlabel(r'$\rm M_g$ [mag]')
plt.legend()
plt.savefig('fit_re_M.png',dpi=300)
plt.show()
