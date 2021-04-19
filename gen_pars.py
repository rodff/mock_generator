from astropy.io import fits
from scipy.special import gamma
import numpy as np
import sys

def get_sampled_continuous(X,CDF):
    a = np.random.uniform(0,1)
    #   The cut using [1:] is to avoid 0 and -inf in CDF and X respectvely
    j = np.argmax(CDF[1:]>=a)-1
    xmin = X[1:][j]
    if(j==(X[1:].shape[0]-1)):
        x = np.random.uniform(X[1:][j-1],xmin)    # Check implications
    else:
        # Sample continuous values
        x = np.random.uniform(xmin,X[1:][j+1])
    return x

def logRe_pc(M):
    '''
    Return random effective radius (in log(re[pc])) based on correlation
    for linear fit to Carlsten+2020 data. M is absolute magnitude.
    '''
    a,b = -0.10526965,1.4593823
    delta = 0.3804598
    f = a*M+b + np.random.normal(0,delta/2.)
    return f

def bn(n):   # Using funtional forms of MacArthur et al. (2003)
    if(n>0.36):
        bnmod = (2*n)-(1./3.)+(46./(25515*n*n))+(131./(1148175*n*n*n))-(2194697./(30690717750*n*n*n*n))
    else:
        bnmod = 0.01945+(-0.8902*n)+(10.95*n*n)+(-19.67*n*n*n)+(13.43*n*n*n*n)
    return bnmod

def f(n):   # Using eq. 8 from Graham and Driver (2005)
    fn = np.exp(bn(n))*n*(bn(n)**(-2*n))*gamma(2*n)
    return fn

def I_e(re,m,zpt,n,ell): # Using same equation of GALFIT Manual (eq. 3)
    Ie = 10**((zpt-m)*2/5.)/(2*np.pi*(re**2)*f(n)*(1-ell))
    return Ie

point = str(sys.argv[1])
bnd = str(sys.argv[2])
ID = str(sys.argv[3])

# INPUT
drt = 'point_{0}_{1}/gal_{2}'.format(point,bnd,ID)
hdr = fits.getheader('{0}/cut_{1}_{2}.fits'.format(drt,ID,bnd))
dt = fits.getdata('{0}/cut_{1}_{2}.fits'.format(drt,ID,bnd))
zpt = float(hdr['MAGZERO'])
dmod = 29.93    # Distance modulus for NGC 3115
px_s = 0.27     # DECAM pixel scale
pc_s = 46.94    # Parsec scale assuming distance modulus 29.93

# SAVING RANDOM STATE
np.random.seed()    # Set random state
arq = open('{0}/seed_parameters.txt'.format(drt),'w')
arq.write('{0}'.format(np.random.get_state())) 
arq.close()

# Generate parameters
X0 = dt.shape[1]/2.
Y0 = dt.shape[0]/2.
PA = np.random.uniform(0,180)   # in degree
ell = np.random.uniform(0,0.75)
n = np.random.uniform(0.4,2)
'''
Res,Re_CDF = np.load('Re_carlsten_2020.npy'),np.load('CDF_Re_carlsten_2020.npy')
Re = get_sampled_continuous(Res,Re_CDF)
'''
Ms,M_CDF = np.load('M_g_carlsten_2020.npy'),np.load('CDF_M_g_carlsten_2020.npy')
M = get_sampled_continuous(Ms,M_CDF)
# Sample effective radius using correlation with magnitude
Re =(10**logRe_pc(M))/px_s/pc_s

m = M+dmod  # Apparent magnitude
Ie = I_e(Re,m,zpt,n,ell)  # Conversion to Ie

# Save true parameters file
arq = open('{0}/true_{1}_{2}.dat'.format(drt,ID,bnd),'w')
arq.write('X0 \t {0} \n'.format(X0))
arq.write('Y0 \t {0} \n'.format(Y0))
arq.write('FUNCTION Sersic\n')
arq.write('PA \t {0} \n'.format(PA))
arq.write('ell \t {0} \n'.format(ell))
arq.write('n \t {0} \n'.format(n))
arq.write('I_e \t {0} \n'.format(Ie))
arq.write('r_e \t {0} \n'.format(Re))
arq.close()
