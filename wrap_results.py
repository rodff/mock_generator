from astropy.io import fits
import pandas as pd
import numpy as np
import os
import sys

point = str(sys.argv[1])
bnd = str(sys.argv[2])
Ni = 0
Nf = int(sys.argv[3])

ncp = '1'
ftnm = 'sersic_a'

rtd = 'point_{0}_{1}'.format(point,bnd) # Root dir

#drts = [f.name for f in os.scandir('./{0}'.format(rtd)) if f.is_dir()]

# Lists for parameters
ID_l,RA_l,DEC_l = [],[],[]
ell_l,n_l,Ie_l,re_l = [],[],[],[]       # Sersic model parameters
ellI_l,nI_l,IeI_l,reI_l = [],[],[],[]   # Initial guesses for model parameters
ellR_l,nR_l,IeR_l,reR_l = [],[],[],[]   # Recovered model parameters
m_l,mu_0_l,mu_e_l,av_mu_e_l = [],[],[],[]
mR_l,mu_0R_l,mu_eR_l,av_mu_eR_l = [],[],[],[]   # Recovered magnitudes
rec_l = []    # Boolean to say if the model was successfuly recovered (0=no,1=yes)

# Loop over galaxies
for N in range(Ni,Nf+1):
    drt = 'gal_{0}'.format(N)
    ID = int(drt[4:])
    print('Galaxy {0}'.format(ID))
    # Coordinates
    hdr = fits.getheader('{0}/{1}/mock_{2}_{3}.fits'.format(rtd,drt,ID,bnd),ext=1)
    RA,DEC = float(hdr['RA']),float(hdr['DEC'])

    # Sersic model fitted parameters
    pars = dict(np.loadtxt('{0}/{5}/pars_{1}_{2}_{3}_{4}.dat'.format(rtd,ID,bnd,ncp,ftnm,drt),dtype=str))
    reR = float(pars['r_e'])
    IeR = float(pars['I_e'])
    nR = float(pars['n'])
    ellR = float(pars['ell'])

    # Sersic model initial guesses
    pars = np.loadtxt('{0}/{5}/conf_{1}_{2}_{3}_{4}.dat'.format(rtd,ID,bnd,ncp,ftnm,drt),dtype=str,comments=['#','FUNCTION'])
    pars = dict(pars[2:-1][:,:-1])
    reI = float(pars['r_e'])
    IeI = float(pars['I_e'])
    nI = float(pars['n'])
    ellI = float(pars['ell'])

    # Sersic model true parameters
    pars = dict(np.loadtxt('{0}/{3}/true_{1}_{2}.dat'.format(rtd,ID,bnd,drt),dtype=str))
    re = float(pars['r_e'])
    Ie = float(pars['I_e'])
    n = float(pars['n'])
    ell = float(pars['ell'])

    # Computed Magnitudes
    f = pd.read_csv('{0}/{3}/mags_{1}_{2}.csv'.format(rtd,ID,bnd,drt))
    mR = float(f['mag_tot'][0])
    mu_0R = float(f['mu_0'][0])
    mu_eR = float(f['mu_e'][0])
    av_mu_eR = float(f['mean_mu_e'][0])

    # True Magnitudes
    f = pd.read_csv('{0}/{3}/true_mags_{1}_{2}.csv'.format(rtd,ID,bnd,drt))
    m = float(f['mag_tot'][0])
    mu_0 = float(f['mu_0'][0])
    mu_e = float(f['mu_e'][0])
    av_mu_e = float(f['mean_mu_e'][0])

    # Save to lists
    ID_l.append(ID)
    RA_l.append(RA)
    DEC_l.append(DEC)
    ell_l.append(ell)
    ellR_l.append(ellR)
    ellI_l.append(ellI)
    re_l.append(re)
    reR_l.append(reR)
    reI_l.append(reI)
    n_l.append(n)
    nR_l.append(nR)
    nI_l.append(nI)
    Ie_l.append(Ie)
    IeR_l.append(IeR)
    IeI_l.append(IeI)
    m_l.append(m)
    mR_l.append(mR)
    mu_0_l.append(mu_0)
    mu_0R_l.append(mu_0R)
    mu_e_l.append(mu_e)
    mu_eR_l.append(mu_eR)
    av_mu_e_l.append(av_mu_e)
    av_mu_eR_l.append(av_mu_eR)

# Create Data Frame
df = pd.DataFrame({'ID':ID_l,'RA':RA_l,'DEC':DEC_l,
'ell':ell_l,'ell_f':ellR_l,'ell_i':ellI_l,'n':n_l,'n_f':nR_l,'n_i':nI_l,
're':re_l,'re_f':reR_l,'re_i':reI_l,'Ie':Ie_l,'Ie_f':IeR_l,'Ie_i':IeI_l,
'm':m_l,'m_f':mR_l,'mu_0':mu_0_l,'mu_0_f':mu_0R_l,'mu_e':mu_e_l,'mu_e_f':mu_eR_l,
'mean_mu_e':av_mu_e_l,'mean_mu_e_f':av_mu_eR_l})

# Save to CSV file
df.to_csv('point_{0}_{1}.csv'.format(point,bnd),index=False)
