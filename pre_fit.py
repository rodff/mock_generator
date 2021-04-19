from astropy.io import fits
from scipy import stats
import numpy as np
import subprocess
import sys
import os

def RunSex(gal_img: str,ID:str, bnd: str, point: str, out_dir=os.getcwd(), sex_dir=os.getcwd()):
    # Defining arguments for Sextractor
    cf_nm = '{0}/default_{2}_{1}.sex'.format(sex_dir,bnd,point)  # Sextractor conf. file
    sg_nm = '{2}/seg_{0}_{1}.fits'.format(ID,bnd,out_dir) # SEGMENTATION image filename
    back_nm = '{2}/back_{0}_{1}.fits'.format(ID,bnd,out_dir) # BACKGROUND image filename

    # Call sextractor
    os.system('sextractor '+gal_img+' -c '+cf_nm+' -CATALOG_NAME '+out_dir+'/'+ID+'_'+bnd+
    '.cat -CHECKIMAGE_NAME '+sg_nm+','+back_nm)

def GenMask(img_fn: str, seg_fn: str, label: str, center=None, out_dir=os.getcwd(), sex_dir=os.getcwd()):

        seg_dt = fits.getdata(seg_fn)  # Segmentation array
        img_dt = fits.getdata(img_fn)   # Data array

        if(center==None):
            # Center of the galaxy in the image (specifically if the galaxy is in the center of the stamp)
            X0,Y0 = int(img_dt.shape[1]/2.),int(img_dt.shape[0]/2.)
        else:
            X0,Y0 = int(center[0]),int(center[1])

        R_I = seg_dt[Y0][X0]

        # Remove galaxy position from the construction of the mask
        # (Turning the intensity of the pixels into zero)

        print('Iterating image... \n')

        # Iteration over the array to remove galaxy from the masked objects
        for cl in np.nditer(seg_dt, op_flags=['readwrite']):
            if cl[...]==R_I:
                cl[...]=0       # If the intensity matches galaxy segmentation, zero it
            else:
                pass            # If not, just go to the next pixel

        print("Saving mask image... \n")

        # Save mask image
        fits.writeto('{0}/msk_{1}.fits'.format(out_dir,label), seg_dt, overwrite=True)


def ParNoise(pmin,pmax,deltaV,val):
    ns = np.random.uniform(-1,1)*deltaV
    n_val = val+ns
    if(n_val<pmin):
        return pmin+np.abs(ns)
    elif(n_val>pmax):
        return pmax-np.abs(ns)
    else:
        return n_val

def SkyStatistics(back_fn):
    '''
    DISCLAIMER: This function assumes that the background is well estimated by
    SExtractor and that the Sextractor mode equation is a good estimator for the
    background level, this would be wrong if there is a bright source influencing
    the background map, for example.
    '''
    m_dt = fits.getdata(back_fn)
    # Mode estimation from https://sextractor.readthedocs.io/en/latest/Background.html
    sex_mode = 2.5*np.median(m_dt)-1.5*np.mean(m_dt)
    sky_mode = stats.mode(m_dt,axis=None)
    sky_median = np.median(m_dt)
    sky_mean = np.mean(m_dt)
    sky_std = np.std(m_dt)
    return sex_mode,sky_mode[0],sky_median,sky_mean,sky_std

point = str(sys.argv[1])
bnd = str(sys.argv[2])
ID = str(sys.argv[3])

drt = 'point_{0}_{1}/gal_{2}'.format(point,bnd,ID)

img_fn = '{0}/mock_{1}_{2}.fits'.format(drt,ID,bnd)
seg_fn = '{0}/seg_{1}_{2}.fits'.format(drt,ID,bnd)
back_fn = '{0}/back_{1}_{2}.fits'.format(drt,ID,bnd)
pars = dict(np.loadtxt('{0}/true_{1}_{2}.dat'.format(drt,ID,bnd),dtype=str))

# Run sextractor
RunSex(img_fn,ID,bnd,point,out_dir=drt)
# Generate Mask
GenMask(img_fn,seg_fn,'{0}_{1}'.format(ID,bnd),out_dir=drt)

# Estimate different values for the sky using background image of Sextractor
sky_sex,sky_mod,sky_med,sky_av,sky_std = SkyStatistics(back_fn)
Isky = sky_med

arq = open('{0}/sky_{1}_{2}.csv'.format(drt,ID,bnd),'w')
arq.write('sky_sex,back_mode,back_median,back_mean,back_std \n')
arq.write('{0},{1},{2},{3},{4}'.format(sky_sex,sky_mod[0],sky_med,sky_av,sky_std))
arq.close()

# Get true values and add perturbation
X0 = float(pars['X0']) + np.random.uniform(-2,2)
Y0 = float(pars['Y0']) + np.random.uniform(-2,2)
re = ParNoise(2,200,float(pars['r_e'])/10.,float(pars['r_e']))
Ie = ParNoise(0.00005,150000,float(pars['I_e'])/5.,float(pars['I_e']))
n = ParNoise(0.35,2.5,0.5,float(pars['n']))
ell = ParNoise(0,0.85,0.3,float(pars['ell']))
PA = ParNoise(0,180,15,float(pars['PA']))

# Create configuration file with perturbed values
arq = open('{0}/conf_{1}_{2}_1_sersic_a.dat'.format(drt,ID,bnd),'w')
arq.write('X0 \t {0} \t {1},{2} \n'.format(X0,X0-3,X0+3))
arq.write('Y0 \t {0} \t {1},{2} \n'.format(Y0,Y0-3,Y0+3))
arq.write('FUNCTION Sersic\n')
arq.write('PA \t {0} \t 0,180 \n'.format(PA))
arq.write('ell \t {0} \t 0,0.85 \n'.format(ell))
arq.write('n \t {0} \t 0.3,2.5 \n'.format(n))
arq.write('I_e \t {0} \t 0.00005,150000 \n'.format(Ie))
arq.write('r_e \t {0} \t 2,200 \n'.format(re))
arq.write('FUNCTION FlatSky\n')
arq.write('I_sky \t {0} \t {1},{2} \n'.format(Isky,Isky-(sky_std/2),Isky+(sky_std/2)))
arq.close()
