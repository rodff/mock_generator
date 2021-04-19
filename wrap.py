import time
import sys
import os
import configparser

cf_nm = str(sys.argv[1])
Ni = int(sys.argv[2])   # Starting ID number for mock galaxies
Nf = int(sys.argv[3])   # Ending ID number for mock galaxies

# CONFIGURATION
cfg = configparser.ConfigParser()
cfg.read(cf_nm)
D = dict(cfg['default'])
D_o = dict(cfg['optional'])   # Optional configurations
D_d = dict(cfg['data'])   # Optional configurations
outPath = str(D['outpath'])    # Output directory to save mock galaxies

ncp = '1'
ftnm = 'sersic_a'

start_time = time.time()

for N in range(Ni,Nf+1):

    drt = outPath+'/gal{0}'.format(N)

    print('#### \t gal{0} \t ####'.format(N))
    if(os.path.exists(drt)):
        os.system('rm -r {0}'.format(drt))
        os.system('mkdir {0}'.format(drt))
    else:
        os.system('mkdir {0}'.format(drt))

    print('Cutting stamp from pointing...\n')
    os.system('python cut.py {0} {1}'.format(cf_nm,N))
    print('#####################################################################')

    print('Generating random parameters for mock galaxy...\n')
    os.system('python gen_pars.py {0} {1}'.format(cf_nm,N))
    print('#####################################################################')
    
    print('Create model image, apply noise, convolve with PSF...\n')
    os.system('python gen_mock.py {0} {1}'.format(cf_nm,N))
    print('#####################################################################')

t = time.time() - start_time
print("--- {0} minutes ---".format(t/60.)) 