import sys
import os

point = str(sys.argv[1])
bnd = str(sys.argv[2])
Ni = int(sys.argv[3])
Nf = int(sys.argv[4])

ncp = '1'
ftnm = 'sersic_a'

for N in range(Ni,Nf+1):

    drt = 'point_{0}_{1}/gal_{2}'.format(point,bnd,N)

    print('#### \t gal_{0} \t ####'.format(N))

    print('Imfit profile fitting...\n')
    os.system('python imfit.py {0} {1} {2} {3} {4}'.format(point,bnd,N,ncp,ftnm))
    print('#####################################################################')
    print('#####################################################################')

    print('Computing total magnitude and surface brightness...\n')
    os.system('python mags.py {0} {1} {2} {3} {4}'.format(point,bnd,N,ncp,ftnm))
    print('#####################################################################')
    print('#####################################################################')
