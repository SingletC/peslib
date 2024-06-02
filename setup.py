from numpy.distutils.core import setup, Extension
import subprocess
import shutil
library_dirs = None #[]

'''
 or
 setenv LIBRARY_PATH "$LIBRARY_PATH":"blas":"lapack"
 setenv LIBRARY_PATH "$LIBRARY_PATH":"/mmfs1/data/tengcc/lib/BLAS-3.11.0/":"/mmfs1/data/tengcc/lib/lapack-3.10.0"                                                                                                                      bck:set
 setenv BLASDIR blas
 or
 module load blas lapack
'''
subprocess.run(['make', 'clean'], cwd='src/CH2OH', stdout=subprocess.PIPE)
subprocess.run(['make'], cwd='src/CH2OH', stdout=subprocess.PIPE)
shutil.copy("src/CH2OH/evalsurf.x","./peslib/")
ext_modules = [
    Extension(name='peslibf.o4_singlet', sources=['./src/O4_singlet.f90', './src/O4_singlet.pyf'], ),
    Extension(name='peslibf.n4_singlet',
              sources=['./src/PES_N4_singlet_umn_v3.f90', './src/PES_N4_singlet_umn_v3.pyf'], ),
    Extension(name='peslibf.o4_triplet', sources=['./src/O4_triplet_v2.f', './src/O4_triplet_v2.pyf']),
    Extension(name='peslibf.n2o2_triplet',
              sources=['./src/N2O2_3A_MB-PIP-MEG2.f90', './src/N2O2_3A_MB-PIP-MEG2.pyf'], ),
    Extension(name='peslibf.ch4oh', sources=['./src/ch4oh.pyf', './src/ch4oh.f', './src/lib/utility.f'], ),
    Extension(name='peslibf.h2o2', sources=['./src/h2o2.pyf', './src/h2o2.f'], ),
    Extension(name='peslibf.phoh', sources=['./src/phoh.pyf', './src/phoh_aprp.f'], libraries=['lapack', 'blas'],
              library_dirs=library_dirs),
    Extension(name='peslibf.phsch3', sources=['./src/phsch3.f90', './src/phsch3.pyf'], libraries=['lapack', 'blas'],
              library_dirs=library_dirs),
    Extension(name='peslibf.oh3', sources=['./src/oh3.f90', './src/oh3.pyf', './src/lib/oh3_util.f90'],
              libraries=['lapack', 'blas'],library_dirs=library_dirs),
    Extension(name='peslibf.nh3', sources=['./src/nh3.f', './src/nh3.pyf']),
    #           ),

]

setup(name='peslib',
          packages=['peslib'],
          package_data={'peslib':['data/**/*','evalsurf.x']},
          ext_modules=ext_modules,
          install_requires=[
              'ase',
              'numpy',
          ]

          )
