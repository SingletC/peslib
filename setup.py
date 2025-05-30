import os.path

from numpy.distutils.core import setup, Extension
import subprocess
import shutil
from setuptools.command.install import install
import os
import stat
import sys

library_dirs = None  #[]

'''
 or
 setenv LIBRARY_PATH "$LIBRARY_PATH":"blas":"lapack"
 setenv LIBRARY_PATH "$LIBRARY_PATH":"/mmfs1/data/tengcc/lib/BLAS-3.11.0/":"/mmfs1/data/tengcc/lib/lapack-3.10.0"
 setenv BLASDIR blas
 or
 module load blas lapack
'''
# note for andromeda bc user,
"""
setenv LIBRARY_PATH "$LIBRARY_PATH":"/mmfs1/data/tengcc/lib/BLAS-3.11.0/"
module load lapack
"""
try:
    print("Cleaning previous build...")
    clean_result = subprocess.run(['make', 'clean'], cwd='src/CH2OH', stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if clean_result.returncode != 0:
        print(f"Warning during make clean: {clean_result.stderr}")
    
    print("Compiling evalsurf.x...")
    make_result = subprocess.run(['make'], cwd='src/CH2OH', stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if make_result.returncode != 0:
        print("Compilation failed with the following errors:")
        print(f"STDOUT: {make_result.stdout}")
        print(f"STDERR: {make_result.stderr}")
        raise subprocess.CalledProcessError(make_result.returncode, 'make')
    
    if not os.path.exists('src/CH2OH/evalsurf.x'):
        print("Compilation appeared to succeed but evalsurf.x was not created.")
        print(f"STDOUT: {make_result.stdout}")
        print(f"STDERR: {make_result.stderr}")
        raise FileNotFoundError('evalsurf.x compile failed')
    
    print("Compilation successful, copying evalsurf.x...")
    shutil.copy2("src/CH2OH/evalsurf.x","./peslib/")
except (subprocess.CalledProcessError, FileNotFoundError) as e:
    print(f"Error during compilation: {str(e)}")
    print("Please check if all dependencies are installed (BLAS, LAPACK, etc.)")
    print("You may need to load required modules or set environment variables:")
    print("  module load blas lapack")
    print("  export LIBRARY_PATH=$LIBRARY_PATH:/path/to/blas:/path/to/lapack")
    sys.exit(1)

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
    Extension(name='peslibf.oh3', sources=['./src/OH3/pes.f90', './src/OH3/oh3.pyf',
                                           './src/OH3/adiabats.f90',
                                           './src/OH3/dcdx_u12.f90',
                                           './src/OH3/dcdx_u13.f90',
                                           './src/OH3/dcdx_u23.f90',
                                           './src/OH3/diabaticcouplings.f90',
                                           './src/OH3/diabats.f90',
                                           ],
              libraries=['lapack', 'blas'], library_dirs=library_dirs),
    Extension(name='peslibf.o2h3', sources=['./src/O2H3/h3o2-pipnn.f', './src/O2H3/o2h3.pyf']),
    Extension(name='peslibf.nh3', sources=['./src/nh3.f', './src/nh3.pyf']),

]


class CustomInstallCommand(install):

    def run(self):
        install.run(self)

        install_dir = self.install_lib
        target_file = os.path.join(install_dir, 'peslib', 'evalsurf.x')
        if os.path.exists(target_file):
            print(f"Modifying permissions for {target_file}")
            os.chmod(target_file, stat.S_IRWXU)
        else:
            print(f"Warning: {target_file} not found.")

setup(name='peslib',
      packages=['peslib'],
      package_data={'peslib': ['data/**/*', 'evalsurf.x']},
      ext_modules=ext_modules,
      install_requires=[
          'ase',
          'numpy',
      ],
      cmdclass={
        'install': CustomInstallCommand,
    },
      )