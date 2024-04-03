from numpy.distutils.core import setup, Extension
import subprocess

result = subprocess.run(['make', 'clean'], cwd='src/CH2OH', stdout=subprocess.PIPE)
result = subprocess.run(['make'], cwd='src/CH2OH', stdout=subprocess.PIPE)
ext_modules = [
    # Extension(name='peslibf.o4_singlet', sources=['./src/O4_singlet.f90', './src/O4_singlet.pyf'], ),
    # Extension(name='peslibf.n4_singlet',
    #           sources=['./src/PES_N4_singlet_umn_v3.f90', './src/PES_N4_singlet_umn_v3.pyf'], ),
    # Extension(name='peslibf.o4_triplet', sources=['./src/O4_triplet_v2.f', './src/O4_triplet_v2.pyf']),
    # Extension(name='peslibf.n2o2_triplet',
    #           sources=['./src/N2O2_3A_MB-PIP-MEG2.f90', './src/N2O2_3A_MB-PIP-MEG2.pyf'], ),
    # Extension(name='peslibf.ch4oh', sources=['./src/ch4oh.pyf', './src/ch4oh.f', './src/lib/utility.f'], ),
    # Extension(name='peslibf.h2o2', sources=['./src/h2o2.pyf', './src/h2o2.f'], ),
    # Extension(name='peslibf.phoh', sources=['./src/phoh.pyf', './src/phoh_aprp.f'], libraries=['lapack', 'blas']),
    # Extension(name='peslibf.phsch3', sources=['./src/phsch3.f90', './src/phsch3.pyf'], libraries=['lapack', 'blas']),
    # Extension(name='peslibf.oh3', sources=['./src/oh3.f90', './src/oh3.pyf','./src/lib/oh3_util.f90'],
    # libraries=['lapack', 'blas']),
    # Extension(name='peslibf.nh3', sources=['./src/nh3.f', './src/nh3.pyf']),
    # Extension(name='peslibf.ch2oh',
    #           sources=['./src/CH2OH/evalsurf.f90', './src/CH2OH/ch2oh.pyf', ],
    #           libraries=['surfgen','lapack', 'blas', 'gomp'], library_dirs=['./src/CH2OH/','/usr/lib/x86_64-linux-gnu/'],
    #           extra_compile_args=['-fdefault-integer-8', '-fopenmp'],extra_c_compile_args=['--debug']
    #           ),

]
if __name__ == "__main__":
    setup(name='peslib',
          packages=['peslib'],
          ext_modules=ext_modules,
          install_requires=[
              'ase',
              'numpy',
          ]

          )
