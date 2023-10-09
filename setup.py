import os

from numpy.distutils.core import setup, Extension


def find_lapack_in_path(path_list):
    for path in path_list:
        lapack_path = os.path.join(path, 'liblapack.so')
        if os.path.exists(lapack_path):
            return lapack_path
    return None


ld_library_path = os.environ.get('LD_LIBRARY_PATH', '').split(':')
lapack_lib_path = find_lapack_in_path(ld_library_path)
ext_modules = [
    Extension(name='peslibf.o4_singlet', sources=['./src/O4_singlet.f90', './src/O4_singlet.pyf'], ),
    Extension(name='peslibf.n4_singlet',
              sources=['./src/PES_N4_singlet_umn_v3.f90', './src/PES_N4_singlet_umn_v3.pyf'], ),
    Extension(name='peslibf.o4_triplet', sources=['./src/O4_triplet_v2.f', './src/O4_triplet_v2.pyf']),
    Extension(name='peslibf.n2o2_triplet',
              sources=['./src/N2O2_3A_MB-PIP-MEG2.f90', './src/N2O2_3A_MB-PIP-MEG2.pyf'], ),
    Extension(name='peslibf.ch4oh', sources=['./src/ch4oh.pyf', './src/ch4oh.f', './src/lib/utility.f'], ),
    Extension(name='peslibf.h2o2', sources=['./src/h2o2.pyf', './src/h2o2.f'], ),
    Extension(name='peslibf.phoh', sources=['./src/phoh.pyf', './src/phoh_aprp.f'], libraries=['lapack'],
              library_dirs=[lapack_lib_path]),
    Extension(name='peslibf.phsch3', sources=['./src/phsch3.f90', './src/phsch3.pyf'], libraries=['lapack'],
              library_dirs=[lapack_lib_path]),
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
