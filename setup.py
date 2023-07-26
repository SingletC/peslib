import glob
import os
from numpy.distutils.core import setup, Extension

ext_modules = [
    Extension(name='peslibf.o4_singlet', sources=['./src/o4_singlet.f90'], ),
    Extension(name='peslibf.n4_singlet', sources=['./src/PES_N4_singlet_umn_v3.f90','./src/PES_N4_singlet_umn_v3.pyf'],
              ),
               ]
if __name__ == "__main__":
    setup(name='peslib',
          packages=['peslib'],
          ext_modules=ext_modules, )
