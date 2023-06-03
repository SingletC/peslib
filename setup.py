import glob
import os
from numpy.distutils.core import setup, Extension

ext_modules = [Extension(name='peslibf.' + os.path.splitext(os.path.basename(f))[0],
                         sources=[f], )
               for f in glob.glob('./src/*')]
if __name__ == "__main__":
    setup(name='peslib',
          packages=['peslib'],
          ext_modules=ext_modules, )
