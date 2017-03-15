from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# extension = Extension("pcUtils",
#                       sources      = ["pcUtils.pyx", "pc_utils.c"])

extension = Extension("pcUtils",
                      sources      = ["pcUtils.pyx"],
                      library_dirs = ['./'],
                      libraries    = ["pc_utils"])

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [extension])

# python setup_pc.py build_ext --inplace
# gcc -shared -O3 -Wall -fPIC -c pc_utils.c -o pc_utils.so
# gcc -shared -O3 -Wall -fPIC -fopenmp -c pc_utils.c -o pc_utils.so
