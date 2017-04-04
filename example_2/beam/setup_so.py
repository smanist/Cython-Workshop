from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

extension = Extension("pcUtils",
                      sources      = ["pcUtils.pyx"],
                      library_dirs = ['./'],
                      libraries    = ["pc_utils"])

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [extension])
