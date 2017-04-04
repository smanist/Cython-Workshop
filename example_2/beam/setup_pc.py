from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

extension = Extension("pcUtils",
                      sources      = ["pcUtils.pyx", "pc_utils.c"])

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [extension])
