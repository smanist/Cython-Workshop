from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = 'cyUtils',
    ext_modules = cythonize('cyUtils.pyx')
    )
