from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = 'test',
    ext_modules = cythonize('prime_cython.pyx')
    )
