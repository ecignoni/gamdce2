from setuptools import setup, Extension
# from numpy.distutils.core import setup
# from numpy.distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

def fortran_compiler(statement):
    print(statement)
    system(statement)

# compile the fortran modules without linking
# compile fortran modules
fortran_compiler('cd src ; gfortran mod_types.f90 -c -o mod_types.o -O3 -fPIC')
fortran_compiler('cd src ; gfortran mod_gamdce2.f90 -c -o mod_gamdce2.o -O3 -fPIC')

# compile fortran shared modules
fortran_compiler('cd src; gfortran mod_c_interface.f90 -c -o mod_c_interface.o -O3 -fPIC')

# editable components
module_name = "gamdce2"
source_file = "src/py_interface.pyx"
compiled_modules = ['src/mod_gamdce2.o', 'src/mod_types.o', 'src/mod_c_interface.o']

# setup
ext_modules = [Extension(# module name:
                         module_name,
                         # source file:
                         [source_file],
                         libraries=["gfortran"],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                         # other files to link to
                         extra_link_args=compiled_modules)]

for e in ext_modules:
    e.cython_directives = {'language_level': "3"} #all are Python-3

setup(name = module_name,
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)

