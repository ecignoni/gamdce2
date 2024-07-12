from os import system
import numpy as np
from setuptools import setup, Extension
from Cython.Build import cythonize

def fortran_compiler(statement):
    print(statement)
    system(statement)

# compile the fortran modules without linking
# compile fortran modules
fortran_compiler('cd src ; gfortran mod_types.f90 -c -o mod_types.o -O3 -fPIC')
fortran_compiler('cd src ; gfortran mod_gamdce2.f90 -c -o mod_gamdce2.o -O3 -fPIC')

# compile fortran shared modules
fortran_compiler('cd src; gfortran mod_c_interface.f90 -c -o mod_c_interface.o -O3 -fPIC')

# setup
ext_modules = [(
    Extension(
        "gamdce2.py_interface",
        ["src/py_interface.pyx"],
        libraries=["gfortran"],
        # other compile args for gcc
        extra_compile_args=["-O3", "-fPIC"],
        extra_f90_compile_args=['-fPIC', '-O3'],
        include_dirs=[np.get_include()],
        # other files to link to
        extra_link_args=['src/mod_gamdce2.o', 'src/mod_types.o', 'src/mod_c_interface.o'],
    )
)]

# for e in ext_modules:
#     e.cython_directives = {'language_level': "3"}

packages = [
    "gamdce2",
]

setup(
    name="gamdce2",
    version="0.1.0",
    author="Edoardo Cignoni",
    author_email="edoardo.cignoni96@gmail.com",
    packages=packages,
    ext_modules = cythonize(ext_modules, compiler_directives={"language_level": "3"}),
    description="Fast reweighting of GaMD simulations with 2nd order Cumulant Expansion",
    long_description=open("README.md").read(),
    setup_required=["numpy", "cython"],
    install_required=["numpy", "cython", "pandas"],
    #cmdclass = {'build_ext': build_ext},
    #include_dirs = [get_include()],
)
