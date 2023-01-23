from setuptools import setup,Extension
from Cython.Build import cythonize
import numpy


cython_directives = {'language_level': "3"}

setup(
   ext_modules=cythonize(
                    Extension(
                             name='Calculation'
                            ,sources=["Calculation.pyx"]
                            ,include_dirs=[numpy.get_include()]
                               )
                            ,include_path=["/home/nehleh/anaconda3/envs/PhiloBacteria/lib/python3.8/site-packages/numpy/core/include/numpy"]
                            ,compiler_directives= cython_directives
                            )
)

