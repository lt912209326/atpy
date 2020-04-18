from distutils.core import setup, Extension
from setuptools import setup, Extension
from Cython.Build import cythonize


extension = Extension(name="cppelement",
               sources=['cppelement.pyx'],
               language="c++")

setup(ext_modules=cythonize(extension,
                            compiler_directives={'language_level':3},
                            include_path=['.','../include']),
                            package_data={'../include': ['cppelement.pxd']},
                            package_dir = {'elements': '../include'})
# setup(
    # name='my_cython_package',
    # packages=['core'],
    # package_data={'core': ['cppelement.pxd']},  # or wherever the files are
    # ...etc...

                            # package_data=['../include':'cppelement.pxd'])