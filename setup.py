# from distutils.core import setup, Extension
from setuptools import setup, Extension,find_packages
from Cython.Build import cythonize
from Cython.Distutils import build_ext


extension = [
            Extension("atpy.src.cython.constants",
                    ['src/cython/constants.pyx'],
                    include_dirs = [".",'../cython','src/cpp'],
                    language="c++",
                    extra_compile_args = ["-O3", "-Wall"]),
            Extension("atpy.src.cython.elements",
                    ['src/cython/elements.pyx'],
                    include_dirs = ["./cython",'src/cpp'],
                    language="c++",
                    extra_compile_args = ["-O3", "-Wall"]),
            Extension("atpy.src.cython.structures",
                    ['src/cython/structures.pyx'],
                    include_dirs = [".",'../cython','src/cpp'],
                    language="c++",
                    extra_compile_args = ["-O3", "-Wall"]),
            Extension("atpy.src.cython.beamline",
                    ['src/cython/beamline.pyx'],
                    include_dirs = ["./cython",'src/cpp'],
                    language="c++",
                    extra_compile_args = ["-O3", "-Wall"]),
            Extension("atpy.src.cython.parser",
                    ['src/cython/parser.pyx'],
                    include_dirs = ["./cython"],
                    language="c++",
                    extra_compile_args = ["-O3", "-Wall"]),
            Extension("atpy.src.cython.lattice",
                    ['src/cython/lattice.pyx'],
                    include_dirs = ["./cython",'src/cpp'],
                    language="c++",
                    extra_compile_args = ["-O3", "-Wall"]),
            ]

setup(name='atpy',
      packages=find_packages(),
      cmdclass = {'build_ext': build_ext},
      ext_modules=cythonize(extension,
                            compiler_directives={'profile':True, 'language_level':3,
                                                'boundscheck':False, 'cdivision':True, 
                                                'initializedcheck':False }),
      package_dir={'':'..'},
      zip_safe = False,
     )
