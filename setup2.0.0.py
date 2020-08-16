# from distutils.core import setup, Extension
from setuptools import setup, Extension,find_packages
from Cython.Build import cythonize
from Cython.Distutils import build_ext


extension = [Extension("atpy.src.cython.elements",
                    ['src/cython/elements.pyx'],
                    include_dirs = ["./cython"],
                    language="c++",
                    extra_compile_args = ["-O3", "-Wall"]),
            Extension("atpy.src.cython.variables",
                    ['src/cython/variables.pyx'],
                    include_dirs = ["./cython"],
                    language="c++",
                    extra_compile_args = ["-O3", "-Wall"]),
            Extension("atpy.src.cython.constraints",
                    ['src/cython/constraints.pyx'],
                    include_dirs = ["./cython"],
                    language="c++",
                    extra_compile_args = ["-O3", "-Wall"]),
            Extension("atpy.src.cython.optima",
                    ['src/cython/optima.pyx'],
                    include_dirs = ["./cython"],
                    language="c++",
                    extra_compile_args = ["-O3", "-Wall"]),
            Extension("atpy.src.cython.lattice",
                    ['src/cython/lattice.pyx'],
                    #'src/cython/elements.pyx','src/cython/variables.pyx','src/cython/constraints.pyx','src/cython/optima.pyx'],
                    include_dirs = ["./cython"],
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
      # packages = ['atpy.src.cpp','atpy.src.cython']
      # package_data={'./include': ['cppelement.pxd']},
      # include_path=['.core','.src'],
      # compiler_directives={'language_level':3},
      # cythonize(extension,
                            # compiler_directives={'language_level':3},
                            # include_path=['./core','../include'])
     )
# setup(
    # name='my_cython_package',
    # packages=['core'],
    # package_data={'core': ['cppelement.pxd']},  # or wherever the files are
    # ...etc...

                            # package_data=['../include':'cppelement.pxd'])
