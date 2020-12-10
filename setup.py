# from distutils.core import setup, Extension
from setuptools import setup, Extension,find_packages
from Cython.Build import cythonize
from Cython.Distutils import build_ext


link_args = ['-static-libgcc',
             '-static-libstdc++',
             '-Wl,-Bstatic,--whole-archive',
             '-lwinpthread',
             '-Wl,--no-whole-archive']

#Extension args
language ='c++'
include_dirs=[".",'atpy/src/','atpy/src/cpp', 'tpsa/include']

extra_compile_args = ["-O3", "-Wall"]

std_cpp11 = ["-O3", "-Wall", "-std=c++11"]

sources=[{'name':"atpy.src.constants", 'source':['atpy/src/constants.pyx'],   },
        {'name':"atpy.src.tpsa.pytpsa"      , 'source':['atpy/src/tpsa/pytpsa.pyx','atpy/src/tpsa/src/da.cpp','atpy/src/tpsa/src/tpsa_extend.cpp'], 
          'extra_compile_args':std_cpp11 },
        {'name':"atpy.src.element.elements"  , 'source':['atpy/src/element/elements.pyx'],    },
        {'name':"atpy.src.beamline.structures", 'source':['atpy/src/beamline/structures.pyx'],  },
        {'name':"atpy.src.beamline.beamline"  , 'source':['atpy/src/beamline/beamline.pyx'],    },
        {'name':"atpy.src.parser.parser"    , 'source':['atpy/src/parser/parser.pyx'],      },
        {'name':"atpy.src.beamline.lattice"   , 'source':['atpy/src/beamline/lattice.pyx'],     },
        {'name':"atpy.test_cimport"           , 'source':['atpy/test_cimport.pyx'],             },
        ]
def set_extension():
    extension=[]
    for pkg in sources:
        extension.append( Extension(name=pkg['name'], sources=pkg['source'], 
                                    language=language, 
                                    include_dirs=include_dirs, 
                                    extra_compile_args=extra_compile_args if 'extra_compile_args' not in pkg.keys() else pkg['extra_compile_args']) 
                        )
    return extension


# extension = [
            # Extension("atpy.src.cython.constants",
                    # ['src/cython/constants.pyx'],
                    # include_dirs = [".",'../cython','src/cpp'],
                    # # language="c++",
                    # extra_compile_args = ["-O3", "-Wall"]),
            # Extension("atpy.src.tpsa.pytpsa",
                    # ['src/tpsa/pytpsa.pyx','src/tpsa/src/da.cpp','src/tpsa/src/tpsa_extend.cpp'],#,'src/tpsa/src/tpsa.cpp'
                    # include_dirs = ['../tpsa/include'],
                    # # language="c++",
                    # link_args=link_args,
                    # extra_compile_args = ["-O3", "-Wall","-std=c++11"]),
            # Extension("atpy.src.cython.elements",
                    # ['src/cython/elements.pyx'],
                    # include_dirs = ["./cython",'src/cpp'],
                    # # language="c++",
                    # extra_compile_args = ["-O3", "-Wall"]),
            # Extension("atpy.src.cython.structures",
                    # ['src/cython/structures.pyx'],
                    # include_dirs = [".",'../cython','src/cpp'],
                    # # language="c++",
                    # extra_compile_args = ["-O3", "-Wall"]),
            # Extension("atpy.src.cython.beamline",
                    # ['src/cython/beamline.pyx'],
                    # include_dirs = ["./cython",'src/cpp'],
                    # # language="c++",
                    # extra_compile_args = ["-O3", "-Wall"]),
            # Extension("atpy.src.cython.parser",
                    # ['src/cython/parser.pyx'],
                    # include_dirs = ["./cython"],
                    # # language="c++",
                    # extra_compile_args = ["-O3", "-Wall"]),
            # Extension("atpy.src.cython.lattice",
                    # ['src/cython/lattice.pyx'],
                    # include_dirs = ["./cython",'src/cpp'],
                    # # language="c++",
                    # extra_compile_args = ["-O3", "-Wall"]),
            # ]

setup(name='atpy',
      version='1.0.0',
      author='LiuTao',
      packages=find_packages(),
      cmdclass = {'build_ext': build_ext},
      ext_modules=cythonize(set_extension(),
                            annotate=True,
                            force=True,
                            compiler_directives={'profile':True, 'language_level':3,
                                                'boundscheck':False, 'cdivision':True, 
                                                'initializedcheck':False,'linetrace':True}),
      package_dir={'':'.'},
      # data_files=['*.pxd','*.h','*.py'],
      # include_package_data = False,
      package_data = {
            "": ["*.pxd","*.in"],
            "atpy.src":["*.pxd"],
            "atpy.src.beamline":["*.pxd","*.pyx"],
            "atpy.src.parser":["*.pxd","*.pyx"],
            "atpy.src.optimize":["*.pxd","include/*.h","src/*cpp","*.pyx"],
            "atpy.src.element":["*.pxd","*.pyx"],
            "atpy.src.tpsa":["*.pxd","include/*.h","src/*cpp","*.pyx"],
            },
      zip_safe = False,
      
      install_requires=[
            'cython',
            'cymem',
            'matplotlib',
            'numpy',
        ],
     )
print(find_packages())
