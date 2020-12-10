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
include_dirs=[".",'atpy/src']

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

#Setup config
setup(name='atpy',
      version='1.0.0',
      author='LiuTao',
      packages=find_packages(),
      cmdclass = {'build_ext': build_ext},
      ext_modules=cythonize(set_extension(),
                            annotate=False,
                            force=True,
                            compiler_directives={'profile':True, 'language_level':3,
                                                'boundscheck':False, 'cdivision':True, 
                                                'initializedcheck':False,'linetrace':True}),
      package_dir={'':'.'},
      package_data = {
            "atpy": ["*.pxd","*.py"],
            "atpy.src":["*.pxd"],
            "atpy.src.beamline":["*.pxd","*.pyx"],
            "atpy.src.parser":["*.pxd","*.pyx"],
            "atpy.src.optimize":["*.pxd","include/*.h","src/*cpp","*.pyx"],
            "atpy.src.element":["*.pxd","*.pyx"],
            "atpy.src.tpsa":["*.pxd","include/*.h","src/*cpp","*.pyx"],
            },
      zip_safe = False,
      
      install_requires=[
            'cython>=3.0.0',
            'cymem>=2.0.2',
            'matplotlib>=3.1.0',
            'numpy>=1.19.1',
        ],
     )
