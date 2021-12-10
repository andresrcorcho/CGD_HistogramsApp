
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy

extensions = [
    Extension(
                "Thread_C",
                sources=["Cfunctions.pyx"],
                libraries=["m"],
                include_dirs=[numpy.get_include()],
                extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-openmp" ],
              extra_link_args=['-openmp']
            )]

setup(
    #ext_modules = cythonize('Cfunctions.pyx'),
    cmdclass = {"build_ext": build_ext},
    ext_modules = extensions
    
    )
