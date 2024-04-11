from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

richdem_include_path = "/path/to/richdem/headers"  # Modify this path accordingly

ext_modules = [
    Extension(
        "pyas",
        sources=["pyas.pyx", "pyasc.cpp"],
        include_dirs=[np.get_include(), richdem_include_path],
        extra_compile_args=["-std=c++11"] if any(f.endswith('.cpp') for f in ["pyas.pyx", "pyasc.cpp"]) else [],
        language="c++"  # Specify the language for the extension
    )
]

setup(
    ext_modules=cythonize(ext_modules)
)
