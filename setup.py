#!/usr/bin/env python
# pylint: disable=invalid-name,E0401
"""Minimal setup.py for MODAPE - handles only Cython extensions"""

from setuptools import setup, Extension
import numpy

USE_CYTHON = "auto"

if USE_CYTHON:
    try:
        from Cython.Distutils import build_ext
        if USE_CYTHON == "auto":
            USE_CYTHON = True
    except ImportError:
        if USE_CYTHON == "auto":
            USE_CYTHON = False
        else:
            raise

cmdclass = {}
ext_modules = []

if USE_CYTHON:
    ext_modules += [
        Extension("modape.whittaker",
                  ["modape/_whittaker.pyx"], 
                  extra_compile_args=["-O3", "-ffast-math"],
                  include_dirs=[numpy.get_include()])]
    cmdclass.update({"build_ext": build_ext})
else:
    ext_modules += [
        Extension("modape.whittaker",
                  ["modape/_whittaker.c"], 
                  extra_compile_args=["-O3", "-ffast-math"],
                  include_dirs=[numpy.get_include()])]

# Set Python 3 for Cython
for ext in ext_modules:
    ext.cython_directives = {"language_level": "3"}

# Only handle extensions - all other config in pyproject.toml
setup(
    cmdclass=cmdclass,
    ext_modules=ext_modules,
)