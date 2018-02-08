#!/usr/bin/env python

#from distutils.core import setup
#from distutils.extension import Extension
from setuptools import setup, Extension, find_packages
import numpy
import _version
USE_CYTHON = 'auto'


if USE_CYTHON:
    try:
        from Cython.Distutils import build_ext
        if USE_CYTHON=='auto':
            USE_CYTHON=True
    except ImportError:
        if USE_CYTHON=='auto':
            USE_CYTHON=False
        else:
            raise

cmdclass = {}
ext_modules = []

if USE_CYTHON:
    ext_modules += [
        Extension("wsmtk.whittaker", ["wsmtk/whittaker/_whittaker.pyx"],extra_compile_args = ["-O3", "-ffast-math"])]
    cmdclass.update({'build_ext': build_ext })
else:
    ext_modules += [
        Extension("wsmtk.whittaker", ["wsmtk/whittaker/_whittaker.c"],extra_compile_args = ["-O3", "-ffast-math"])]


setup(
    name='wsmkt',
    description='Whittaker Smoothing Toolkit',
    version=_version.__version__,
    author='Valentin Pesendorfer',
    author_email='valentin.pesendorfer@wfp.org',
    long_description=open('README.rst').read(),
    include_dirs=[numpy.get_include()],
    #packages=['wsmtk'],
    packages=find_packages(),
    cmdclass = cmdclass,
    ext_modules=ext_modules,
    classifiers=[
    # How mature is this project? Common values are
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
    'Development Status :: 3 - Alpha',

    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3',
],
    install_requires=['numpy'],
)
