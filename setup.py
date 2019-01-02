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
        Extension("wsmtk.whittaker", ["wsmtk/_whittaker.pyx"],extra_compile_args = ["-O3", "-ffast-math"])]
    cmdclass.update({'build_ext': build_ext })
else:
    ext_modules += [
        Extension("wsmtk.whittaker", ["wsmtk/_whittaker.c"],extra_compile_args = ["-O3", "-ffast-math"])]


setup(
    name='wsmtk',
    description='Whittaker Smoothing Toolkit',
    version=_version.__version__,
    author='Valentin Pesendorfer',
    author_email='valentin.pesendorfer@wfp.org',
    long_description=open('README.rst').read(),
    include_dirs=[numpy.get_include()],
    entry_points={
    'console_scripts':[
    'downloadMODIS=wsmtk.downloadMODIS:main',
    'processMODIS=wsmtk.processMODIS:main',
    'windowMODIS=wsmtk.windowMODIS:main',
    'smoothMODIS=wsmtk.smoothMODIS:main',
    'infoMODIS=wsmtk.infoMODIS:main',
    'producttableMODIS=wsmtk.producttableMODIS:main',
    'smoothCSV=wsmtk.smoothCSV:main',
    'smoothRTS=wsmtk.smoothRTS:main',
    ]
    },
    packages=find_packages(),
    cmdclass = cmdclass,
    ext_modules=ext_modules,
    include_package_data=True,
    classifiers=[
    'Development Status :: 3 - Alpha',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
],
    install_requires=['numpy','gdal>=2','h5py','beautifulsoup4','requests','progress','pandas', "cryptography"],
    python_requires='>=2.7.11, <4',
)
