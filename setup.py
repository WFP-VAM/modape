#!/usr/bin/env python
# pylint: disable=invalid-name
"""setup.py for MODAPE"""

from setuptools import setup, Extension, find_packages

import numpy
import _version
USE_CYTHON = 'auto'

if USE_CYTHON:
    try:
        from Cython.Distutils import build_ext
        if USE_CYTHON == 'auto':
            USE_CYTHON = True
    except ImportError:
        if USE_CYTHON == 'auto':
            USE_CYTHON = False
        else:
            raise

cmdclass = {}
ext_modules = []

if USE_CYTHON:
    ext_modules += [
        Extension("modape.whittaker",
                  ["modape/_whittaker.pyx"], extra_compile_args=["-O3", "-ffast-math"])]
    cmdclass.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("modape.whittaker",
                  ["modape/_whittaker.c"], extra_compile_args=["-O3", "-ffast-math"])]

setup(
    name='modape',
    description='MODIS Assimilation and Processing Engine',
    version=_version.__version__,
    author='Valentin Pesendorfer',
    author_email='valentin.pesendorfer@wfp.org',
    url='http://wfp-vam.github.io/modape',
    long_description='''HDF5 based processing chain optimized for MODIS data combined with a State-of-the art whittaker smoother, implemented as fast C-extension through Cython and including a V-curve optimization of the smoothing parameter.\n\nFor more information, please visit: http://github.com/WFP-VAM/modape''',
    include_dirs=[numpy.get_include()],
    entry_points={
        'console_scripts':[
            'modis_download=modape.scripts.modis_download:main',
            'modis_collect=modape.scripts.modis_collect:main',
            'modis_smooth=modape.scripts.modis_smooth:main',
            'modis_window=modape.scripts.modis_window:main',
            'modis_info=modape.scripts.modis_info:main',
            'modis_product_table=modape.scripts.modis_product_table:main',
            'csv_smooth=modape.scripts.csv_smooth:main',
            'rts_smooth=modape.scripts.rts_smooth:main',
        ]
    },
    packages=find_packages(),
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    include_package_data=True,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    install_requires=[
        'numpy>=1.15.1',
        'gdal>=2',
        'h5py>=2.9',
        'beautifulsoup4>=4.7',
        'requests>=2',
        'progress>=1.5',
        'pandas>=0.24',
        'cryptography>=2.6',
        'mock;python_version<"3.0"'
    ],
    python_requires='>=2.7.11, <4',
)
