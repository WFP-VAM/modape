#!/usr/bin/env python
# pylint: disable=invalid-name,E0401
"""setup.py for MODAPE"""

from setuptools import setup, Extension, find_packages

import numpy
import _version
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
                  ["modape/_whittaker.pyx"], extra_compile_args=["-O3", "-ffast-math"])]
    cmdclass.update({"build_ext": build_ext})
else:
    ext_modules += [
        Extension("modape.whittaker",
                  ["modape/_whittaker.c"], extra_compile_args=["-O3", "-ffast-math"])]

# set py 3:

for ext in ext_modules:
    ext.cython_directives = {"language_level": "3"}


setup(
    name="modape",
    description="MODIS Assimilation and Processing Engine",
    version=_version.__version__,
    author="Valentin Pesendorfer",
    author_email="valentin.pesendorfer@wfp.org",
    url="http://wfp-vam.github.io/modape",
    long_description="""HDF5 based processing chain optimized for MODIS data combined with a State-of-the art whittaker smoother, implemented as fast C-extension through Cython and including a V-curve optimization of the smoothing parameter.\n\nFor more information, please visit: http://github.com/WFP-VAM/modape""",
    include_dirs=[numpy.get_include()],
    entry_points={
        "console_scripts":[
            "modis_download=modape.scripts.modis_download:cli_wrap",
            "modis_collect=modape.scripts.modis_collect:cli_wrap",
            "modis_smooth=modape.scripts.modis_smooth:cli_wrap",
            "modis_window=modape.scripts.modis_window:cli_wrap",
            "csv_smooth=modape.scripts.csv_smooth:cli_wrap",
            "modis_info=modape.scripts.modis_info:cli_wrap",
            "modape_version=_version.version_info:main",
        ]
    },
    packages=find_packages(),
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
    ],
    setup_requires=["pytest-runner"],
    tests_require=["pytest", "click"],
    install_requires=[
        "click<8",
        "numpy>=1.16.1",
        "gdal>=2",
        "h5py>=2.9",
        "python-cmr>=0.4",
        "requests>=2",
        "pandas>=0.24",
        "pycksum>=0.4.3",
    ],
    python_requires=">=3, <4",
)
