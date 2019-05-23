"""
_version

Version information for MODAPE

"""

__all__ = [
    'VERSION_MAJOR',
    'VERSION_MINOR',
    'VERSION_PATCH',
    'VERSION_STRING',
    '__version__',
]

VERSION_MAJOR = 0
VERSION_MINOR = 1
VERSION_PATCH = 5

VERSION_STRING = '%s.%s.%s' % (VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH)

__version__ = VERSION_STRING
