"""Custom exceptions for MODAPE"""
#pylint: disable=W0107
from typing import List

class DownloadError(Exception):
    """Exception for failed download of MODIS data"""

    def __init__(self, fails: List) -> None:
        """Init custom DownloadError Exception.

        Args:
            fails (List): List of failed file IDs.
        """


        message = """

ERROR downloading MODIS data!

Failed downloads:

"""
        for failed in fails:
            message += f"{failed}\n"

        super(DownloadError, self).__init__(message)


class TargetNotEmptyError(Exception):
    """Exception when target directory is not empty
    when it's required to be."""
    pass

class HDF5CreationError(Exception):
    """Exception when creating HDF5 file"""
    pass

class HDF5WriteError(Exception):
    """Exception when writing to HDF5 file"""
    pass

class SgridNotInitializedError(Exception):
    """Exception when requiring non existing sgrid"""
    pass

class HDF5MosaicError(Exception):
    """Exception when mosaicing incompatible HDF5s"""
