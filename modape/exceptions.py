'''Custom exceptions for MODAPE'''

from typing import List, Tuple

class DownloadError(Exception):
    '''Exception for failed download of MODIS data'''

    def __init__(self, fails: List[Tuple]) -> None:
        """Init custom DownloadError Exception.

        Args:
            fails (List[Tuple]): List of failed downloads. Each is tuple with (URI, Error).
        """


        message = '''

ERROR downloading MODIS data!

Failed downloads:

'''
        for failed, error in fails:
            message += f"{failed}: {error}\n"

        super(DownloadError, self).__init__(message)
