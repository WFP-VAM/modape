class DownloadError(Exception):
    def __init__(self, fails):

        message = '''

ERROR downloading MODIS data!

Failed downloads:

'''
        for failed, error in fails:
            message += f"{failed}: {error}\n"

        super(DownloadError, self).__init__(message)
