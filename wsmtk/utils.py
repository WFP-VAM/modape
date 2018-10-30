import numpy as np
import datetime
import requests
import gdal
import ctypes
import multiprocessing
from .whittaker import ws2d, ws2d_vc, ws2d_vc_asy

def dtype_GDNP(dt):
    '''GDAL/NP DataType helper.

    Parses datatype in str or GDAL int.

    Agrs:
        dt (str/int): DataType

    Returns:
        Tuple with DataType as GDAL integer and string
    '''

    dt_dict={
    1: 'uint8',
    2: 'uint16',
    3: 'int16',
    4: 'uint32',
    5: 'int32',
    6: 'float32',
    7: 'float64'
    }

    dt_tuple = [(k,v) for k,v in dt_dict.items() if k == dt or v == dt]
    return(dt_tuple[0])

def LDOM(x):
    '''Get last day of month.

    Args:
        x (datetime): Input date

    Returns:
        Datetime object
    '''


    yr = x.year
    mn = x.month
    if mn is 12:
        mn = 1
        yr +=1
    else:
        mn += 1
    return(datetime.date(yr,mn,1) - datetime.timedelta(days=1))

class SessionWithHeaderRedirection(requests.Session):
    ''' Session class for MODIS query.

    Overriding requests.Session.rebuild_auth to mantain headers when redirected.
    From https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
    '''

    AUTH_HOST = 'urs.earthdata.nasa.gov'

    def __init__(self, username, password):
        '''Create SessionWithHeaderRedirection instance.

        Args:
            username (str): Earthdata username
            password (str): Earthdata password
        '''

        super(SessionWithHeaderRedirection,self).__init__()

        self.auth = (username, password)

   # Overrides from the library to keep headers when redirected to or from

   # the NASA auth host.

    def rebuild_auth(self, prepared_request, response):

        headers = prepared_request.headers

        url = prepared_request.url

        if 'Authorization' in headers:

            original_parsed = requests.utils.urlparse(response.request.url)

            redirect_parsed = requests.utils.urlparse(url)

            if (original_parsed.hostname != redirect_parsed.hostname) and redirect_parsed.hostname != self.AUTH_HOST and original_parsed.hostname != self.AUTH_HOST:
                del headers['Authorization']

        return

class FileHandler:
    ''' Filehandler class to handle GDAL file references'''

    def __init__(self, files, sds):
        ''' Create handler instance

        Args:
            files (list): List of total file paths
            sds (str): Subdataset to extract
        '''

        self.files = files
        self.sds = sds

    def open(self):
        '''Open the files and store handle'''

        self.handles = []

        for ix, f in enumerate(self.files):

            # try statement to catch problem with reading file

            try:

                fl_o = gdal.Open(f)

                val_sds = [x[0] for x in fl_o.GetSubDatasets() if self.sds in x[0]][0]

                self.handles.append(gdal.Open(val_sds))

                fl_o = None


            except AttributeError:

                self.handles.append(None)

    def close(self):
        '''Open the files'''

        for ii in range(len(self.handles)):
            self.handles[ii] = None


def txx(x):
    '''Converts tempint integer to flag.'''

    if x:
        if int(x) == 5:
            return 'p'
        elif int(x) == 10:
            return 'd'
        else:
            return 'c'
    else:
        return 'n'


def fromjulian(x):
    '''Parses julian date string to datetime object.'''

    return datetime.datetime.strptime(x,'%Y%j').date()

def init_shared(ncell):
    '''Create shared value array for smoothing.'''
    shared_array_base = multiprocessing.Array(ctypes.c_float,ncell,lock=False)
    return(shared_array_base)

def tonumpyarray(shared_array):
    '''Create numpy array from shared memory.'''
    nparray= np.frombuffer(shared_array,dtype=ctypes.c_float)
    assert nparray.base is shared_array
    return nparray

def init_parameters(**kwargs):
    '''Initialize parameters for smoothing in workers.'''

    params = dict(zip(['shared_array','shared_sarr','nd','s','srange','p','dim'],[None for x in range(7)]))

    for key, value in kwargs.items():
        params[key] = value
    return params

def init_worker(shared_array_,parameters_):
    '''Initialize worker for smoothing.

    Args:
        shared_array_: Object returned by init_shared
        parameters_: Dictionary returned by init_parameters
    '''

    global shared_array
    global nd
    global s
    global srange
    global p
    global dim
    shared_array = tonumpyarray(shared_array_)
    nd = parameters_['nd']
    s = parameters_['s']
    srange = parameters_['srange']
    p = parameters_['p']
    dim = parameters_['dim']

    try:
        global shared_sarr
        shared_sarr = tonumpyarray(parameters_['shared_sarr'])
    except:
        shared_sarr = None


def execute_ws2d(ix):
    '''Execute whittaker smoother with fixed s in worker.

    Args:
        ix ([int]): List of indices as integer
    '''
    arr = tonumpyarray(shared_array)
    arr.shape = dim
    for ii in ix:
        if (arr[ii,] != nd ).any():
            arr[ii,] = ws2d(y = arr[ii,], lmda = s, w = np.array((arr[ii,] != nd) * 1,dtype='float32'))

def execute_ws2d_sgrid(ix):
    '''Execute whittaker smoother with s from grid in worker.'''
    arr = tonumpyarray(shared_array)
    sarr = tonumpyarray(shared_sarr)
    arr.shape = dim

    for ii in ix:
        if (arr[ii,] != nd ).any():
            arr[ii,] = ws2d(y = arr[ii,], lmda = 10**sarr[ii], w = np.array((arr[ii,] != nd ) * 1,dtype='float32'))

def execute_ws2d_vc(ix):
    '''Execute whittaker smoother with V-curve optimization of s in worker.'''
    arr = tonumpyarray(shared_array)
    sarr = tonumpyarray(shared_sarr)
    arr.shape = dim

    for ii in ix:
        if (arr[ii,] != nd ).any():
            arr[ii,], sarr[ii] =  ws2d_vc(y = arr[ii], w = np.array((arr[ii,] != nd) * 1,dtype='float32'), llas = srange)

def execute_ws2d_vc_asy(ix):
    '''Execute asymmetric whittaker smoother with V-curve optimization of s in worker.'''
    arr = tonumpyarray(shared_array)
    sarr = tonumpyarray(shared_sarr)
    arr.shape = dim

    for ii in ix:
        if (arr[ii,] != nd ).any():
            arr[ii,], sarr[ii] = ws2d_vc_asy(y = arr[ii], w = np.array((arr[ii,] != nd) * 1,dtype='float32'), llas = srange, p = p)
