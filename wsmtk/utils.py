from numpy.lib.stride_tricks import as_strided as ast
import numpy as np
import datetime
import requests
import gdal
import ctypes
import multiprocessing

def dtype_GDNP(dt):
    '''GDAL/NP DataType helper'''
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


def block_view(A, block= (3, 3)):
    ## Credit to http://stackoverflow.com/a/5078155/1828289
    """Provide a 2D block view to 2D array. No error checking made.
    Therefore meaningful (as implemented) only for blocks strictly
    compatible with the shape of A."""
    shape= (A.shape[0]// block[0], A.shape[1]// block[1])+ block
    strides= (block[0]* A.strides[0], block[1]* A.strides[1])+ A.strides
    return ast(A, shape= shape, strides= strides)

def LDOM(x):
    '''get last day of month'''
    yr = x.year
    mn = x.month
    if mn is 12:
        mn = 1
        yr +=1
    else:
        mn += 1
    return(datetime.date(yr,mn,1) - datetime.timedelta(days=1))


def aoi2ix(ref,aoi,res):

    '''Extract indices for intersection over reference'''

    isect = [max(ref[0],aoi[0]),min(ref[1],aoi[1]),min(ref[2],aoi[2]),max(ref[3],aoi[3])]

    xoff = int(round((isect[0] - ref[0])/res))
    yoff = int(round((ref[1] - isect[1])/res))

    xd = int(round((isect[2] - isect[0])/res))#+1
    yd = int(round((isect[1] - isect[3])/res))#+1

    return((xoff,xd,yoff,yd))


## from https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python

# overriding requests.Session.rebuild_auth to mantain headers when redirected

class SessionWithHeaderRedirection(requests.Session):

    AUTH_HOST = 'urs.earthdata.nasa.gov'

    def __init__(self, username, password):

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

def txx(x):
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
    return datetime.datetime.strptime(x,'%Y%j').date()


def init_array(cshp, ncell):
    '''Create shared value array for smoothing'''
    shared_array_base = multiprocessing.Array(ctypes.c_float,ncell,lock=False)
    main_nparray = np.frombuffer(shared_array_base, dtype=ctypes.c_float)
    main_nparray = main_nparray.reshape(cshp[0]*cshp[1],cshp[2])

    assert main_nparray.base.base is shared_array_base
    return(main_nparray)

def init_lamarray(cshp):
    '''Create shared lambda array for smoothing'''
    nc = cshp[0] * cshp[1]
    shared_array_base = multiprocessing.Array(ctypes.c_float,nc,lock=False)
    main_nparray = np.frombuffer(shared_array_base, dtype=ctypes.c_float)

    assert main_nparray.base is shared_array_base
    return(main_nparray)

def worker_fn(ix):
    #worker function for parallel smoothing
    if (wts[ix,...].sum().item() != 0.0):
        arr[ix,...] = ws2d(arr[ix,...],lmda = l,w = wts[ix,...])
