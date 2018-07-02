from numpy.lib.stride_tricks import as_strided as ast
import numpy as np
import datetime
import requests
import gdal
import ctypes
import multiprocessing
from .whittaker import ws2d, ws2d_vc, ws2d_vc_asy

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

def init_shared(ncell):
    '''Create shared value array for smoothing'''
    shared_array_base = multiprocessing.Array(ctypes.c_float,ncell,lock=False)
    return(shared_array_base)

def tonumpyarray(shared_array):

    nparray= np.frombuffer(shared_array,dtype=ctypes.c_float)
    assert nparray.base is shared_array
    return nparray

def init_worker(shared_arr_,nd_,l_ = None,llas_ = None,p_ = None):
    global shared_arr
    global nd
    global l
    global llas
    global p
    shared_arr = tonumpyarray(shared_arr_)
    nd = nd_
    l = l_
    llas = llas_
    p = p_

def execute_ws2d(ix):
    #worker function for parallel smoothing using whittaker 2d with fixed lambda
    arr = tonumpyarray(shared_array)
    if (arr[ix] != nd ).any():
        arr[ix] = ws2d(y = arr[ix], lmda = l, w = np.array((arr[ix] != nd) * 1,dtype='float32'))

def execute_ws2d_lgrid(ix):
    #worker function for parallel smoothing using whittaker 2d with existing lambda grid
    if (arr[ix] != nd ).any():
        arr[ix] = ws2d(y = arr[ix], lmda = 10**lamarr[ix], w = np.array((arr[ix] != nd ) * 1,dtype='float32'))

def execute_ws2d_vc(ix):
    #worker function for parallel smoothing using whittaker 2d with v-curve optimization
    if (arr[ix] != nd ).any():
        arr[ix], lamarr[ix,] =  ws2d_vc(y = arr[ix], w = np.array((arr[ix] != nd ) * 1,dtype='float32'), llas = llas)

def execute_ws2d_vc_asy(ix):
    #worker function for parallel asymmetric smoothing using whittaker 2d with v-curve optimization
    if (arr[ix] != nd ).any():
        arr[ix], lamarr[ix] = ws2d_vc_asy(y = arr[ix], w = np.array((arr[ix] != nd ) * 1,dtype='float32'), llas = llas, p = p)
