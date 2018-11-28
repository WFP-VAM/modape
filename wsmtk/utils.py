import numpy as np
import datetime
import requests
import gdal
import ctypes
import multiprocessing
import multiprocessing.pool
import array
from .whittaker import ws2d, ws2d_vc, ws2d_vc_asy

# assign xrange to range if py2
try:
    range = xrange
except NameError:
    pass

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

# adapted from https://stackoverflow.com/questions/17223301/python-multiprocessing-is-it-possible-to-have-a-pool-inside-of-a-pool/17229030#17229030

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class Pool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess


class DateHelper:

    def __init__(self,rawdates,rtres,stres,start=None,tshift = 0,nupdate=0):

        if start:

            stop = (fromjulian(rawdates[-1]) + datetime.timedelta(rtres)).strftime('%Y%j')

            tdiff = (fromjulian(stop) - fromjulian(rawdates[0])).days

            self.daily = [(fromjulian(rawdates[0]) + datetime.timedelta(x)).strftime('%Y%j') for x in range(tdiff+1)]

            self.target = [self.daily[x] for x in range(self.daily.index(start),len(self.daily),stres)]
            self.target = self.target[-nupdate:]


        else:

            yrmin = int(min([x[:4] for x in rawdates]))
            yrmax = int(max([x[:4] for x in rawdates]))

            daily_tmp = [y for x in range(yrmin,yrmax+1,1) for y in tvec(x,1)]

            stop = (fromjulian(rawdates[-1]) + datetime.timedelta(rtres)).strftime('%Y%j')

            self.daily = daily_tmp[daily_tmp.index(rawdates[0]):daily_tmp.index(stop)+1]

            if stres == 5:

                target_temp = [y for x in range(yrmin,yrmax+1,1) for y in pentvec(x)]

            elif stres == 10:

                target_temp = [y for x in range(yrmin,yrmax+1,1) for y in dekvec(x)]

            else:

                target_temp = [y for x in range(yrmin,yrmax+1,1) for y in tvec(x,stres)]

            target_temp.sort()

            for sd in self.daily:
                if sd in target_temp:
                    start_target = sd
                    del sd
                    break

            for sd in reversed(self.daily):
                if sd in target_temp:
                    stop_target = sd
                    del sd
                    break


            self.target = target_temp[target_temp.index(start_target):target_temp.index(stop_target)+1]
            self.target = self.target[-nupdate:]

    def getDV(self,nd):

        return(np.full(len(self.daily),nd,dtype='float32'))

    def getDIX(self):
        return([self.daily.index(x) for x in self.target])

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

def tvec(yr,step):
    '''Create MODIS-like date vector with given timestep.'''

    start = fromjulian('{}001'.format(yr)) + datetime.timedelta()
    tdiff = fromjulian('{}001'.format(yr+1)) - start
    tv = [(start + datetime.timedelta(x)).strftime('%Y%j') for x in range(0,tdiff.days,step)]
    return(tv)

def pentvec(yr):
    '''Create pentadal date vector for given year with fixed days.'''

    t = []

    for m in range(1,13):
        for d in ['03','08','13','18','23','28']:
            try:
                t.append(datetime.datetime.strptime('{}{:02d}{}'.format(yr,m,d),'%Y%m%d').date().strftime('%Y%j'))
            except ValueError:
                pass
    return(t)

def dekvec(yr):
    '''Create dekadal date vector for given year with fixed days.'''

    return([datetime.datetime.strptime(str(yr)+y+x,'%Y%m%d').date().strftime('%Y%j') for x in ['05','15','25'] for y in [str(z).zfill(2) for z in range(1,13)]])


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

    params = {}

    for key, value in kwargs.items():
        params[key] = value
    return params

def init_worker(shared_array_,parameters_):
    '''Initialize worker for smoothing.

    Args:
        shared_array_: Object returned by init_shared
        parameters_: Dictionary returned by init_parameters
    '''

    global arr_raw
    global arr_smooth
    global arr_sgrid
    global parameters

    arr_raw = tonumpyarray(shared_array_)
    parameters = parameters_
    arr_raw.shape = parameters['rdim']


    try:
         arr_sgrid = tonumpyarray(parameters['shared_array_sgrid'])
    except (KeyError, AttributeError):
         arr_sgrid = None

    try:
         arr_smooth = tonumpyarray(parameters['shared_array_smooth'])
         arr_smooth.shape = parameters['sdim']
    except (KeyError, AttributeError):
         arr_smooth = None


def execute_ws2d(ix):
    '''Execute whittaker smoother with fixed s in worker.

    Args:
        ix ([int]): List of indices as integer
    '''

    arr_raw[ix,:] = ws2d(y = arr_raw[ix,:], lmda = 10**parameters['s'], w = np.array((arr_raw[ix,:] != parameters['nd']) * 1, dtype='float32'))

    if parameters['shared_array_smooth']:

        z2 = parameters['vec_dly'].copy()
        z2[ z2 != parameters['nd'] ] = arr_raw[ix,:]
        z2[...] = ws2d(y = z2, lmda = 0.0001, w = np.array((z2 != parameters['nd']) * 1,dtype='float32'))
        arr_smooth[ix,:] = z2[parameters['dix']]

def execute_ws2d_sgrid(ix):
    '''Execute whittaker smoother with s from grid in worker.'''

    arr_raw[ix,:] = ws2d(y = arr_raw[ix,:], lmda = 10**arr_sgrid[ix], w = np.array((arr_raw[ix,:] != parameters['nd']) * 1,dtype='float32'))

    if parameters['shared_array_smooth']:

        z2 = parameters['vec_dly'].copy()
        z2[ z2 != parameters['nd'] ] = arr_raw[ix,:]
        z2[...] = ws2d(y = z2, lmda = 0.0001, w = np.array((z2 != parameters['nd']) * 1,dtype='float32'))
        arr_smooth[ix,:] = z2[parameters['dix']]

def execute_ws2d_vc(ix):
    '''Execute whittaker smoother with V-curve optimization of s in worker.'''

    if parameters['p']:

        arr_raw[ix,:], arr_sgrid[ix] = ws2d_vc_asy(y = arr_raw[ix,:], w = np.array((arr_raw[ix,:] != parameters['nd']) * 1,dtype='float32'), llas = array.array('f',parameters['srange']), p = parameters['p'])
    else:

        arr_raw[ix,:], arr_sgrid[ix] = ws2d_vc(y = arr_raw[ix,:], w = np.array((arr_raw[ix,:] != parameters['nd']) * 1,dtype='float32'), llas = array.array('f',parameters['srange']))

    if parameters['shared_array_smooth']:

        z2 = parameters['vec_dly'].copy()
        z2[ z2 != parameters['nd'] ] = arr_raw[ix,:]
        z2[...] = ws2d(y = z2, lmda = 0.0001, w = np.array((z2 != parameters['nd']) * 1,dtype='float32'))
        arr_smooth[ix,:] = z2[parameters['dix']]


def execute_ws2d_vcOpt(ix):
    '''Execute whittaker smoother with V-curve 2step-optimization of s in worker.'''

    z, lopt =  ws2d_vc(y = arr_raw[ix,:], w = np.array((arr_raw[ix,:] != parameters['nd']) * 1,dtype='float32'), llas = array.array('f',parameters['srange']))

    srange_lim = parameters['srange'][parameters['srange'] <= np.log10(lopt)]

    if len(srange_lim)==1:
        srange_lim = np.concatenate([srange_lim-0.2,srange_lim])

    if parameters['p']:

        arr_raw[ix,:], arr_sgrid[ix] = ws2d_vc_asy(y = arr_raw[ix,:], w = np.array((arr_raw[ix,:] != parameters['nd']) * 1,dtype='float32'), llas = array.array('f',srange_lim), p = parameters['p'])
    else:

        arr_raw[ix,:], arr_sgrid[ix] = ws2d_vc(y = arr_raw[ix,:], w = np.array((arr_raw[ix,:] != parameters['nd']) * 1,dtype='float32'), llas = array.array('f',srange_lim))

    if parameters['shared_array_smooth']:

        z2 = parameters['vec_dly'].copy()
        z2[ z2 != parameters['nd'] ] = arr_raw[ix,:]
        z2[...] = ws2d(y = z2, lmda = 0.0001, w = np.array((z2 != parameters['nd']) * 1,dtype='float32'))
        arr_smooth[ix,:] = z2[parameters['dix']]
