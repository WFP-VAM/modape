# pylint: disable=invalid-name, global-variable-undefined, broad-except
"""
Utility classes and functions used in MODIS processing chain.

Author: Valentin Pesendorfer, April 2019
"""
#pylint: disable=E0401
import datetime
from typing import List

import numpy as np
import requests
from base64 import b64encode


__all__ = [
    "SessionWithHeaderRedirection",
    "DateHelper",
    "fromjulian",
    "tvec",
    "pentvec",
    "dekvec",
    ]

class SessionWithHeaderRedirection(requests.Session):
    """Session class for MODIS query.

    Overriding requests.Session.rebuild_auth to mantain headers when redirected.
    From https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
    """

    AUTH_HOST = "urs.earthdata.nasa.gov"

    def __init__(self, username, password):
        """Create SessionWithHeaderRedirection instance.

        Args:
            username: Earthdata username
            password: Earthdata password
        """

        super(SessionWithHeaderRedirection, self).__init__()
        self.__username = username
        self.__password = password

    def rebuild_auth(self, prepared_request, response):
        """Overrides from the library to keep headers when redirected
           to or from the NASA auth host.
        """
        headers = prepared_request.headers
        original_parsed = requests.utils.urlparse(response.request.url)
        redirect_parsed = requests.utils.urlparse(prepared_request.url)
        if (original_parsed.hostname != redirect_parsed.hostname) and \
            redirect_parsed.hostname != self.AUTH_HOST and \
            original_parsed.hostname != self.AUTH_HOST:
            if 'Authorization' in headers:
                del headers['Authorization']
        elif redirect_parsed.hostname == self.AUTH_HOST and 'Authorization' not in headers:
            headers['Authorization'] = 'Basic %s' %  \
                b64encode(bytes(f"{self.__username}:{self.__password}", 'utf-8')).decode("ascii")
        return


class DateHelper(object):
    """Helper class for handling dates in temporal interpolation."""

    def __init__(self, rawdates, rtres, stres, start=None):
        """Creates the date lists from input.

        Args:
             rawdates: list of dates from raw file(s)
             rtres: raw temporal resolution
             stres: smooth temporal resolution
             start: start date for custom interpolation
             nupdate: number of points in time to be updated in file (backwards)
            """

        if start:
            stop = (fromjulian(rawdates[-1]) + datetime.timedelta(rtres)).strftime("%Y%j")
            tdiff = (fromjulian(stop) - fromjulian(rawdates[0])).days
            self.daily = [(fromjulian(rawdates[0]) + datetime.timedelta(x)).strftime("%Y%j") for x in range(tdiff+1)]
            self.target = [self.daily[x] for x in range(self.daily.index(start), len(self.daily), stres)]
        else:
            yrmin = int(min([x[:4] for x in rawdates]))
            yrmax = int(max([x[:4] for x in rawdates]))
            daily_tmp = [y for x in range(yrmin, yrmax+2, 1) for y in tvec(x, 1)]
            stop = (fromjulian(rawdates[-1]) + datetime.timedelta(rtres)).strftime("%Y%j")
            self.daily = daily_tmp[daily_tmp.index(rawdates[0]):daily_tmp.index(stop)+1]

            if stres != rtres:

                if stres == 5:
                    target_temp = [y for x in range(yrmin, yrmax+1, 1) for y in pentvec(x)]
                elif stres == 10:
                    target_temp = [y for x in range(yrmin, yrmax+1, 1) for y in dekvec(x)]
                else:
                    target_temp = [y for x in range(yrmin, yrmax+1, 1) for y in tvec(x, stres)]
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

            else:
                self.target = rawdates
        self.target_length = len(self.target)

    def getDV(self, nd):
        """Gets an array of no-data values in daily timesteps.

        Args:
            nd: no-data value

        Returns:
            numpy array with no-data values in daily steps
        """

        return np.full(len(self.daily), nd, dtype="double")

    def getDIX(self, nupdate=0):
        """Gets indices of target dates in daily no-data array.

        Returns:
            list with indices of target dates in no-data array
        """

        return [self.daily.index(x) for x in self.target[-nupdate:]]

def check_sequential(
        reference: List[str],
        update: List[str]
) -> bool:
    """Helper to check if update continues sequence.

    Helps to check if updated dates are all in sequence.

    Args:
        reference (List[str]): Reference list of dates to be updated.
        update (List[str]): New dates to update.

    Returns:
        bool: True if sequential, False if not.

    """

    index = [i for i, item in enumerate(reference) if item in set(update)]
    index_it = iter(index)
    first_item = next(index_it)
    is_sequential = all(a == b for a, b in enumerate(index_it, first_item + 1))

    return is_sequential

def fromjulian(x):
    """Parses julian date string to datetime object.

    Args:
        x: julian date as string YYYYJJJ

    Returns:
        datetime object parsed from julian date
    """

    return datetime.datetime.strptime(x, "%Y%j").date()

def tvec(yr, step):
    """Create MODIS-like date vector with given timestep.

    Args:
        yr: year
        step: timestep

    Returns:
        list with dates
    """

    start = fromjulian("{}001".format(yr)) + datetime.timedelta()
    tdiff = fromjulian("{}001".format(yr+1)) - start
    tv = [(start + datetime.timedelta(x)).strftime("%Y%j") for x in range(0, tdiff.days, step)]
    return tv

def pentvec(yr):
    """Create pentadal date vector for given year with fixed days.

    Args:
        yr: year

    Returns:
        list of dates
    """

    t = []
    for m in range(1, 13):
        for d in ["03", "08", "13", "18", "23", "28"]:
            try:
                t.append(datetime.datetime.strptime("{}{:02d}{}".format(yr, m, d), "%Y%m%d").date().strftime("%Y%j"))
            except ValueError:
                pass
    return t

def dekvec(yr):
    """Create dekadal date vector for given year with fixed days.

    Args:
        yr: year

    Returns:
        list of dates
    """

    return([
        datetime.datetime.strptime(str(yr)+y+x, "%Y%m%d").date().strftime("%Y%j")
        for x in ["05", "15", "25"] for y in [str(z).zfill(2)
                                              for z in range(1, 13)]
    ])
