"""
MODIS query & download class.

This file contains the class used for querying and downloading
raw MODIS products from NASA's servers.

Author: Valentin Pesendorfer, April 2019
"""
from __future__ import absolute_import, division, print_function

from datetime import datetime
import os
from os.path import basename
try:
    from pathlib2 import Path
except ImportError:
    from pathlib import Path
import re
from subprocess import Popen
import sys
import time
import uuid
import warnings

import numpy as np
import requests
from bs4 import BeautifulSoup # pylint: disable=import-error
from modape.utils import fromjulian

# turn off BeautifulSoup warnings
warnings.filterwarnings('ignore', category=UserWarning, module='bs4')

class ModisQuery(object):
    """Class for querying and downloading MODIS data."""

    def __init__(self, url, begindate,
                 enddate, username=None, password=None,
                 targetdir=os.getcwd(), global_flag=None,
                 tile_filter=None, strict_begindate=False):
        """Creates a ModisQuery object.

        Args:
            url: Query URL as created by modis_download.py
            begindate: Begin of period for query as ISO 8601 date string (YYYY-MM-DD)
            enddate: End of period for query as ISO 8601 date string (YYYY-MM-DD)
            username: Earthdata username (only required for download)
            password: Earthdata password (only required for download)
            targetdir: Path to target directory for downloaded files (default cwd)
            global_flag: Boolean flag indictaing queried product is global file instead of tiled product
            tile_filter: List of MODIS files to query and optionally download
        """

        self.query_url = url
        self.username = username
        self.password = password
        self.targetdir = Path(targetdir)
        self.files = []
        self.modis_urls = []
        self.begin = datetime.strptime(begindate, '%Y-%m-%d').date()
        self.end = datetime.strptime(enddate, '%Y-%m-%d').date()
        self.global_flag = global_flag

        # query for products using session object
        with requests.Session() as sess:

            print('Checking for MODIS products ...', end='')
            try:
                response = sess.get(self.query_url)
                self.statuscode = response.status_code
                response.raise_for_status()
            except requests.exceptions.RequestException as e:
                print(e)
                sys.exit(1)

            soup = BeautifulSoup(response.content, features='html.parser')

            # results for global products are date directories on server, so need to query separately
            date_regexp = re.compile(r'.+A(\d{7}).+')
            if self.global_flag:
                regex = re.compile('.*.hdf$')
                dates = np.array([x.getText() for x in soup.findAll('a', href=True) if re.match(r'\d{4}\.\d{2}\.\d{2}', x.getText())])
                dates_parsed = [datetime.strptime(x, '%Y.%m.%d/').date() for x in dates]
                dates_ix = np.flatnonzero(np.array([self.begin <= x < self.end for x in dates_parsed]))

                for date_sel in dates[dates_ix]:
                    try:
                        response = sess.get(self.query_url + date_sel)
                    except requests.exceptions.RequestException as e:
                        print(e)
                        print('Error accessing {} - skipping.'.format(self.query_url + date_sel))

                    soup = BeautifulSoup(response.content, features='html.parser')
                    hrefs = soup.find_all('a', href=True)
                    hdf_file = [x.getText() for x in hrefs if re.match(regex, x.getText())]

                    # Issue warning if there are multiple HDF
                    if len(hdf_file) > 1:
                        warnings.warn("More than 1 HDF files for specific date found for URL: {}".format(self.query_url + date_sel), Warning)

                    try:
                        fname = hdf_file[0]
                        fname = fname[fname.rfind('/') + 1:]
                        if (strict_begindate) and fromjulian(re.findall(date_regexp, fname)[0]) < self.begin:
                            continue
                        self.modis_urls.append(self.query_url + date_sel + hdf_file[0])
                    except IndexError:
                        print('No HDF file found in {} - skipping.'.format(self.query_url + date_sel))
                        continue
            else:
                regex = re.compile(r'.+(h\d+v\d+).+')
                urls = [x.getText() for x in soup.find_all('url')]

                if strict_begindate:
                    # exclude time steps that temporally covers the begin/end interval but are time-stamped before
                    urls = [x for x in urls if
                            fromjulian(re.findall(date_regexp, x[x.rfind('/') + 1:])[0]) >= self.begin]
                if tile_filter:
                    tiles = [x.lower() for x in tile_filter]
                    self.modis_urls = [x for x in urls if any(t in x for t in tiles)]
                else:
                    self.modis_urls = urls

                self.tiles = list({regex.search(x).group(1) for x in self.modis_urls})

        self.results = len(self.modis_urls)
        print('... done.\n')

        # check for results
        if self.results > 0:
            print('{} results found.\n'.format(self.results))
        else:
            print('0 results found. Please check query!')

    def set_credentials(self, username, password):
        """Set Earthdata credentials.

        Sets Earthdata username and password in created ModisQuery object.

        Args:
            username (str): Earthdata username
            password (str): Earthdata password
        """

        self.username = username
        self.password = password

    def download(self):
        """Downloads MODIS products.

        Download of files found through query, Earthdata username and password required!
        """

        if self.username is None or self.password is None:
            raise ValueError('No credentials found. Please run .setCredentials(username,password)!')

        print('[{}]: Downloading products to {} ...\n'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()), self.targetdir))

        # if targetdir doesn't exist, create
        self.targetdir.mkdir(parents=True, exist_ok=True)

        flist = self.targetdir.joinpath(str(uuid.uuid4())).as_posix()

        # write URLs of query resuls to disk for WGET
        with open(flist, 'w') as thefile:
            for item in self.modis_urls:
                thefile.write('%s\n' % item)

        args = [
            'aria2c',
            '--file-allocation=none',
            '-m', '50',
            '--retry-wait', '2',
            '-c',
            '-x', '10',
            '-s', '10',
            '--http-user', self.username,
            '--http-passwd', self.password,
            '-d', self.targetdir.as_posix(),
        ]

        # execute subprocess
        process_output = Popen(args + ['-i', flist])
        process_output.wait()

        # remove filelist.txt if all downloads are successful
        if process_output.returncode != 0:
            print('\nError (error code {}) occured during download, please check files against MODIS URL list ({})!\n'.format(process_output.returncode, flist))
            # remove incoplete files
            for incomplete_file in self.targetdir.glob('*.aria2'):
                incomplete_file.unlink()
                self.targetdir.joinpath(incomplete_file.stem).unlink()
        else:
            os.remove(flist)

        self.files = [self.targetdir.joinpath(basename(x)) for x in self.modis_urls if self.targetdir.joinpath(basename(x)).exists()]

        print('\n[{}]: Downloading finished.'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))