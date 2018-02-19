from __future__ import print_function
import requests
from bs4 import BeautifulSoup
import re
import sys, os
import glob
import time
from subprocess import Popen, check_output
try:
    import gdal
except ImportError:
    from osgeo import gdal


class MODISquery:

    def __init__(self,url,username=None,password=None,rawdir=os.getcwd(),targetdir=os.getcwd()):

        self.queryURL = url
        self.username = username
        self.password = password
        self.rawdir = rawdir
        self.targetdir = targetdir
        self.files = []

        r = re.compile(".+(h\d+v\d+).+")

        print('Checking for MODIS products ...')
        try:
            response = requests.get(url)
            self.statuscode = response.status_code
            response.raise_for_status()

        except requests.exceptions.RequestException as e:
            print(e)
            sys.exit(1)

        soup = BeautifulSoup(response.content,"html5lib")

        self.modisURLs = [x.getText() for x in soup.find_all('url')]
        self.results = len(self.modisURLs)
        self.tiles = list(set([r.search(x).group(1) for x in self.modisURLs]))

        print('... done.\n')

        if self.results > 0:
            print('%s results found.' % self.results)
        else:
            print('0 results found. Please check query!')


    def setCredentials(self,username,password):
        self.username=username
        self.password=password

    def download(self):

        try:
            temp = check_output(['wget', '--version'])
        except:
            print("Download needs wget to be available in PATH! Please make sure it's installed and available in PATH!")
            sys.exit(1)

        if self.username is None or self.password is None:
            print('No credentials found. Please run .setCredentials(username,password)!')
            sys.exit(1)


        args = ['wget','--user',self.username,'--password',self.password,'-P',self.rawdir]

        print('[%s]: Downloading products to %s ...\n' % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),self.rawdir))

        for ix,u in enumerate(self.modisURLs):
            print('%s of %s' %(ix+1,self.results))
            p = Popen(args + [u])
            p.wait()
            if p.returncode is not 0:
                print("Couldn't download %s - continuing." % u)
                continue
            self.files = self.files + [self.rawdir + os.path.basename(u)]

        print('\n[%s]: Downloading finished.' % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))


    def process(self):

        print('[%s]: Starting processing ...\n' % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))











def MODISprocess(modisHDFs,dstorage=None,rdstorage=None):
    pass
