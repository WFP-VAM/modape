from __future__ import print_function,division
import requests
from bs4 import BeautifulSoup
import re
import sys


class MODISquery:

    def __init__(self,url):

        r = re.compile(".+(h\d+v\d+).+")

        self.queryURL = url

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
