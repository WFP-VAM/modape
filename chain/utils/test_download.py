

import requests
from base64 import b64encode
 
# ***********************
# overriding requests.Session.rebuild_auth to maintain headers when redirected
# ***********************
class SessionWithHeaderRedirection(requests.Session):
    AUTH_HOST = 'urs.earthdata.nasa.gov'

    def __init__(self, username, password):
        super().__init__()
        self.__username = username
        self.__password = password

    # Overrides from the library to keep headers when redirected to or from the NASA auth host.
    def rebuild_auth(self, prepared_request, response):
        
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
 
# create session with the user credentials that will be used to authenticate access to the data

username = "africanriskcapacity"
password= "Nasa4ARC!"
 
session = SessionWithHeaderRedirection(username, password)
 
# the url of the file we wish to retrieve
url = 'https://e4ftl01.cr.usgs.gov/DP131/MOLT/MOD13A2.061/2021.07.28/MOD13A2.A2021209.h17v07.061.2021226040020.hdf'
# url = "https://e4ftl01.cr.usgs.gov/MOLA/MYD17A3H.006/2009.01.01/MYD17A3H.A2009001.h12v05.006.2015198130546.hdf.xml"
url = 'https://e4ftl01.cr.usgs.gov/DP131/MOLT/MOD13A2.061/2021.07.28/MOD13A2.A2021209.h17v07.061.2021226040020.hdf.xml'
 
# extract the filename from the url to be used when saving the file
filename = url[url.rfind('/')+1:]  

r = session.get(url, stream=True, allow_redirects=True) #, auth=(username, password))

print(r.content)

