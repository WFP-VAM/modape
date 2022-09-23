#!/usr/bin/python
from http.cookiejar import CookieJar
# from urllib import urlencode
import urllib
import urllib.request
 
 
# The user credentials that will be used to authenticate access to the data
 
username = "africanriskcapacity"
password = "Nasa4ARC!"
  
 
# The url of the file we wish to retrieve
 
# url = "http://e4ftl01.cr.usgs.gov/MOLA/MYD17A3H.006/2009.01.01/MYD17A3H.A2009001.h12v05.006.2015198130546.hdf.xml"
url = 'https://e4ftl01.cr.usgs.gov/DP131/MOLT/MOD13A2.061/2021.07.28/MOD13A2.A2021209.h17v07.061.2021226040020.hdf.xml'
 
 
# Create a password manager to deal with the 401 reponse that is returned from
# Earthdata Login
 
password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
password_manager.add_password(None, "https://urs.earthdata.nasa.gov", username, password)
 
 
# Create a cookie jar for storing cookies. This is used to store and return
# the session cookie given to use by the data server (otherwise it will just
# keep sending us back to Earthdata Login to authenticate).  Ideally, we
# should use a file based cookie jar to preserve cookies between runs. This
# will make it much more efficient.
 
cookie_jar = CookieJar()
  
 
# Install all the handlers.
 
opener = urllib.request.build_opener(
    urllib.request.HTTPBasicAuthHandler(password_manager),
    #urllib2.HTTPHandler(debuglevel=1),    # Uncomment these two lines to see
    #urllib2.HTTPSHandler(debuglevel=1),   # details of the requests/responses
    urllib.request.HTTPCookieProcessor(cookie_jar))
# urllib.request.install_opener(opener)
 
 
# Create and submit the request. There are a wide range of exceptions that
# can be thrown here, including HTTPError and URLError. These should be
# caught and handled.
 
request = urllib.request.Request(url)
# response = urllib.request.urlopen(request)
response = opener.open(request)

print(response.read())
