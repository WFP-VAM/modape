## Fake responses for testing

tiled = b'''<?xml version="1.0"?>
<inventory>
<url>https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.02.18/MOD13A2.A2000049.h18v06.006.2015136104646.hdf</url>
<url>https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.03.05/MOD13A2.A2000065.h18v06.006.2015136022922.hdf</url>
<url>https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.03.21/MOD13A2.A2000081.h18v06.006.2015136035955.hdf</url>
<url>https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.04.06/MOD13A2.A2000097.h18v06.006.2015136035959.hdf</url>
<url>https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.04.22/MOD13A2.A2000113.h18v06.006.2015137034359.hdf</url>
</inventory>
'''

glob =  b'''<img src="/icons/blank.gif" alt="Icon "> <a href="?C=N;O=D">Name</a>                    <a href="?C=M;O=A">Last modified</a>      <a href="?C=S;O=A">Size</a>  <a href="?C=D;O=A">Description</a><hr><img src="/icons/back.gif" alt="[PARENTDIR]"> <a href="/MOLA/">Parent Directory</a>                             -
<img src="/icons/folder.gif" alt="[DIR]"> <a href="2002.07.04/">2002.07.04/</a>             2016-03-22 15:43    -
<img src="/icons/folder.gif" alt="[DIR]"> <a href="2002.07.12/">2002.07.12/</a>             2016-03-22 15:42    -
<img src="/icons/folder.gif" alt="[DIR]"> <a href="2002.07.20/">2002.07.20/</a>             2016-03-22 15:44    -
<img src="/icons/folder.gif" alt="[DIR]"> <a href="2002.07.28/">2002.07.28/</a>             2016-03-22 15:42    -
<img src="/icons/folder.gif" alt="[DIR]"> <a href="2002.08.05/">2002.08.05/</a>             2016-03-22 15:43    -
</hr>
'''

mola = {
'http://global-test.query/2002.07.04/':
b'''<img src="/icons/blank.gif" alt="Icon "> <a href="?C=N;O=D">Name</a>                                            <a href="?C=M;O=A">Last modified</a>      <a href="?C=S;O=A">Size</a>  <a href="?C=D;O=A">Description</a><hr><img src="/icons/back.gif" alt="[PARENTDIR]"> <a href="/MOLA/MYD11C2.006/">Parent Directory</a>                                                     -
<img src="/icons/image2.gif" alt="[IMG]"> <a href="BROWSE.MYD11C2.A2002185.006.2015168205556.1.jpg">BROWSE.MYD11C2.A2002185.006.2015168205556.1.jpg</a> 2015-06-17 16:40  3.0M
<img src="/icons/image2.gif" alt="[IMG]"> <a href="BROWSE.MYD11C2.A2002185.006.2015168205556.2.jpg">BROWSE.MYD11C2.A2002185.006.2015168205556.2.jpg</a> 2015-06-17 16:40  2.8M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="MYD11C2.A2002185.006.2015168205556.hdf">MYD11C2.A2002185.006.2015168205556.hdf</a>          2015-06-17 16:39   57M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="MYD11C2.A2002185.006.2015168205556.hdf.xml">MYD11C2.A2002185.006.2015168205556.hdf.xml</a>      2015-06-17 16:39  6.0K
<hr>''',
'http://global-test.query/2002.07.12/':
b'''<img src="/icons/blank.gif" alt="Icon "> <a href="?C=N;O=D">Name</a>                                            <a href="?C=M;O=A">Last modified</a>      <a href="?C=S;O=A">Size</a>  <a href="?C=D;O=A">Description</a><hr><img src="/icons/back.gif" alt="[PARENTDIR]"> <a href="/MOLA/MYD11C2.006/">Parent Directory</a>                                                     -
<img src="/icons/image2.gif" alt="[IMG]"> <a href="BROWSE.MYD11C2.A2002193.006.2015149021321.1.jpg">BROWSE.MYD11C2.A2002193.006.2015149021321.1.jpg</a> 2015-05-29 00:07  3.1M
<img src="/icons/image2.gif" alt="[IMG]"> <a href="BROWSE.MYD11C2.A2002193.006.2015149021321.2.jpg">BROWSE.MYD11C2.A2002193.006.2015149021321.2.jpg</a> 2015-05-29 00:07  2.9M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="MYD11C2.A2002193.006.2015149021321.hdf">MYD11C2.A2002193.006.2015149021321.hdf</a>          2015-05-29 00:06   57M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="MYD11C2.A2002193.006.2015149021321.hdf.xml">MYD11C2.A2002193.006.2015149021321.hdf.xml</a>      2015-05-29 00:06  6.0K
<hr>''',
'http://global-test.query/2002.07.20/':
b'''<img src="/icons/blank.gif" alt="Icon "> <a href="?C=N;O=D">Name</a>                                            <a href="?C=M;O=A">Last modified</a>      <a href="?C=S;O=A">Size</a>  <a href="?C=D;O=A">Description</a><hr><img src="/icons/back.gif" alt="[PARENTDIR]"> <a href="/MOLA/MYD11C2.006/">Parent Directory</a>                                                     -
<img src="/icons/image2.gif" alt="[IMG]"> <a href="BROWSE.MYD11C2.A2002201.006.2015149021446.1.jpg">BROWSE.MYD11C2.A2002201.006.2015149021446.1.jpg</a> 2015-05-29 01:53  3.1M
<img src="/icons/image2.gif" alt="[IMG]"> <a href="BROWSE.MYD11C2.A2002201.006.2015149021446.2.jpg">BROWSE.MYD11C2.A2002201.006.2015149021446.2.jpg</a> 2015-05-29 01:53  2.9M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="MYD11C2.A2002201.006.2015149021446.hdf">MYD11C2.A2002201.006.2015149021446.hdf</a>          2015-05-29 01:53   57M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="MYD11C2.A2002201.006.2015149021446.hdf.xml">MYD11C2.A2002201.006.2015149021446.hdf.xml</a>      2015-05-29 01:53  6.0K
<hr>''',
'http://global-test.query/2002.07.28/':
b'''<img src="/icons/blank.gif" alt="Icon "> <a href="?C=N;O=D">Name</a>                                            <a href="?C=M;O=A">Last modified</a>      <a href="?C=S;O=A">Size</a>  <a href="?C=D;O=A">Description</a><hr><img src="/icons/back.gif" alt="[PARENTDIR]"> <a href="/MOLA/MYD11C2.006/">Parent Directory</a>                                                     -
<img src="/icons/image2.gif" alt="[IMG]"> <a href="BROWSE.MYD11C2.A2002209.006.2015149021240.1.jpg">BROWSE.MYD11C2.A2002209.006.2015149021240.1.jpg</a> 2015-05-29 00:40  3.8M
<img src="/icons/image2.gif" alt="[IMG]"> <a href="BROWSE.MYD11C2.A2002209.006.2015149021240.2.jpg">BROWSE.MYD11C2.A2002209.006.2015149021240.2.jpg</a> 2015-05-29 00:40  3.5M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="MYD11C2.A2002209.006.2015149021240.hdf">MYD11C2.A2002209.006.2015149021240.hdf</a>          2015-05-29 00:40   47M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="MYD11C2.A2002209.006.2015149021240.hdf.xml">MYD11C2.A2002209.006.2015149021240.hdf.xml</a>      2015-05-29 00:40  5.6K
<hr>''',
'http://global-test.query/2002.08.05/':
b'''<img src="/icons/blank.gif" alt="Icon "> <a href="?C=N;O=D">Name</a>                                            <a href="?C=M;O=A">Last modified</a>      <a href="?C=S;O=A">Size</a>  <a href="?C=D;O=A">Description</a><hr><img src="/icons/back.gif" alt="[PARENTDIR]"> <a href="/MOLA/MYD11C2.006/">Parent Directory</a>                                                     -
<img src="/icons/image2.gif" alt="[IMG]"> <a href="BROWSE.MYD11C2.A2002217.006.2015149021044.1.jpg">BROWSE.MYD11C2.A2002217.006.2015149021044.1.jpg</a> 2015-05-29 00:40  3.4M
<img src="/icons/image2.gif" alt="[IMG]"> <a href="BROWSE.MYD11C2.A2002217.006.2015149021044.2.jpg">BROWSE.MYD11C2.A2002217.006.2015149021044.2.jpg</a> 2015-05-29 00:40  3.0M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="MYD11C2.A2002217.006.2015149021044.hdf">MYD11C2.A2002217.006.2015149021044.hdf</a>          2015-05-29 00:40   55M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="MYD11C2.A2002217.006.2015149021044.hdf.xml">MYD11C2.A2002217.006.2015149021044.hdf.xml</a>      2015-05-29 00:40  5.8K
<hr>'''
}
