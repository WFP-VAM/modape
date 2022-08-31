#!/bin/bash
if [ ! -d "/var/modape" ] 
then
	cd /var
	git clone https://github.com/WFP-VAM/modape.git
	cd modape && git checkout tags/v1.0 -b modape-v1.0
	python -m pip install .

	mkdir /var/storage/arc_modape_ndvi
	cd /var/storage/arc_modape_ndvi
	ln -s /var/modape/modape/scripts/modis_download.py .
	ln -s /var/modape/modape/scripts/modis_collect.py .
	ln -s /var/modape/modape/scripts/modis_smooth.py .
	ln -s /var/modape/modape/scripts/modis_window.py .
fi
cd /var/storage/arc_modape_ndvi
# A1
export TILES="h15v07,h16v06,h16v07,h16v08,h17v05,h17v06,h17v07,h17v08,\
h18v05,h18v06,h18v07,h18v08,h18v09,h19v05,h19v06,h19v07,h19v08,h19v09,\
h19v10,h19v11,h19v12,h20v05,h20v06,h20v07,h20v08,h20v09,h20v10,h20v11,\
h20v12,h21v06,h21v07,h21v08,h21v09,h21v10,h21v11,h22v07,h22v08,h22v09,\
h22v10,h22v11,h23v07,h23v08"
export CMR_USERNAME=africanriskcapacity
export CMR_PASSWORD=Nasa4ARC!
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2002-07-04 -e 2003-06-26 M?D13A2
python modis_collect.py --interleave --cleanup .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2003-07-04 -e 2004-06-25 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2003177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2004-07-03 -e 2005-06-26 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2004177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2005-07-04 -e 2006-06-26 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2005177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2006-07-04 -e 2007-06-26 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2006177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2007-07-04 -e 2008-06-25 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2007177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2008-07-03 -e 2009-06-26 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2008177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2009-07-04 -e 2010-06-26 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2009177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2010-07-04 -e 2011-06-26 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2010177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2011-07-04 -e 2012-06-25 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2011177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2012-07-03 -e 2013-06-26 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2012177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2013-07-04 -e 2014-06-26 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2013177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2014-07-04 -e 2015-06-26 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2014177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2015-07-04 -e 2016-06-25 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2015177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2016-07-03 -e 2017-06-26 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2016177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2017-07-04 -e 2018-06-26 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2017177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2018-07-04 -e 2019-06-26 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2018177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2019-07-04 -e 2020-06-25 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2019177 .
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin \
  --tile-filter $TILES -b 2020-07-03 -e 2020-12-31 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2020177 .

# A2:
python modis_smooth.py --soptimize --tempint 10 \
  --last-collected 2020361 -d ./VIM/SMOOTH ./VIM
python modis_window.py -b 2002-07-15 -e 2020-12-25 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --last-smoothed 2020361 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH

# B1:
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin --tile-filter $TILES \
  -b 2021-01-01 -e 2021-01-01 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2020361 .
python modis_smooth.py --nsmooth 64 --nupdate 6 --tempint 10 \
  --last-collected 2021001 -d ./VIM/SMOOTH ./VIM
python modis_window.py -b 2021-01-05 -e 2021-01-05 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021001 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH
python modis_window.py -b 2020-12-25 -e 2020-12-25 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021001 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH

# B2:
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin --tile-filter $TILES \
  -b 2021-01-09 -e 2021-01-09 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2021001 .
python modis_smooth.py --nsmooth 64 --nupdate 6 --tempint 10 \
  --last-collected 2021009 -d ./VIM/SMOOTH ./VIM
python modis_window.py -b 2021-01-15 -e 2021-01-15 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021009 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH
python modis_window.py -b 2021-01-05 -e 2021-01-05 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021009 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH

# B3:
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin --tile-filter $TILES \
  -b 2021-01-17 -e 2021-01-17 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2021009 .
python modis_smooth.py --nsmooth 64 --nupdate 6 --tempint 10 \
  --last-collected 2021017 -d ./VIM/SMOOTH ./VIM
python modis_window.py -b 2021-01-15 -e 2021-01-15 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021017 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH
python modis_window.py -b 2021-01-05 -e 2021-01-05 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021017 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH

# B4:
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin --tile-filter $TILES \
  -b 2021-01-25 -e 2021-01-25 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2021017 .
python modis_smooth.py --nsmooth 64 --nupdate 6 --tempint 10 \
  --last-collected 2021025 -d ./VIM/SMOOTH ./VIM
python modis_window.py -b 2021-01-25 -e 2021-01-25 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021025 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH
python modis_window.py -b 2021-01-15 -e 2021-01-15 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021025 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH

# B5:
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin --tile-filter $TILES \
  -b 2021-02-02 -e 2021-02-02 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2021025 .
python modis_smooth.py --nsmooth 64 --nupdate 6 --tempint 10 \
  --last-collected 2021033 -d ./VIM/SMOOTH ./VIM
python modis_window.py -b 2021-02-05 -e 2021-02-05 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021033 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH
python modis_window.py -b 2021-01-25 -e 2021-01-25 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021033 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH

# B6:
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin --tile-filter $TILES \
  -b 2021-02-10 -e 2021-02-10 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2021033 .
python modis_smooth.py --nsmooth 64 --nupdate 6 --tempint 10 \
  --last-collected 2021041 -d ./VIM/SMOOTH ./VIM
python modis_window.py -b 2021-02-15 -e 2021-02-15 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021041 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH
python modis_window.py -b 2021-02-05 -e 2021-02-05 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021041 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH
python modis_window.py -b 2021-01-25 -e 2021-01-25 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021041 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH

# B7:
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin --tile-filter $TILES \
  -b 2021-02-18 -e 2021-02-18 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2021041 .
python modis_smooth.py --nsmooth 64 --nupdate 6 --tempint 10 \
  --last-collected 2021049 -d ./VIM/SMOOTH ./VIM
python modis_window.py -b 2021-02-25 -e 2021-02-25 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021049 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH
python modis_window.py -b 2021-02-15 -e 2021-02-15 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021049 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH
python modis_window.py -b 2021-02-05 -e 2021-02-05 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021049 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH

# B8:
python modis_download.py --download --multithread \
  --username=$CMR_USERNAME --password=$CMR_PASSWORD \
  --robust --target-empty --match-begin --tile-filter $TILES \
  -b 2021-02-26 -e 2021-02-26 M?D13A2
python modis_collect.py --interleave --cleanup --last-collected 2021049 .
python modis_smooth.py --nsmooth 64 --nupdate 6 --tempint 10 \
  --last-collected 2021057 -d ./VIM/SMOOTH ./VIM
python modis_window.py -b 2021-03-05 -e 2021-03-05 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021057 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH
python modis_window.py -b 2021-02-25 -e 2021-02-25 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021057 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH
python modis_window.py -b 2021-02-15 -e 2021-02-15 --clip-valid --round-int 2 \
  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \
  --overwrite --last-smoothed 2021057 -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH
