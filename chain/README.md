The ARC MODIS Filtered NDVI Production Chain is developed in Python 3, building on the WFP VAM MODAPE toolkit.

For development, basically the same environment is needed as for production. Open a command line window / terminal in the directory 
you intend to work in; then, as outlined in the deployment SOP, issue the following commands:

  mkvirtualenv arc-modis-ndvi
  pip install numpy
  pip install cython gunicorn flask gdal==3.2.0
  
  git clone https://github.com/WFP-VAM/modape.git
  cd modape
  git checkout -b v1.0rc origin/v1.0rc
  pip install .

Python FILES
------------
- arc_modis_ndvi.py -- actual processing chain (if you don't know where to start working on the MODIS processing, start here)
- modape_calendar.py -- utility script to produce the Release Calendar and MODAPE CLI parameters

Scripts
-------
- arc_modape_run.sh -- bootstrapper script to run another script (in $PWD) in a docker container
- policy_dataset_test.sh -- script to reproduce the policy dataset (starting from line 19, this script is in the policy annex)
- arc_modape_bash.sh -- start a terminal session in the docker container
- 