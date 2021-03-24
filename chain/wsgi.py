"""
    Run as Gunicorn service from within the container: gunicorn --workers=1 --threads=1 --bind 0.0.0.0:5001 wsgi:arc_modis_ndvi_server
"""

from arc_modis_ndvi import app_setup
arc_modis_ndvi_server = app_setup('production.json')
