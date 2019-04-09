"""test_cli.py: Test command line scripts."""
from __future__ import absolute_import, division, print_function
import unittest

from modape.scripts.csv_smooth import main as csv_smooth_main
from modape.scripts.modis_download import main as modis_download_main
from modape.scripts.modis_collect import main as modis_collect_main
from modape.scripts.modis_smooth import main as modis_smooth_main
from modape.scripts.modis_window import main as modis_window_main
from modape.scripts.modis_info import main as modis_info_main
from modape.scripts.modis_product_table import main as modis_product_table_main
from modape.scripts.rts_smooth import main as rts_smooth_main

class TestConsoleScripts(unittest.TestCase):
    """Test class for console scripts."""

    def test_csv_smooth(self):
        """Test csv_smooth.py."""
        try:
            csv_smooth_main()
        except SystemExit as system_exit:
            if system_exit.code == 0:
                pass
            else:
                raise

    def test_modis_download(self):
        """Test modis_download.py."""
        try:
            modis_download_main()
        except SystemExit as system_exit:
            if system_exit.code == 0:
                pass
            else:
                raise

    def test_modis_collect(self):
        """Test modis_collect.py."""
        try:
            modis_collect_main()
        except SystemExit as system_exit:
            if system_exit.code == 0:
                pass
            else:
                raise

    def test_modis_smooth(self):
        """Test modis_smooth.py."""
        try:
            modis_smooth_main()
        except SystemExit as system_exit:
            if system_exit.code == 0:
                pass
            else:
                raise

    def test_modis_window(self):
        """Test modis_window.py."""
        try:
            modis_window_main()
        except SystemExit as system_exit:
            if system_exit.code == 0:
                pass
            else:
                raise


    def test_modis_info(self):
        """Test modis_info.py."""
        try:
            modis_info_main()
        except SystemExit as system_exit:
            if system_exit.code == 0:
                pass
            else:
                raise

    def test_modis_product_table(self):
        """Test modis_product_table.py."""
        try:
            modis_product_table_main()
        except SystemExit as system_exit:
            if system_exit.code == 0:
                pass
            else:
                raise

    def test_rts_smooth(self):
        """Test rts_smooth.py."""
        try:
            rts_smooth_main()
        except SystemExit as system_exit:
            if system_exit.code == 0:
                pass
            else:
                raise
