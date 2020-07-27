"""test_cli.py: Test command line scripts."""
#pylint: disable=E0401
from pathlib import Path
import shutil
import unittest
from unittest.mock import patch, Mock, MagicMock
from click.testing import CliRunner

from modape.scripts.modis_download import cli as modis_download_cli

class TestConsoleScripts(unittest.TestCase):
    """Test class for console scripts."""

    @classmethod
    def setUpClass(cls):
        '''Set up testing class'''

        cls.testpath = Path(__name__).parent

    @classmethod
    def tearDownClass(cls):
        try:
            shutil.rmtree('__pycache__')
        except: #pylint: disable=W0702
            pass

    def setUp(self):
        """Set up a test"""
        self.runner = CliRunner()

    def test_modis_download(self):
        """Test modis_download.py."""

        test_query = MagicMock()
        test_query.nresults = 4
        test_query.download = Mock(return_value='Mocking download!')

        with patch('modape.scripts.modis_download.ModisQuery', return_value=test_query) as mocked_query:

            result = self.runner.invoke(modis_download_cli, ["MOD13A2"])

            mocked_query.assert_called()
            args = mocked_query.call_args_list[0][1]
            self.assertEqual(args['products'], ["MOD13A2"])
            assert result.exit_code == 0

            mocked_query.reset_mock()
            #
            result = self.runner.invoke(modis_download_cli, ["m?d13a2"])
            mocked_query.assert_called()
            args = mocked_query.call_args_list[0][1]
            self.assertEqual(args['products'], ["MOD13A2", 'MYD13A2'])
            assert result.exit_code == 0

            mocked_query.reset_mock()

            result = self.runner.invoke(modis_download_cli, ["MOD13A2", "--download"])
            mocked_query.assert_not_called()
            assert result.exit_code == 1

            result = self.runner.invoke(modis_download_cli, ["MOD13A2", "--username", "test", "--password", "test", "--download"])
            mocked_query.assert_called()

            mocked_query.reset_mock()

            fake_hdf = self.testpath.joinpath("MOD13A2.A2002193.h18v06.006.2019256103823.hdf")
            fake_hdf.touch()

            result = self.runner.invoke(modis_download_cli, ["MOD13A2", "--targetdir", str(self.testpath)])
            mocked_query.assert_called()
            assert result.exit_code == 0

            mocked_query.reset_mock()

            result = self.runner.invoke(modis_download_cli, ["MOD13A2", "--targetdir", str(self.testpath), "--target-empty"])
            mocked_query.assert_not_called()

            assert result.exit_code == 1

            fake_hdf.unlink()

            result = self.runner.invoke(modis_download_cli, ["MOD13A2", "--targetdir", str(self.testpath), "--target-empty"])
            mocked_query.assert_called()

            assert result.exit_code == 0

            mocked_query.reset_mock()

            result = self.runner.invoke(modis_download_cli, ["MOD13A2", "--roi", "10,10"])
            mocked_query.assert_called()
            args = mocked_query.call_args_list[0][1]
            self.assertEqual(args['aoi'], (10.0, 10.0))
            assert result.exit_code == 0

            mocked_query.reset_mock()

            result = self.runner.invoke(modis_download_cli, ["MOD13A2", "--roi", "10,10,20,20"])
            mocked_query.assert_called()
            args = mocked_query.call_args_list[0][1]
            self.assertEqual(args['aoi'], (10.0, 10.0, 20.0, 20.0))
            assert result.exit_code == 0

            mocked_query.reset_mock()

            result = self.runner.invoke(modis_download_cli, ["MOD13A2", "--roi", "10,20,20"])
            mocked_query.assert_not_called()
            assert result.exit_code == 1

            mocked_query.reset_mock()

            result = self.runner.invoke(modis_download_cli, ["MOD13A2", "--roi", "10,20,20,30,50"])
            mocked_query.assert_not_called()
            assert result.exit_code == 1
