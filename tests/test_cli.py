"""test_cli.py: Test command line scripts."""
#pylint: disable=E0401
from pathlib import Path
import shutil
import unittest
from unittest.mock import patch, Mock, MagicMock, call
from click.testing import CliRunner

from modape.scripts.modis_download import cli as modis_download_cli
from modape.scripts.modis_collect import cli as modis_collect_cli

class TestConsoleScripts(unittest.TestCase):
    """Test class for console scripts."""

    @classmethod
    def setUpClass(cls):
        '''Set up testing class'''

        cls.testpath = Path('/tmp')

        cls.lst_files = ['MYD11A2.A2002193.h18v06.006.2015146152945.hdf',
                         'MOD11A2.A2002209.h18v06.006.2015145152020.hdf',
                         'MYD11A2.A2002201.h18v06.006.2015146153241.hdf',
                         'MYD11A2.A2002185.h18v06.006.2015146152642.hdf',
                         'MOD11A2.A2002177.h18v06.006.2015144183717.hdf',
                         'MYD11A2.A2002209.h18v06.006.2015152152813.hdf',
                         'MOD11A2.A2002185.h18v06.006.2015145002847.hdf',
                         'MOD11A2.A2002193.h18v06.006.2015145055806.hdf',
                         'MOD11A2.A2002201.h18v06.006.2015145105749.hdf']

        cls.vim_files = ['MYD13A2.A2002201.h18v06.006.2015149071105.hdf',
                         'MYD13A2.A2002185.h18v06.006.2015149071113.hdf',
                         'MOD13A2.A2002177.h18v06.006.2015149001129.hdf',
                         'MOD13A2.A2002209.h18v06.006.2015149180726.hdf',
                         'MOD13A2.A2002193.h18v06.006.2015149022847.hdf']

    @classmethod
    def tearDownClass(cls):
        try:
            shutil.rmtree('__pycache__')
        except: #pylint: disable=W0702
            pass

        try:
            shutil.rmtree('/tmp/data')
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

    def test_modis_collect(self):
        """Test modis_collect.py"""
        result = self.runner.invoke(modis_collect_cli, ["/not_an_exist_dir"])
        self.assertEqual(result.exit_code, 1)

        testfile = Path('/tmp/file.txt')
        testfile.touch()
        result = self.runner.invoke(modis_collect_cli, [str(testfile)])
        self.assertEqual(result.exit_code, 1)
        testfile.unlink()

        result = self.runner.invoke(modis_collect_cli, ["/tmp"])
        self.assertEqual(result.exit_code, 1)

        data_dir = self.testpath.joinpath('data')
        data_dir.mkdir()

        for file in self.lst_files:
            file_path = data_dir.joinpath(file)
            file_path.touch()

        calls = []
        for product in ['MOD11', 'MYD11']:
            product_files = [str(x) for x in data_dir.glob("*hdf") if product in x.name]
            product_files.sort()
            calls.append(
                call(product_files, data_dir, None, False, "gzip"),
            )

        with patch("modape.scripts.modis_collect._worker") as mocked_worker:
            result = self.runner.invoke(modis_collect_cli, ["/tmp/data"])
            mocked_worker.assert_called()
            self.assertEqual(mocked_worker.call_count, 2)
            mocked_worker.assert_has_calls(calls, any_order=False)

        _ = [x.unlink() for x in data_dir.glob("*hdf")]

        for file in self.vim_files:
            file_path = data_dir.joinpath(file)
            file_path.touch()

        with patch("modape.scripts.modis_collect._worker") as mocked_worker:
            result = self.runner.invoke(modis_collect_cli, ["/tmp/data", "--interleave"])
            mocked_worker.assert_called_once()
            product_files = [str(x) for x in data_dir.glob("*hdf")]
            product_files.sort()
            mocked_worker.assert_called_with(product_files, data_dir, None, True, "gzip")

        _ = [x.unlink() for x in data_dir.glob("*hdf")]
