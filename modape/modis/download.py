"""
MODIS query & download class.

This file contains the class used for querying and downloading
raw MODIS products from NASA's servers.

Author: Valentin Pesendorfer, July 2020
"""

# pylint: disable=E0401, E0611

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
import re
import shutil
from typing import List, Tuple, Union

from cmr import GranuleQuery
import pandas as pd
from requests.exceptions import HTTPError


from utils import SessionWithHeaderRedirection

class DownloadError(Exception):
    def __init__(self, fails):

        message = '''

ERROR downloading MODIS data!

Failed downloads:

'''
        for failed, error in fails:
            message += f"{failed}: {error}\n"

        super(DownloadError, self).__init__(message)


class ModisQuery(object):
    """Class for querying and downloading MODIS data."""

    def __init__(self,
                 products: List[str],
                 #targetdir: Path,
                 aoi: List[Union[float, int]] = None,
                 begindate: datetime = None,
                 enddate: datetime = None,
                 tile_filter: List[str] = None,
                 version: str = "006") -> None:

        assert products, 'No product IDs supplied!'

        self.begin = begindate
        self.end = enddate
        self.tile_filter = tile_filter
        self.api = GranuleQuery()

        # construct query

        self.api.parameters(
            short_name=products,
            version=version,
            temporal=(begindate, enddate)
        )

        if aoi is not None:
            if len(aoi) == 2:
                self.api.point(*aoi)
            elif len(aoi) == 4:
                self.api.bounding_box(*aoi)
            else:
                raise ValueError("Expected point or bounding box as AOI")

    def query(self, strict_dates: bool = True) -> None:

        self.results = []

        if self.begin is None and self.end is None:
            strict_dates = False

        results_all = self.api.get_all()

        # TODO: check for duplicates?

        for result in self._parse_response(results_all):

            if self.tile_filter and result['tile']:
                if result['tile'] not in self.tile_filter:
                    continue

            if strict_dates:
                if self.begin is not None:
                    if result['time_start'] < self.begin.date():
                        continue

                if self.end is not None:
                    if result['time_stop'] > self.end.date():
                        continue

            self.results.append(result)

        self.nresults = len(self.results)


    @staticmethod
    def _parse_response(query: List[dict]) -> dict:

        tile_regxp = re.compile(r'.+(h\d+v\d+).+')

        for entry in query:

            entry_parsed = dict(
                file_id=entry['producer_granule_id'],
                time_start=pd.Timestamp(entry['time_start']).date(),
                time_end=pd.Timestamp(entry['time_end']).date(),
                updated=entry['updated'],
                link=entry['links'][0]['href'],
            )

            try:
                tile = tile_regxp.search(entry_parsed['file_id']).group(1)
            except AttributeError:
                tile = None

            entry_parsed.update({"tile": tile})

            yield entry_parsed

    @staticmethod
    def _fetch_hdf(session: SessionWithHeaderRedirection,
                   url: str,
                   destination: Path
                  ) -> Tuple[str, bool]:

        try:

            filename = destination.joinpath(url.split('/')[-1])

            with session.get(url, stream=True) as response:
                response.raise_for_status()

                with open(filename, "wb") as openfile:
                    shutil.copyfileobj(response.raw, openfile)

            assert filename.exists()

        except (HTTPError, AssertionError) as e:
            return (url, e)

        return (filename, None)


    def download(self,
                 targetdir: Path,
                 username: str,
                 password: str,
                 multithread: bool = False
                ) -> None:

        assert targetdir.exists()
        assert targetdir.is_dir()

        with SessionWithHeaderRedirection(username, password) as session:

            if multithread:

                with ThreadPoolExecutor(4) as pool:

                    futures = [pool.submit(self._fetch_hdf, session, x['link'], targetdir)
                               for x in self.results]

                downloaded_files = [x.result() for x in futures]

            else:

                downloaded_files = []

                for result in self.results:

                    downloaded_files.append(
                        self._fetch_hdf(session, result['link'], targetdir)
                    )

        errors = []

        for file, err in downloaded_files:

            if err is not None:
                errors.append((file, err))
            else:
                assert file.exists(), "Downloaded file is missing! No download error was reported"

        if errors:
            raise DownloadError(errors)
