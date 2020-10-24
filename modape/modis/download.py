"""
MODIS query & download class.

This file contains the class used for querying and downloading
raw MODIS products from NASA's servers.

Author: Valentin Pesendorfer, July 2020
"""

# pylint: disable=E0401, E0611

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import logging
from os.path import exists
from pathlib import Path
import re
import shutil
from typing import List, Tuple, Union

from cmr import GranuleQuery
import pandas as pd
from requests.adapters import HTTPAdapter
from requests.exceptions import HTTPError
from requests.packages.urllib3.util.retry import Retry

from modape.exceptions import DownloadError
from modape.utils import SessionWithHeaderRedirection

log = logging.getLogger(__name__)

class ModisQuery(object):
    """Class for querying and downloading MODIS data."""

    def __init__(self,
                 products: List[str],
                 aoi: List[Union[float, int]] = None,
                 begindate: datetime = None,
                 enddate: datetime = None,
                 tile_filter: List[str] = None,
                 version: str = "006") -> None:
        """Initialize ModisQuery instance.

        This class is used for querying and downloading MODIS data
        from MODIS CMR.

        Args:
            products (List[str]): List of product codes to be queried / downloaded.
            aoi (List[Union[float, int]]): Area of interes (point as lat/lon or bounding box as xmin, ymin, xmax, ymax).
            begindate (datetime): Start date for query.
            enddate (datetime): End date for query.
            tile_filter (List[str]): List of tiles to be queried / downloaded (refines initial results).
            version (str): MODIS collection version.
        """

        assert products, "No product IDs supplied!"

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

    def search(self, strict_dates: bool = True) -> None:
        """Query MODIS data.

        Performs query based on inut to class instance.

        Args:
            strict_dates (bool): Flag for strict date enforcement (no data with timestamp outside
                                 of begindate and enddate are allowed).
        """

        # init results list
        self.results = []

        # if no dates supplied, we can't be strict
        if self.begin is None and self.end is None:
            strict_dates = False

        log.debug("Starting query")

        # get all results
        results_all = self.api.get_all()

        log.debug("Query complete, filtering results")

        for result in self._parse_response(results_all):

            # skip tiles outside of filter
            if self.tile_filter and result["tile"]:
                if result["tile"] not in self.tile_filter:
                    continue

            # enforce begin date if required

            if strict_dates:
                if self.begin is not None:
                    if result["time_start"] < self.begin.date():
                        continue

                if self.end is not None:
                    if result["time_start"] > self.end.date():
                        continue

            self.results.append(result)

        # final results
        self.nresults = len(self.results)

        log.debug("Search complete. Total results: %s, filtered: %s", len(results_all), self.nresults)


    @staticmethod
    def _parse_response(query: List[dict]) -> dict:
        """Generator for parsing API response.

        Args:
            query (List[dict]): Query returned by CMR API.

        Returns:
            dict: Parsed query as dict.

        """

        tile_regxp = re.compile(r".+(h\d+v\d+).+")

        for entry in query:

            entry_parsed = dict(
                file_id=entry["producer_granule_id"],
                time_start=pd.Timestamp(entry["time_start"]).date(),
                time_end=pd.Timestamp(entry["time_end"]).date(),
                updated=entry["updated"],
                link=entry["links"][0]["href"],
            )

            try:
                tile = tile_regxp.search(entry_parsed["file_id"]).group(1)
            except AttributeError:
                tile = None

            entry_parsed.update({"tile": tile})

            yield entry_parsed

    @staticmethod
    def _fetch_hdf(session: SessionWithHeaderRedirection,
                   url: str,
                   destination: Path,
                   overwrite: bool,
                   ) -> Tuple[str, Union[None, Exception]]:
        """Helper function to fetch HDF files

        Args:
            session (SessionWithHeaderRedirection): requests session to fetch file.
            url (str): URL for file.
            destination (Path): target directory.
            overwrite (bool): Overwrite existing.

        Returns:
            Tuple[str, Union[None, Exception]]: Returns tuple with
                either (filename, None) for success and (URL, Exception) for error.

        """

        filename = destination.joinpath(url.split("/")[-1])

        if not exists(filename) or overwrite:

            try:

                with session.get(url, stream=True, allow_redirects=True) as response:
                    response.raise_for_status()

                    with open(filename, "wb") as openfile:
                        shutil.copyfileobj(response.raw, openfile, length=16*1024*1024)#

                assert filename.exists(), "File not on disk after download!"

            except (HTTPError, AssertionError) as e:
                return (url, e)
        else:
            log.info("%s exists in target. Please set overwrite to True.", filename)

        return (filename, None)

    def download(self,
                 targetdir: Path,
                 username: str,
                 password: str,
                 overwrite: bool = False,
                 multithread: bool = False,
                 nthreads: int = 4,
                ) -> None:
        """Download queried MODIS files.

        Args:
            targetdir (Path): target directory for file being downloaded.
            username (str): Earthdata username.
            password (str): Earthdata password.
            overwrite (bool): Replace existing.
            multithread (bool): use multiple threads for downloading.
            nthreads (int): Number of threads.

        """

        # make sure target directory is dir and exists
        assert targetdir.is_dir()
        assert targetdir.exists()

        with SessionWithHeaderRedirection(username, password) as session:

            retries = Retry(total=5, backoff_factor=1, status_forcelist=[502, 503, 504])
            session.mount(
                "https://",
                HTTPAdapter(pool_connections=nthreads, pool_maxsize=nthreads*2, max_retries=retries)
            )

            if multithread:
                log.debug("Multithreaded download using %s threads. Warming up connection pool.", nthreads)
                # warm up pool
                _ = session.get(self.results[0]["link"], stream=True, allow_redirects=True)

                with ThreadPoolExecutor(nthreads) as executor:

                    futures = [executor.submit(self._fetch_hdf, session, x["link"], targetdir, overwrite)
                               for x in self.results]

                downloaded_files = [x.result() for x in futures]

            else:
                log.debug("Serial download")
                downloaded_files = []

                for result in self.results:

                    downloaded_files.append(
                        self._fetch_hdf(session, result["link"], targetdir, overwrite)
                    )

        errors = []

        # check if downloads are OK
        for file, err in downloaded_files:

            if err is not None:
                errors.append((file, err)) # append to error list
            else:
                # if no error, make sure file is on disk
                assert file.exists(), "Downloaded file is missing! No download error was reported"

        if errors:
            raise DownloadError(errors)
