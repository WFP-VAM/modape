"""
MODIS query & download class.

This file contains the class used for querying and downloading
raw MODIS products from NASA's servers.
"""

# pylint: disable=E0401, E0611

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import logging
from os.path import exists
from pathlib import Path
import re
import shutil
import hashlib
from typing import List, Tuple, Union
from xml.etree import ElementTree

from cmr import GranuleQuery
import pandas as pd
from pycksum import cksum
from requests.adapters import HTTPAdapter
from requests.exceptions import HTTPError, ConnectionError #pylint: disable=W0622
from urllib3.util import Retry

from modape.exceptions import DownloadError
from modape.utils import SessionWithHeaderRedirection

log = logging.getLogger(__name__)

class ModisQuery(object):
    """Class for querying and downloading MODIS data.

    The user can create a query and send it to NASA's CMR servers.
    The response can be either just printed to console or passed to
    the `download` method, to fetch the resulting HDF files to local disk.
    """

    def __init__(self,
                 products: List[str],
                 aoi: List[Union[float, int]] = None,
                 begindate: datetime = None,
                 enddate: datetime = None,
                 tile_filter: List[str] = None,
                 version: str = "006") -> None:
        """Initialize instance ModisQuery class.

        This creates an instance of `ModisQuery` with the basic query parameters.
        The `aoi` needs to be a list if either 2 or 4 coordinates as `float` or `int`, depending
        on if a point or bounding box is requested.
        The tile_filter needs to be specified as list of MODIS tile IDs in format hXXvXX (e.g. `h20v08`).


        Args:
            products (List[str]): List of product codes to be queried / downloaded.
            aoi (List[Union[float, int]]): Area of interes (point as lon/lat or bounding box as xmin, ymin, xmax, ymax).
            begindate (datetime): Start date for query.
            enddate (datetime): End date for query.
            tile_filter (List[str]): List of tiles to be queried / downloaded (refines initial results).
            version (str): MODIS collection version.

        Raises:
            AssertionError: If no product code is supplied or `if len(aoi) not in [2, 4]`.
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

    def search(self, match_begin: bool = True) -> None:
        """Send quert to MODIS CMR servers.

        Constructs the query from parameters passed to `__init__`
        and sends the query to the NASA servers. The returned results
        will be stored in a class variable.
        To deal with overlapping date ranges of composite products,
        the specified start date can be matched to the MODIS native timestamp.

        Args:
            match_begin (bool): Flag to match begin date with native MODIS timestamp (no data with timestamp earlier than begindate is allowed).
        """

        # init results dict
        self.results = {}

        # if no dates supplied, we can't be match
        if self.begin is None and self.end is None:
            match_begin = False

        log.debug("Starting query")

        # get all results
        results_all = self.api.get_all()

        log.debug("Query complete, filtering results")

        for result in self._parse_response(results_all):

            # skip tiles outside of filter
            if self.tile_filter and result["tile"]:
                if result["tile"] not in self.tile_filter:
                    continue

            # enforce dates if required

            if match_begin:
                if self.begin is not None:
                    if result["time_start"] < self.begin.date():
                        continue

                if self.end is not None:
                    if result["time_start"] > self.end.date():
                        continue

            filename = result["filename"]
            del result["filename"]

            self.results.update({filename: result})

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
                filename=entry["producer_granule_id"],
                time_start=pd.Timestamp(entry["time_start"]).date(),
                time_end=pd.Timestamp(entry["time_end"]).date(),
                updated=entry["updated"],
                link=entry["links"][0]["href"],
            )

            try:
                tile = tile_regxp.search(entry_parsed["filename"]).group(1)
            except AttributeError:
                tile = None

            entry_parsed.update({"tile": tile})

            yield entry_parsed

    @staticmethod
    def _parse_hdfxml(response):
        result = {}
        tree = ElementTree.fromstring(response.content)
        for entry in tree.iter(tag='GranuleURMetaData'):
            for datafile in entry.iter(tag="DataFiles"):
                for datafilecont in datafile.iter(tag="DataFileContainer"):
                    for content in datafilecont:
                        if content.tag in ("FileSize", "ChecksumType", "Checksum"):
                            result.update({content.tag: content.text.strip()})
        return result
    
    @staticmethod
    def _parse_cmrxml(response, hdf_filename):
        result = {}
        tree = ElementTree.fromstring(response.content)
        entry = tree.find(f"DataGranule/AdditionalFile[Name = '{hdf_filename}']")
        result.update({"FileSize": entry.find("SizeInBytes").text})
        result.update({"ChecksumType": entry.find("Checksum/Algorithm").text})
        result.update({"Checksum": entry.find("Checksum/Value").text})
        return result


    def _fetch(self,
               session: SessionWithHeaderRedirection,
               url: str,
               destination: Path,
               overwrite: bool,
               check: bool,
               ) -> Tuple[str, Union[None, Exception]]:
        """Helper function to fetch HDF files

        Args:
            session (SessionWithHeaderRedirection): requests session to fetch file.
            url (str): URL for file.
            destination (Path): Target directory.
            overwrite (bool): Overwrite existing.
            check (bool): Check file size and checksum.

        Returns:
            Tuple[str, Union[None, Exception]]: Returns tuple with
                either (filename, None) for success and (URL, Exception) for error.

        """
        filename = url.split("/")[-1]
        filename_full = destination.joinpath(filename)

        if not exists(filename_full) or overwrite:

            filename_temp = filename_full.with_suffix(".modapedl")

            try:

                with session.get(url, stream=True, allow_redirects=True) as response:
                    response.raise_for_status()
                    with open(filename_temp, "wb") as openfile:
                        shutil.copyfileobj(response.raw, openfile, length=16*1024*1024)

                if check:

                    with session.get(url + ".xml", allow_redirects=True) as hdfxml:
                        if hdfxml.status_code == 404:
                            with session.get(url[:-4] + ".cmr.xml", allow_redirects=True) as cmrxml:
                                cmrxml.raise_for_status()
                                file_metadata = self._parse_cmrxml(cmrxml, url.split("/")[-1])
                        else:
                            hdfxml.raise_for_status()
                            file_metadata = self._parse_hdfxml(hdfxml)

                    # check filesize
                    assert str(filename_temp.stat().st_size).strip() == file_metadata["FileSize"], \
                        f'Size: {filename_temp.stat().st_size} != {file_metadata["FileSize"]}'
                    with open(filename_temp, "rb") as openfile:
                        if file_metadata["ChecksumType"] == "CKSUM":
                            checksum = str(cksum(openfile))
                        elif file_metadata["ChecksumType"] == "MD5":
                            md5_hash = hashlib.md5()
                            chunk = openfile.read(65536)
                            while chunk:
                                md5_hash.update(chunk)
                                chunk = openfile.read(65536)
                            checksum = md5_hash.hexdigest().lower()
                        elif file_metadata["ChecksumType"] == "SHA256":
                            sha256_hash = hashlib.sha256()
                            chunk = openfile.read(65536)
                            while chunk:
                                sha256_hash.update(chunk)
                                chunk = openfile.read(65536)
                            checksum = sha256_hash.hexdigest().lower()
                        else:
                            raise ValueError(f'Unknown Checksum Type: {file_metadata["ChecksumType"]}')
                    # check checksum
                    assert checksum == file_metadata["Checksum"], f'Hash: {checksum} != {file_metadata["Checksum"]}'

                shutil.move(filename_temp, filename_full)

            except (HTTPError, ConnectionError, AssertionError, FileNotFoundError) as e:
                try:
                    filename_temp.unlink()
                except FileNotFoundError:
                    pass
                return (filename, e)
        else:
            log.info("%s exists in target. Please set overwrite to True.", filename_full)

        return (filename, None)

    def download(self,
                 targetdir: Path,
                 username: str,
                 password: str,
                 overwrite: bool = False,
                 multithread: bool = False,
                 nthreads: int = 4,
                 max_retries: int = -1,
                 robust: bool = False,
                ) -> List:
        """Download MODIS HDF files.

        This method downloads the MODIS HDF files contained in the
        server response to the `search` call. This requires
        NASA earthdata credentials. To speed up download, multiple
        threads can be used.

        Args:
            targetdir (Path): Target directory for files being downloaded.
            username (str): Earthdata username.
            password (str): Earthdata password.
            overwrite (bool): Replace existing files.
            multithread (bool): Use multiple threads for downloading.
            nthreads (int): Number of threads.
            max_retries (int): Maximum number of retries for failed downloads (for no max, set -1).
            robust (bool): Perform robust downloading (checks file size and checksum).

        Raises:
            DownloadError: If one or more errors were encountered during downloading.
        Returns:
            List of downloaded MODIS HDF filenames.
        """

        # make sure target directory is dir and exists
        assert targetdir.is_dir()
        assert targetdir.exists()

        retry_count = 0
        to_download = self.results.copy()
        downloaded = []

        while True:

            with SessionWithHeaderRedirection(username, password) as session:

                backoff = min(450, 2**retry_count)

                retries = Retry(total=5, backoff_factor=backoff, status_forcelist=[502, 503, 504])
                session.mount(
                    "https://",
                    HTTPAdapter(pool_connections=nthreads, pool_maxsize=nthreads*2, max_retries=retries)
                )

                if multithread:
                    log.debug("Multithreaded download using %s threads. Warming up connection pool.", nthreads)
                    # warm up pool
                    _ = session.get(list(to_download.values())[0]["link"], stream=True, allow_redirects=True)

                    with ThreadPoolExecutor(nthreads) as executor:

                        futures = [executor.submit(self._fetch, session, values["link"], targetdir, overwrite, robust)
                                   for key, values in to_download.items()]

                    downloaded_temp = [x.result() for x in futures]

                else:
                    log.debug("Serial download")
                    downloaded_temp = []

                    for _, values in to_download.items():

                        downloaded_temp.append(
                            self._fetch(session, values["link"], targetdir, overwrite, robust)
                        )

            # check if downloads are OK
            for fid, err in downloaded_temp:
                if err is None:
                    try:
                        del to_download[fid]
                    except KeyError:
                        del to_download[fid[:-4]]
                    downloaded.append(fid)

            if to_download:
                if retry_count < max_retries or max_retries == -1:
                    retry_count += 1
                    log.debug("Retrying downloads! Files left: %s", len(to_download))
                    if max_retries > 0:
                        log.debug("Try %s of %s", retry_count, max_retries)
                    continue

                raise DownloadError(list(to_download.keys()))

            break

        return downloaded
