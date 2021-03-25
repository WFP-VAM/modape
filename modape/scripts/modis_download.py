#!/usr/bin/env python
"""modis_download.py: Query and download MODIS HDF files."""
# pylint: disable=E0401
import datetime
import os
import pathlib
import re
import sys
from typing import List

import click
from modape.exceptions import TargetNotEmptyError
from modape.modis import ModisQuery

@click.command()
@click.argument("products", nargs=-1, type=click.STRING)
@click.option("--roi", type=click.STRING, help="Region of interest. Either LON,LAT or xmin,ymin,xmax,ymax")
@click.option("-b", "--begin-date", type=click.DateTime(formats=["%Y-%m-%d"]), help="Start date for query")
@click.option("-e", "--end-date", type=click.DateTime(formats=["%Y-%m-%d"]), help="End date for query")
@click.option("-d", "--targetdir", type=click.Path(dir_okay=True, writable=True, resolve_path=True),
              help="Destination directory for downloaded files")
@click.option("--target-empty", is_flag=True, help="Fail if there are hdf files in the target directory")
@click.option("--tile-filter", type=click.STRING, help="Filter tiles - supplied as csv list")
@click.option("--username", type=click.STRING, help="Earthdata username")
@click.option("--password", type=click.STRING, help="Earthdata password")
@click.option("--match-begin", is_flag=True, help="Don't allow files with timestamps outside of provided date(s)")
@click.option("--print-results", is_flag=True, help="Print results to console")
@click.option("--download", is_flag=True, help="Download data")
@click.option("--overwrite", is_flag=True, help="Overwrite existing files")
@click.option("--robust", is_flag=True, help="Perform robust download")
@click.option("--max-retries", type=click.INT, help="Max number of retries for downloading", default=-1)
@click.option("--multithread", is_flag=True, help="Use multiple threads for downloading")
@click.option("--nthreads", type=click.INT, help="Number of threads to use", default=4)
@click.option("-c", "--collection", type=click.STRING, default="006", help="MODIS collection")
def cli(products: List[str],
        begin_date: datetime.datetime,
        end_date: datetime.datetime,
        targetdir: str,
        roi: str,
        target_empty: bool,
        tile_filter: str,
        username: str,
        password: str,
        match_begin: bool,
        print_results: bool,
        download: bool,
        overwrite: bool,
        robust: bool,
        max_retries: int,
        multithread: bool,
        nthreads: int,
        collection: str,
        ) -> List:
    """Query and download MODIS products.

    This function allows for querying and downloading MODIS products in bulk.
    Multiple products can be queried and downloaded with one
    function call. For downloading data, valid earthdata credentials are required
    (to register, visit https://urs.earthdata.nasa.gov/users/new).

    To query for both MODIS AQUA and TERRA, replace MOD/MYD with M?D.
    Product IDs also accepted in lowercase.

    Args:
        products (List[str]): List of MODIS product codes to download.
        begin_date (datetime.date): Start date for query / download.
        end_date (datetime.date): End date for query / download.
        targetdir (pathlib.Path): Target directory.
        roi (str): Region of interest (point or bbox as csv string).
        tile_filter (str): MODIS tile filter (as csv string of tile IDs).
        username (str): Earthdata username.
        password (str): Earthdata password.
        match_begin (bool): Match native MODIS timestamp.
        download (bool): Download data.
        overwrite (bool): Replace existing.
        robust (bool): Perform robust downloading (checks file size and checksum).
        max_retries (int): Maximum number of retries for failed downloads (default is -1, infinite).
        multithread (bool): Use multiple threads for downloading.
        nthreads (int): Number of threads for multithread.
        collection (str): MODIS collection version.

    Returns:
        List of downloaded HDF filenames (if overwrite is False, also skipped downloads are included if already existing in targetdir)
    """

    click.echo("\nSTART download_modis.py!")

    if download:

        if username is None or password is None:
            raise ValueError("Download was requested, but credentials are missing.")

    if targetdir is None:
        targetdir = pathlib.Path(os.getcwd())
    else:
        targetdir = pathlib.Path(targetdir)

    # handle targetdir
    targetdir.mkdir(exist_ok=True)

    assert targetdir.exists()
    assert targetdir.is_dir()

    if target_empty:
        for _ in targetdir.glob("M*hdf"):
            try:
                raise TargetNotEmptyError("Found HDF files in target directory with flag --target-empty set!")
            except StopIteration:
                pass

    products_parsed = []

    for product_code in products:
        if "M?D" in product_code.upper():
            products_parsed.append(product_code.replace("?", "O").upper())
            products_parsed.append(product_code.replace("?", "Y").upper())
        else:
            products_parsed.append(product_code.upper())

    if roi is not None:

        try:
            coords = list(map(float, roi.split(',')))
        except ValueError:
            click.echo("Error parsing ROI coordinates. Expected numeric!")

        assert len(coords) == 2 or len(coords) == 4, "Only point or bbox as roi allowed!"

        roi = tuple(coords)

    # parse tile filter and check if OK
    if tile_filter is not None:
        tile_regxp = re.compile(r'^h\d{2}v\d{2}$')
        tiles = []

        for tile_sel in tile_filter.split(','):
            assert re.match(tile_regxp, tile_sel.lower())
            tiles.append(tile_sel.lower())

        tile_filter = tiles

    click.echo('Querying NASA CMR ...')

    query = ModisQuery(
        products=products_parsed,
        aoi=roi,
        begindate=begin_date,
        enddate=end_date,
        tile_filter=tile_filter,
        version=collection,
    )

    query.search(match_begin=match_begin)

    if query.nresults == 0:
        click.echo("No results found! Please check query or make sure CMR is available / reachable.")
        return []

    click.echo(f'Done! Found {query.nresults} results!')

    if print_results:
        click.echo("\n")
        for key, values in query.results.items():
            click.secho(key, bold=True)
            click.echo(values)
            click.echo("\n")

    downloaded = []

    if download:

        click.echo('Downloading!')

        if query.nresults > 0:

            downloaded = query.download(
                targetdir=targetdir,
                username=username,
                password=password,
                overwrite=overwrite,
                multithread=multithread,
                nthreads=nthreads,
                robust=robust,
                max_retries=max_retries,
            )

    click.echo('modis_download.py COMPLETED! Bye! \n')
    return downloaded

def cli_wrap():
    """Wrapper for cli"""

    if len(sys.argv) == 1:
        cli.main(['--help'])
    else:
        cli() #pylint: disable=E1120

if __name__ == '__main__':
    cli_wrap()
