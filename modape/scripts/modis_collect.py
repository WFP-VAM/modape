#!/usr/bin/env python
# pylint: disable=broad-except,E0401
"""modis_collect.py: Collect raw MODIS data into HDF5 file."""
from concurrent.futures import ProcessPoolExecutor, wait
from datetime import datetime
import logging
import multiprocessing as mp
from pathlib import Path
import re
import sys

import click
from modape.constants import REGEX_PATTERNS
from modape.modis import ModisRawH5

@click.command()
@click.argument("src_dir", type=click.Path(dir_okay=True, resolve_path=True))
@click.option("-d", "--targetdir", type=click.Path(dir_okay=True, writable=True, resolve_path=True),
              help="Destination for raw HDF5 files")
@click.option('-x', '--compression', type=click.STRING, default='gzip', help='Compression for HDF5 files')
@click.option('--vam-code', type=click.STRING, help='VAM code for dataset to process')
@click.option('--interleave', is_flag=True, help='Interleave MOD13 & MYD13 products to MXD (only works for VIM!)')
@click.option('--parallel-tiles', type=click.INT, default=1, help='Number of tiles processed in parallel (default = 1)')
@click.option('--cleanup', is_flag=True, help='Remove collected HDF files')
@click.option('--force', is_flag=True, help='Force collect process not failing on corrupt inputs')
@click.option('--last-collected', type=click.DateTime(formats=['%Y%j']), help='Last collected date in julian format (YYYYDDD - %Y%j)')
@click.option("--tiles-required", type=click.STRING, help="Required tiles - supplied as csv list")
def cli(src_dir: str,
        targetdir: str,
        compression: str,
        vam_code: str,
        interleave: bool,
        parallel_tiles: int,
        cleanup: bool,
        force: bool,
        last_collected: datetime,
        tiles_required: str,
        ) -> None:
    """Collect raw MODIS hdf files into a raw MODIS HDF5 file.

    All MODIS HDF files within srcdir will be collected into a raw MODIS HDF5 file, corresponding to product type and tile (if not global).
    If the respective HDF5 file does not exists in the target directory, it will be created. Otherwhise, the file will be
    updated and the data will be appended.

    If an HDF file cannot be opened / read, an IOException will be raised and the collection process fails. If the `--force` flag is used,
    no exception is raised and a NoData fill value is used, keeping the collection process running.

    16-day MOD13* and MYD13* products can be interleaved into an 8-day product with the new product ID MXD* by adding the `--interleave` flag.

    Args:
        src_dir (str): Path to source directory with HDF files.
        targetdir (str): Targetdir for HDF5 file(s).
        compression (str): Compression filter for HDF5 file.
        vam_code (str): VAM product code to process.
        interleave (bool): Interleave 16-day NDVI/EVI products to 8-day.
        parallel_tiles (int): Process tiles in parallel (number can't exceed ncores - 1).
        cleanup (bool): Remove collected HDF files.
        force (bool): When corrupt HDF input is encountered, use nodata instead of raising exception.
        last_collected (datetime.datetime): Julian day of last collected raw timestep.
        tiles_required (str): MODIS tiles required to be in input (as csv string of tile IDs).
    """

    log = logging.getLogger(__name__)

    input_dir = Path(src_dir)
    assert input_dir.exists(), "Source directory (src_dir) does not exist!"
    assert input_dir.is_dir(), "Source directory (src_dir) is not a directory!"

    if targetdir is None:
        targetdir = input_dir
    else:
        targetdir = Path(targetdir)

    targetdir.mkdir(exist_ok=True)
    assert targetdir.exists(), "Target directory (targetdir) doesn't exist!"

    if last_collected is not None:
        last_collected = last_collected.strftime("%Y%j")

    # parse tiles required and check if OK
    if tiles_required is not None:
        tile_regxp = re.compile(r'^h\d{2}v\d{2}$')
        _tiles = []

        for tile_req in tiles_required.split(','):
            assert re.match(tile_regxp, tile_req.lower())
            _tiles.append(tile_req.lower())

        tiles_required = frozenset(_tiles)

    click.echo("Starting MODIS COLLECT!")

    hdf_files = list(input_dir.glob('*hdf'))

    if not hdf_files:
        raise ValueError(f"NO HDF files found in src_dir {src_dir}!")

    log.debug("Found %s hdf files.", len(hdf_files))

    products = []
    tiles = []
    versions = []

    # Seperate input files into group
    for file in hdf_files:

        product = REGEX_PATTERNS["product"].findall(file.name)
        if not product:
            raise ValueError("Could not extract product from Filename!")
        products.append(*product)

        version = REGEX_PATTERNS["version"].findall(file.name)
        if not version:
            raise ValueError("Could not extract version from Filename!")
        versions.append(*version)

        tile = REGEX_PATTERNS["tile"].findall(file.name)
        if not tile:
            tile = ['']
        tiles.append(*tile)

    if tiles_required is not None:
        log.debug("Checking for required tiles ...")

        for product in set(products):
            log.debug("Product: %s", product)
            datetiles = [".".join(x.name.split(".")[1:3]) for x in hdf_files if re.match("^" + product, x.name)]
            product_tiles = {x.split(".")[1] for x in datetiles}

            timestamps = []
            if product_tiles != tiles_required:
                raise ValueError("Tiles for product %s don't match tiles-required!" % product)
            for reqtile in tiles_required:
                timestamps.append({x.split(".")[0] for x in datetiles if reqtile in x})

            iterator = iter(timestamps)
            reference = next(iterator)
            if not all(x == reference for x in iterator):
                raise ValueError("Not all tiles have same timesteps for product %s!" % product)
            log.debug("Timestamps OK for %s", product)

    groups = [".*".join(x) for x in zip(products, tiles, versions)]
    groups = list({re.sub('(M.{1})(D.+)', 'M.'+'\\2', x) if REGEX_PATTERNS["VIM"].match(x) else x for x in groups}) # Join MOD13/MYD13
    groups.sort()
    log.debug("Parsed groups: %s", groups)

    to_process = {}

    collected = []

    for group in groups:
        group_pattern = re.compile(group + '.*hdf')
        group_files = [str(x) for x in hdf_files if group_pattern.match(x.name)]
        group_files.sort()

        _raw_h5 = ModisRawH5(files=group_files,
                            targetdir=targetdir,
                            vam_product_code=vam_code,
                            interleave=interleave)

        if last_collected is not None:

            if not _raw_h5.exists:
                raise ValueError("Output H5 %s does not exist! Can't check last-collected!" % _raw_h5.filename)

            last_collected_infile = _raw_h5.last_collected

            if not last_collected_infile:
                raise ValueError(f"No last_collected recorded in {_raw_h5.filename}")

            if not last_collected == last_collected_infile:
                raise ValueError(f"Last collected date in file is {last_collected_infile} not {last_collected}")

        to_process.update({
            group: {
                "raw_h5": _raw_h5,
                "compression": compression,
                "force": force,
            }
        })

    log.debug("Start processing!")

    if parallel_tiles > 1:
        log.debug("Processing %s parallel tiles!", parallel_tiles)
        available_cores = mp.cpu_count() - 1

        if parallel_tiles > available_cores:
            log.warning("Number of parallel tiles bigger than CPUs available!")

        futures = []

        with ProcessPoolExecutor(parallel_tiles) as executor:

            for group, parameters in to_process.items():

                log.debug("Submitting %s", group)

                futures.append(
                    executor.submit(_worker, **parameters)
                )

            _ = wait(futures)

        for future in futures:
            assert future.done()
            if future.exception() is not None:
                raise future.exception()
            collected.extend(future.result())

    else:

        log.debug("Processing groups sequentially")

        for group, parameters in to_process.items():

            log.debug("Processing %s", group)

            collected.extend(
                _worker(**parameters)
            )

    if cleanup and collected:
        tracefile = targetdir.joinpath(".collected")
        with open(str(tracefile), "a") as tf_open:
            log.debug("Cleaning up collected files")
            for to_remove in collected:
                to_remove_obj = Path(to_remove)
                tf_open.write(to_remove_obj.name + "\n")
                log.debug("Removing %s", to_remove)
                to_remove_obj.unlink()

    click.echo("MODIS COLLECT completed!")

def _worker(raw_h5, compression, force):

    if not raw_h5.exists:
        raw_h5.create(compression=compression)

    collected = raw_h5.update(force=force)

    return collected

def cli_wrap():
    """Wrapper for cli"""

    if len(sys.argv) == 1:
        cli.main(['--help'])
    else:
        cli() #pylint: disable=E1120

if __name__ == '__main__':
    cli_wrap()
