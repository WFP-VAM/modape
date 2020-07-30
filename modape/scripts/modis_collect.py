#!/usr/bin/env python
# pylint: disable=broad-except
"""modis_collect.py: Collect raw MODIS data into HDF5 file."""
from concurrent.futures import ProcessPoolExecutor, wait
import logging
import multiprocessing as mp
from pathlib import Path
import re
import sys

import click
from modape.constants import REGEX_PATTERNS
from modape.modis import ModisRawH5

# TODO: Increase verbosity
# TODO: Implement checks requested by RM

@click.command()
@click.argument("src_dir", type=click.Path(dir_okay=True, resolve_path=True))
@click.option("-d", "--targetdir", type=click.Path(dir_okay=True, writable=True, resolve_path=True),
              help="Destination for raw HDF5 files")
@click.option('-x', '--compression', type=click.STRING, default='gzip', help='Compression for HDF5 files')
@click.option('--vam-code', type=click.STRING, help='VAM code for dataset to process')
@click.option('--interleave', is_flag=True, help='Interleave MOD13 & MYD13 products to MXD (only works for VIM!)')
@click.option('--parallel-tiles', type=click.INT, default=1, help='Number of tiles processed in parallel (default = None)')
def cli(src_dir: str,
        targetdir: str,
        compression: str,
        vam_code: str,
        interleave: bool,
        parallel_tiles: int
        ) -> None:
    """Collect raw MODIS hdf files into a raw MODIS HDF5 file.

    All MODIS hdf files within srcdir will be collected into a raw MODIS HDF5 file, corresponding to product type and tile (if not global).
    If the respective HDF5 file does not exists in the target directory, it will be created. Otherwhise, the file will be
    updated and the data will be appended.

    16-day MOD13* and MYD13* products can be interleaved into an 8-day product with the new product ID MXD* by adding the `--interleave` flag.

    Args:
        src_dir (str): Path to source directory with HDF files.
        targetdir (str): Targetdir for HDF5 file(s).
        compression (str): Compression filter for HDF5 file.
        vam_code (str): VAM product code to process.
        interleave (bool): Interleave 16-day NDVI/EVI products to 8-day.
        parallel_tiles (int): Process tiles in parallel (number can't exceed ncores - 1).
    """

    log = logging.getLogger(__name__)

    input_dir = Path(src_dir)
    assert input_dir.exists(), "Source directory (src_dir) does not exist!"
    assert input_dir.is_dir(), "Source directory (src_dir) is not a directory!"

    if targetdir is None:
        targetdir = input_dir
    else:
        targetdir = Path(targetdir)

    assert targetdir.is_dir(), "Target directory (targetdir) not a direcory!"
    targetdir.mkdir(exist_ok=True)
    assert targetdir.exists(), "Target directory (targetdir) doesn't exist!"

    hdf_files = list(input_dir.glob('*hdf'))

    if not hdf_files:
        raise ValueError(f"NO HDF files found in src_dir {src_dir}!")

    products = []
    tiles = []
    versions = []

    # Seperate input files into group
    for file in hdf_files:
        products.append(*REGEX_PATTERNS["product"].findall(file.name))
        tiles.append(*REGEX_PATTERNS["tile"].findall(file.name))
        versions.append(*REGEX_PATTERNS["version"].findall(file.name))

    groups = [".*".join(x) for x in zip(products, tiles, versions)]
    groups = list({re.sub('(M.{1})(D.+)', 'M.'+'\\2', x) if REGEX_PATTERNS["VIM"].match(x) else x for x in groups}) # Join MOD13/MYD13

    log.debug("Parsed groups: %s", groups)

    processing_dict = {}

    for group in groups:
        group_pattern = re.compile(group + '.*hdf')

        parameters = dict(
            targetdir=targetdir,
            files=[str(x) for x in hdf_files if group_pattern.match(x.name)],
            interleave=interleave,
            group_id=group.replace('M.', 'MX'),
            compression=compression,
            vam_product_codes=vam_code
        )

        processing_dict.update({group: parameters})

    log.debug("Start processing!")

    if parallel_tiles > 1:
        log.debug("Processing %s parallel tiles!", parallel_tiles)
        available_cores = mp.cpu_count() - 1

        if parallel_tiles > available_cores:
            log.warning("Number of parallel tiles bigger than CPUs available. Scaling down.")
            parallel_tiles = available_cores

        futures = []

        with ProcessPoolExecutor(parallel_tiles) as executor:

            for group, parameters in processing_dict.items():

                log.debug("Submitting %s", group)

                futures.append(
                    executor.submit(
                        _worker,
                        parameters["files"],
                        parameters["targetdir"],
                        parameters["vam_product_code"],
                        parameters["interleave"],
                        parameters["compression"],
                    )
                )

            _ = wait(futures)

        for future in futures:
            assert future.done()
            assert future.exception() is None, f"Received exception {future.exception()}"

    else:

        log.debug("Processing groups sequentially")

        for group, parameters in processing_dict.items():

            log.debug("Processing %s", group)

            _ = _worker(
                parameters["files"],
                parameters["targetdir"],
                parameters["vam_product_code"],
                parameters["interleave"],
                parameters["compression"],
            )

def _worker(files, targetdir, vam_code, interleave, compression):

    raw_h5 = ModisRawH5(
        files=files,
        targetdir=targetdir,
        vam_product_code=vam_code,
        interleave=interleave,
    )

    if not raw_h5.exists:
        raw_h5.create(compression=compression)

    raw_h5.update()

    return str(raw_h5.filename)

def cli_wrap():
    """Wrapper for cli"""

    if len(sys.argv) == 1:
        cli.main(['--help'])
    else:
        cli() #pylint: disable=E1120

if __name__ == '__main__':
    cli_wrap()
