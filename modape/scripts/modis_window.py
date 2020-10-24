#!/usr/bin/env python
# pylint: disable=broad-except,C0103
"""modis_window.py: Create mosaics from smooth MODIS HDF5 files and
   save them as GeoTIFFs.
"""
import datetime
import logging
from pathlib import Path
import re
import sys
from typing import Tuple

import click
from modape.constants import REGEX_PATTERNS
from modape.modis import ModisMosaic

__modis_window = None


def register(f):
    global __modis_window
    __modis_window = f
    return f


def modis_window(**kwargs):
    global __modis_window
    return __modis_window(**kwargs)


log = logging.getLogger(__name__)


@click.command()
@click.argument("src")
@click.option("-d", "--targetdir", type=click.Path(dir_okay=True, resolve_path=True, writable=True), help="Target directory for Tiffs")
@click.option("-b", "--begin-date", type=click.DateTime(formats=["%Y-%m-%d"]), help="Begin date for Tiffs")
@click.option("-e", "--end-date", type=click.DateTime(formats=["%Y-%m-%d"]), help="End date for Tiffs")
@click.option("--roi", type=click.STRING, help="AOI for clipping (as ULX,ULY,LRX,LRY)")
@click.option("--region", type=click.STRING, help="Region prefix for Tiffs (default is reg)", default="reg")
@click.option("--sgrid", is_flag=True, help="Extract sgrid instead of data")
@click.option("--force-doy", is_flag=True, help="Force DOY filenaming")
@click.option("--filter-product", type=click.STRING, help="Filter by product")
@click.option("--filter-vampc", help="Filter by VAM parameter code")
@click.option("--target-srs", help="Target spatial reference for warping", default="EPSG:4326")
@click.option('--co', multiple=True, help="GDAL creationOptions", default=["COMPRESS=LZW", "PREDICTOR=2"])
@click.option("--clip-valid", is_flag=True, help="clip values to valid range for product")
@click.option("--round-int", type=click.INT, help="Round to integer palces (either decimals or exponent of 10)")
@click.option("--gdal-kwarg", type=click.STRING, multiple=True, help="Addition kwargs for GDAL in form KEY=VALUE (multiple allowed)")
@click.option("--overwrite", is_flag=True, help="Overwrite existsing Tiffs")
@register
def cli(src: str,
        targetdir: str,
        begin_date: datetime.date,
        end_date: datetime.date,
        roi: str,
        region: str,
        sgrid: bool,
        force_doy: bool,
        filter_product: str,
        filter_vampc: str,
        target_srs: str,
        co: Tuple[str],
        clip_valid: bool,
        round_int: int,
        gdal_kwarg: Tuple[str],
        overwrite: bool) -> None:
    """Creates GeoTiff Mosaics from HDF5 files.

    The input can be either raw or smoothed HDF5 files. With the latter,
    the S-grid can also be mosaiced using the `--sgrid` flag.
    If no ROI is passed, then the full extent of the input files will be mosaiced, otherwise
    the GeoTiffs will be clipped to the ROI after warping.
    By default, the MODIS data will be warped to WGS1984 (EPSG:4326), but a custom spatial reference
    can be passed in with `--target-srs`, in wich case the target resolution has to be manually defined too.
    If required, the output data can be clipped to the valid data range of the input data using the `--clip-valid` flag.
    Also, the output data can be rounded, if it's float (eg. sgrid) to defined precision, or if its integer to the defined
    exponent of 10 (round_int will be multiplied by -1 and passed to np.round!!!)
    Specific creation options can be passed to gdalwarp and gdaltranslate using the `--co` flag. The flag can be used multiple times,
    each input needs to be in the gdal format for COs, e.g. `KEY=VALUE`.
    Additional options can be passed to gdal.Translate (and with restrictions to warp) using `--gdal-kwarg`,
    e.g. `--gdal-kwarg xRes=10 --gdal-kwarg yRes=10`. The keywords are sensitive to how gdal expects them,
    as they are directly passed to gdal.TranlsateOptions. For details, please check the documentation of gdal.TranslateOptions.

    Args:
        ctx (click.core.Context): Context for kwargs.
        src (str): Input directory (or file).
        targetdir (str): Target directory.
        begin_date (datetime.date): Start date for tiffs.
        end_date (datetime.date): End date for tiffs.
        roi (str): ROI for clipping.
        region (str): Region for filename.
        sgrid (bool): Extract sgrid instead of data.
        force_doy (bool): Force DOY in filename.
        filter_product (str): Filter input by product code.
        filter_vampc (str): Filter inpout by vam parameter code.
        target_srs (str): Target spatial reference (in format GDAL can process).
        co (Tuple[str]): Creation options passed to gdal.Translate.
        clip_valid (bool): Clip data to valid range.
        round_int (int): Round integer.
        overwrite (bool): Overwrite existsing Tiffs.

    """

    src_input = Path(src)

    if not src_input.exists():
        msg = "src_dir does not exist."
        log.error(msg)
        raise ValueError(msg)

    if src_input.is_dir():
        files = list(src_input.glob("*.h5"))
    else:
        files = [src_input]

    if filter_product is not None:
        product = filter_product.upper()
        files = [x for x in files if product in x.name]

    if filter_vampc:
        vampc = filter_vampc.upper()
        files = [x for x in files if vampc in x.name]

    if not files:
        msg = "No files found to process! Please check src and/or adjust filters!"
        log.error(msg)
        raise ValueError(msg)

    groups = [REGEX_PATTERNS["tile"].sub("*", x.name) for x in files]
    group_check = {".".join(x.split(".")[:-2]) for x in groups}
    if len(group_check) > 1:
        raise ValueError("Multiple product groups in input. Please filter or use separate directories!")

    groups = list(set(groups))

    if roi is not None and not isinstance(roi, list):
        roi = [float(x) for x in roi.split(',')]
        if len(roi) != 4:
            raise ValueError("ROI for clip needs to be bounding box in format ULX,ULY,LRX,LRY")

    if targetdir is None:
        if src_input.is_dir():
            targetdir = src_input
        else:
            targetdir = src_input.parent
    else:
        targetdir = Path(targetdir)

    if not targetdir.exists():
        targetdir.mkdir()

    if not targetdir.is_dir():
        msg = "Target directory needs to be a valid path!"
        log.error(msg)
        raise ValueError(msg)

    if begin_date:
        begin_date = begin_date.date()

    if end_date:
        end_date = end_date.date()

    if sgrid:
        dataset = "sgrid"
        clip_valid = False
    else:
        dataset = "data"

    if round_int is not None:
        round_int = round_int * -1

    gdal_kwargs = {}
    if gdal_kwarg:
        gdal_kwargs.update(
            {key:value for x in gdal_kwarg for key, value in [x.split("=")]}
        )

    click.echo("\nSTARTING modis_window.py!")

    mosaics = []
    for group in groups:
        log.debug("Processing group %s", group)

        group_pattern = re.compile(group)
        group_files = [str(x) for x in files if group_pattern.match(x.name)]

        mosaic = ModisMosaic(group_files)

        mosaics.extend(
            mosaic.generate_mosaics(
                dataset=dataset,
                targetdir=targetdir,
                target_srs=target_srs,
                aoi=roi,
                overwrite=overwrite,
                force_doy=force_doy,
                prefix=region,
                start=begin_date,
                stop=end_date,
                clip_valid=clip_valid,
                round_int=round_int,
                creationOptions=list(co),
                **gdal_kwarg
            )
        )

    click.echo("\nCOMPLETED modis_window.py!")
    return mosaics


if __name__ == '__main__':
    cli()
