#!/usr/bin/env python3
# pylint: disable=E0401
"""modis_smooth.py: Smooth raw MODIS HDF5 file."""

from concurrent.futures import ProcessPoolExecutor, wait
import logging
import multiprocessing as mp
from pathlib import Path
import sys
from typing import Tuple

import click
import numpy as np
from modape.exceptions import SgridNotInitializedError
from modape.modis import ModisSmoothH5

log = logging.getLogger(__name__)

@click.command(name="Smooth, gapfill and interpolate processed raw MODIS HDF5 files")
@click.argument("src", type=click.Path(dir_okay=True, resolve_path=True))
@click.option("-d", "--targetdir",
              type=click.Path(dir_okay=True, resolve_path=True),
              help='Target directory for smoothed output'
             )
@click.option("-s", "--svalue", type=click.FLOAT, help="S value for smoothing (in log10)")
@click.option("-S", "--srange", type=click.FLOAT, nargs=3, help="S value range for V-curve (float log10(s) values as smin smax sstep - default -1 1 0)")
@click.option("-p", "--pvalue", type=click.FLOAT, help="P value for asymmetric smoothing")
@click.option("-t", "--tempint", type=click.INT, help="Value for temporal interpolation (integer required - default is native temporal resolution i.e. no interpolation)")
@click.option("--tempint-start", type=click.DateTime(formats=["%Y-%m-%d", "%Y%j"]), help="Startdate for temporal interpolation")
@click.option("-n", "--nsmooth", type=click.INT, default=0, help="Number of raw timesteps used for smoothing")
@click.option("-u", "--nupdate", type=click.INT, default=0, help="Number of smoothed timesteps to be updated in HDF5 file")
@click.option("--soptimize", is_flag=True, help="Use V-curve for s value optimization")
@click.option("--parallel-tiles", type=click.INT, help="Number of tiles processed in parallel", default=1)
@click.option('--last-collected', type=click.DateTime(formats=['%Y%j']), help='Last collected date in julian format (YYYYDDD - %Y%j)')
def cli(src: str,
        targetdir: str,
        svalue: float,
        srange: Tuple[float],
        pvalue: float,
        tempint: int,
        tempint_start: str,
        nsmooth: int,
        nupdate: int,
        soptimize: bool,
        parallel_tiles: int,
        last_collected: str,
        ) -> None:
    """Smooth, gapfill and interpolate processed raw MODIS HDF5 files.

        The smoothing function takes a previously created raw MODIS HDF file (as created by modis_collect) as input.
        The raw data can be smoothed with eiter a fixed s value, a pixel-by-pixel s value read from a previously computed grid or
        V-curve optimization of s (creates or updates the s-grid)

        The desired temporal resolution of the ouput file can be defined with tempint, allowing for seamless interpolation.
        Alternatively, tempint can be specified together with tempint_start, returning a timeseries with given frequency
        from the start.

        If a smooth MODIS HDF5 file for a given product, tile (if not global) and temporal interpolation is already in
        the targetdir, it will be updated.

        By default, the entire temporal range of the raw data is used for smoothing, and the entire smoothed data is updated.
        The parameters nsmooth and nupdate can modify this behaviour.

        To speed up processing time, HDF5s can be processed in parallel. By default, all HDF5 files are processed
        sequentially.

    Args:
        src (str): Source (either HDF5 file or directory).
        targetdir (str): Target directory for smooth HDF5 file.
        svalue (float): Smoothing value for fixed smoothing.
        srange (Tuple[float]): Srange for V-Curve optimization as (smin, smax, sstep).
        pvalue (float): Pvalue for asymmetric smoothing.
        tempint (int): Timestep for temporal interpolation.
        tempint_start (str): Start date for custom temporal interpolation.
        nsmooth (int): Number of raw timesteps used for smoothing".
        nupdate (int): Number of smoothed timesteps to be updated in HDF5 file.
        soptimize (bool): Flag for V-Curve optimization of S.
        parallel_tiles (int): Number of paralell HDF5s being processed.
        last_collected (str): Last collected rawdate on which smoothing is performed on.

    """

    input_raw = Path(src)

    if not input_raw.exists():
        msg = "SRC does not exist."
        log.error(msg)
        raise ValueError(msg)

    if input_raw.is_dir():
        files = [str(x) for x in input_raw.glob("*.h5")]
    else:
        files = [str(input_raw)]

    if not files:
        msg = "No files found to process"
        log.error(msg)
        raise ValueError(msg)

    if (nsmooth != 0) and (nupdate != 0):
        if nsmooth < nupdate:
            raise ValueError('nsmooth must be bigger or equal (>=) to nupdate!')

    if targetdir is None:
        if input_raw.is_dir():
            targetdir = input_raw
        else:
            targetdir = input_raw.parent
    else:
        targetdir = Path(targetdir)

    if not targetdir.exists():
        targetdir.mkdir()

    if not targetdir.is_dir():
        msg = "Target directory needs to be a valid path!"
        log.error(msg)
        raise ValueError(msg)

    if tempint_start is not None:
        tempint_start = tempint_start.strftime("%Y%j")

    if last_collected is not None:
        last_collected = last_collected.strftime("%Y%j")

    click.echo("Starting MODIS SMOOTH!")

    if len(srange) != 0:
        assert len(srange) == 3, "Expected 3 values for s-range"

        srange = np.arange(srange[0],
                           srange[1] + srange[2],
                           srange[2],
                           ).round(2)
    else:
        srange = None

    smoothing_parameters = dict(
        svalue=svalue,
        srange=srange,
        p=pvalue,
        nsmooth=nsmooth,
        nupdate=nupdate,
        soptimize=soptimize,
    )

    if parallel_tiles > 1:
        log.debug("Processing %s parallel tiles!", parallel_tiles)
        available_cores = mp.cpu_count() - 1

        if parallel_tiles > available_cores:
            log.warning("Number of parallel tiles bigger than CPUs available!")

        futures = []

        with ProcessPoolExecutor(parallel_tiles) as executor:

            for rawfile in files:

                log.debug("Submitting %s", rawfile)

                futures.append(
                    executor.submit(
                        _worker,
                        rawfile,
                        targetdir,
                        tempint,
                        tempint_start,
                        last_collected,
                        **smoothing_parameters
                    )
                )

            _ = wait(futures)

            for future in futures:
                assert future.done()
                assert future.exception() is None, f"Received exception {future.exception()}"
                assert future.result()

    else:
        log.debug("Processing files sequentially")
        for rawfile in files:
            log.debug("Processing %s", rawfile)

            result = _worker(
                rawfile,
                targetdir,
                tempint,
                tempint_start,
                last_collected,
                **smoothing_parameters
            )

            assert result

    click.echo("MODIS SMOOTH completed!")

def _worker(rawfile: str,
            targetdir: str,
            tempint: int,
            tempint_start: str,
            last_collected: str,
            **kwargs: dict):

    smt_h5 = ModisSmoothH5(
        rawfile=str(rawfile),
        targetdir=str(targetdir),
        tempint=tempint,
        startdate=tempint_start
    )

    if not smt_h5.exists:
        if not kwargs["soptimize"] and not kwargs["svalue"]:
            msg = "Smoothing requires Sgrid which has not been initialized. Please run --soptimize first!"
            raise SgridNotInitializedError(msg)
        smt_h5.create()

    if last_collected is not None:
        last_collected_in_rawfile = smt_h5.last_collected

        if not last_collected_in_rawfile:
            raise ValueError(f"No last_collected recorded in {smt_h5.rawfile}")

        assert last_collected == last_collected_in_rawfile, \
         f"Last collected date in {smt_h5.rawfile} is {last_collected_in_rawfile} not {last_collected}"

    smt_h5.smooth(**kwargs)

    return True

def cli_wrap():
    """Wrapper for cli"""

    if len(sys.argv) == 1:
        cli.main(['--help'])
    else:
        cli() #pylint: disable=E1120

if __name__ == '__main__':
    cli_wrap()
