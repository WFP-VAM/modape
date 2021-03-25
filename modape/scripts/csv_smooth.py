#!/usr/bin/env python
"""csv_smooth.py: Smooth timeseries in a CSV file."""
# pylint: disable=E0401,C0103
from array import array
from pathlib import Path
import sys
from typing import Tuple

import click
from modape.whittaker import ws2d, ws2doptv, ws2doptvp # pylint: disable=E0611
import numpy as np
import pandas as pd


@click.command()
@click.argument("csv_file", type=click.Path(exists=True, file_okay=True, resolve_path=True))
@click.option('-s', '--svalue', type=click.FLOAT, help='S value for smoothing (has to be log10(s))')
@click.option("-S", "--srange", type=click.FLOAT, nargs=3, help="S value range for V-curve (float log10(s) values as smin smax sstep - default 0 4 0.1)")
@click.option("-p", "--pvalue", type=click.FLOAT, help="P value for asymmetric smoothing")
@click.option('-n', '--nodata', type=click.FLOAT, help="nodata value", default=0)
@click.option('--skip-header-rows', type=click.INT, help="Number of header rows to skip", default=0)
@click.option('--data-start', type=click.INT, help="Number of row with first datapoint", default=0)
@click.option('--index-column', type=click.INT, help="Index column number", default=None)
@click.option("--soptimize", is_flag=True, help="Use V-curve for s value optimization")
def cli(csv_file: str,
        svalue: float,
        srange: Tuple[float],
        pvalue: float,
        nodata: float,
        skip_header_rows: int,
        data_start: int,
        index_column: int,
        soptimize: bool) -> None:
    """Smooth timeseries in a CSV file.

    Smoothing of CSV file with one timeseries per column. The number of rows skipped (header) can be controled with the `--skip-header-rows` flag.
    Additionally, the `--data-start` flag sets the index for the first data value in each column (0 by default).
    If the csv includes and index column, it can be set using the `--index-column` flag.

    By default, the nodata value will be assumed to be 0. This is crucial to be set correctly, as the smoothing weights
    are based on it.

    The defined or determined sopt and log10(sopt) values are appended to the end of the smoothed timeseries.

    The output filename contains a suffix which indicates the smoothing method used:

        - fixed s: _filt0.csv\n
        - V-curve: _filtoptv.csv\n
        - asymmetric V-curve: _filtoptvp.csv\n

    The resulting CSV is created in the directory the input file is located.

    Args:
        csv_file (str): Path to input csv.
        svalue (float): log10 s value for fixed smoothing.
        srange (Tuple[float]): Srange for V-Curve optimization as (smin, smax, sstep).
        pvalue (float):  Pvalue for asymmetric smoothing.
        nodata (float): Nodata value.
        skip_header_rows (int): Number of rows to be skipped when reading.
        data_start (int): Index of first datapoint in column.
        index_column (int): Index for index column.
        soptimize (bool): Apply V-curve optimization.

    """

    input_csv = Path(csv_file)

    if input_csv.name[-3:] != "csv":
        raise ValueError("Input needs to be CSV file!")

    if not soptimize and svalue is None:
        raise ValueError("Need either soptimize or svalue!")

    click.echo('\nStarting smoothCSV.py ... \n')

    df = pd.read_csv(str(input_csv), header=skip_header_rows, index_col=index_column) # Read input
    resdf = df.copy()
    resdf = resdf.append(
        pd.Series(name='Sopt', dtype="float32")
        ).append(
            pd.Series(name='logSopt', dtype="float32")
        )

    if soptimize:

        if srange:
            assert len(srange) == 3, "Expected 3 values for s-range"

            srange = np.arange(srange[0],
                               srange[1] + srange[2],
                               srange[2],
                               ).round(2)
        else:
            srange = np.arange(0, 4.1, 0.1).round(2)

        smin, smax, sstep = (srange[0], srange[-1], srange[0]-srange[1])

        if pvalue:
            outname = str(input_csv).replace(".csv", "_filtoptvp.csv")
            resdf = resdf.append(pd.Series(name='pvalue', dtype="float32"))
            msg = f"Smoothing using asymmetric V-curve optimization with smin:{smin}, smax:{smax}, sstep:{sstep} and pvalue:{pvalue}.\n\nWriting to file: {outname}\n"

        else:
            outname = str(input_csv).replace(".csv", "_filtoptv.csv")
            msg = f"\nSmoothing using V-curve optimization with smin:{smin}, smax:{smax}, sstep:{sstep}.\n\nWriting to file: {outname}\n"

        click.echo(msg)

        srange = array("d", srange)

        for col in df.columns:
            y = np.array(df[col].values[data_start:], dtype="double")
            w = np.array((y != nodata)*1, dtype='double')
            if pvalue:
                z, sopt = ws2doptvp(y, w, srange, pvalue)
            else:
                z, sopt = ws2doptv(y, w, srange)

            result = pd.concat([pd.Series(z, dtype='float32'),
                                pd.Series([sopt, np.log10(sopt)])],
                               ignore_index=True)

            if pvalue:
                result = pd.concat([result, pd.Series(pvalue)])

            resdf[col] = result.values

    else:
        outname = str(input_csv).replace(".csv", "_filt0.csv")
        s = 10**svalue

        click.echo(f'\nSmoothing using fixed S value {s}. Writing to file: {outname}\n')

        # Iterate columns (skip 1st)
        for col in df.columns:
            y = np.array(df[col].values[data_start:], dtype="double")
            w = np.array((y != nodata)*1, dtype='double')
            z = ws2d(y, s, w)
            result = pd.concat([pd.Series(z, dtype='float32'),
                                pd.Series([s, np.log10(s)])],
                               ignore_index=True)
            resdf[col] = result.values

    # Write to disk
    resdf.to_csv(outname)
    click.echo('\ncsv_smooth.py finished.\n')

def cli_wrap():
    """Wrapper for cli"""

    if len(sys.argv) == 1:
        cli.main(['--help'])
    else:
        cli() #pylint: disable=E1120

if __name__ == '__main__':
    cli_wrap()
