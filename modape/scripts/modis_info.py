#!/usr/bin/env python
"""modis_info.py: Return metadata stored in MODIS HDF5 files."""
#pylint: disable=E0401
import sys
import time

import click
import h5py

@click.command()
@click.argument("file", type=click.Path(exists=True, file_okay=True, resolve_path=True))
def cli(file):
    """Info tool for processed MODIS HDF5 files.

    Returns metadata on processed MODIS HDF5 files, both for raw and smoothed files.
    """

    # Read metadata
    try:
        with h5py.File(file, 'r') as h5f:
            dset = h5f.get('data')
            dates = h5f.get('dates')
            rawdates = h5f.get('rawdates')
            dim = dset.shape
            attrs = dict(dset.attrs)

            # get dates
            startdate = dates[0].decode()
            enddate = dates[-1].decode()

            if rawdates:
                last_collected = rawdates[-1].decode()
            else:
                last_collected = enddate

        temporalresolution = attrs['temporalresolution']
        resolution = attrs['resolution']
        nodata_value = attrs['nodata']
        processing_timestamp = attrs['processingtimestamp']
        ncols = attrs['RasterXSize']
        nrows = attrs['RasterYSize']

        try:
            last_run = attrs['lastrun']
            hdf5_type = "smoothed"
        except KeyError:
            last_run = None
            hdf5_type = "raw"

    except:
        raise SystemExit('Error reading file information.')

    message_head = '\n\n MODAPE info tool - {}'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))

    click.secho(message_head, bold=True, fg="red")

    lines = [
        " \n --- \n"
        f"\n File: {file}\n",
        f"\n Type: MODIS {hdf5_type} HDF5\n",
        "\n Dimensions:\n",
        f"\n     -{nrows} rows\n",
        f"\n     -{ncols} columns\n",
        f"\n     -{dim[1]} timesteps\n",
        f"\n Start date: {startdate}\n",
        f"\n End date: {enddate}\n",
        f"\n Last collected: {last_collected}\n",
        f"\n Temporal resolution: {temporalresolution}\n",
        f"\n Spatial resolution: {resolution}\n",
        f"\n NoData value: {nodata_value}\n",
        f"\n Last modified: {processing_timestamp}\n",
    ]

    if last_run is not None:
        lines.append(f"\n Last smoothing run: Whittaker smoother with {last_run}\n")

    click.echo_via_pager(lines)
    click.echo(" --- \n")

def cli_wrap():
    """Wrapper for cli"""

    if len(sys.argv) == 1:
        cli.main(['--help'])
    else:
        cli() #pylint: disable=E1120

if __name__ == '__main__':
    cli_wrap()
