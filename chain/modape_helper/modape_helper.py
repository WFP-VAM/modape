import os
import h5py
import re

try:
    from pathlib2 import Path
except ImportError:
    from pathlib import Path
from datetime import datetime
from dateutil.relativedelta import relativedelta
from modape.utils import fromjulian
from .timeslicing import ModisInterleavedOctad


def get_date_from_raw_h5_modis_tile(f, idx):
    with h5py.File(f, 'r') as h5f:
        dates = h5f.get('dates')
        return fromjulian(dates[idx].decode())


def get_first_date_in_raw_h5_modis_tile(f):
    return get_date_from_raw_h5_modis_tile(f, 0)


def get_last_date_in_raw_h5_modis_tile(f):
    return get_date_from_raw_h5_modis_tile(f, -1)


def get_first_date_in_raw_modis_tiles(folder):
    files = Path(folder).glob('*.h5')
    first_dates = []
    for f in files:
        first_dates.append(get_first_date_in_raw_h5_modis_tile(f))
    if len(first_dates) > 0:
        return max(first_dates)
    else:
        return None


def get_last_date_in_raw_modis_tiles(folder):
    files = Path(folder).glob('*.h5')
    last_dates = []
    for f in files:
        last_dates.append(get_last_date_in_raw_h5_modis_tile(f))
    if len(last_dates) > 0:
        return min(last_dates)
    else:
        return None


def has_collected_dates(h5_file, dates):

    r = False
    with h5py.File(h5_file, "r") as h5_open:
        collected_dates = [d.decode() for d in h5_open.get("dates")]
        r = dates == collected_dates
        if not r:
            for d in dates:
                if d not in collected_dates:
                    print("Missing date: {}".format(d))
    return r


def curate_downloads(folder, tiles, begin_date, end_date) -> bool:
    """
    :param folder: download folder with MODUS HDFs
    :param tiles: list of tile-IDs, like: [r"h21v08", r"h22v08", r"h21v09", r"h22v09"]
    :param begin_date: first date that should be available in HDFs
    :param end_date: last date that should be available in HDFs
    :return: bool
    """
    re_hdf_slice_date = re.compile(r'^.+A(\d{7})\.h\d{2}v\d{2}\.\d{3}\.\d{13}\.hdf$')
    re_hdf_production_date = re.compile(r'^.+(\d{13})\.hdf$')
    re_hdf_tile = re.compile(r'^.+(h\d{2}v\d{2})\.\d{3}\.\d{13}\.hdf$')

    # Slices we're looking for:
    time_slices = [ModisInterleavedOctad(begin_date)]
    next_slice = time_slices[-1].next()
    while next_slice.startsBeforeDate(end_date + relativedelta(days=1)):
        time_slices.append(next_slice)
        next_slice = time_slices[-1].next()
    time_slices = set([str(s) for s in time_slices])

    # Locate HDFs:
    hdf_files = [x.as_posix() for x in Path(folder).glob('*.hdf')]
    if len(hdf_files) < 1:
        return False

    # Dates for HDFs:
    hdf_dates = set(filter(lambda dte: re.compile(r'^\d{7}$').match(dte),
                          [re.sub(re_hdf_slice_date, r'\1', x) for x in hdf_files]))

    for sDte in sorted(hdf_dates):
        re_cur_dte = re.compile(r'^.+A' + sDte + r'\.h\d{2}v\d{2}\.\d{3}\.\d{13}\.hdf$')
        # Create a dictionary for the tiles we are looking for:
        tile2file = {tile: '' for tile in set(tiles)}
        for file in filter(lambda f: re_cur_dte.match(f), hdf_files):
            if not sDte in time_slices:
                print('We do not want this date: {}; removing this file: {}'.format(sDte, file))
                os.remove(file)
            else:
                hdfTile = re.sub(re_hdf_tile, r'\1', file)
                if hdfTile not in tile2file:
                    print('We do not want this tile: {}; removing this file: {}'.format(hdfTile, file))
                    os.remove(file)
                else:
                    if len(tile2file[hdfTile]) == 0:
                        # first file we are hitting for this tile on this date; remember it:
                        tile2file[hdfTile] = file
                    elif int(re.sub(re_hdf_production_date, r'\1', file)) \
                            > int(re.sub(re_hdf_production_date, r'\1', tile2file[hdfTile])):
                        print('Removing this older file: {}'.format(tile2file[hdfTile]))
                        os.remove(tile2file[hdfTile])
                        tile2file[hdfTile] = file
                    else:
                        print('Removing this older file: {}'.format(file))
                        os.remove(file)
        if sDte in time_slices:
            missing = [tile for tile, file in tile2file.items() if len(file) == 0]
            if len(missing):
                print('Missing tile(s) for {}: {}'.format(sDte, ', '.join(missing)))
                return False

    return True
