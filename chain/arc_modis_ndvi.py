#!/usr/bin/env python
"""
  arc_modis_ndvi.py: Flask Service app for collecting, processing and disseminating filtered NDVI.
           Production the time series leverages the WFP VAM MODAPE toolkit: https://github.com/WFP-VAM/modape

  Dependencies: arc-modape (1.0), Numpy, ...

  Author: Rob Marjot, (c) ARC 2020

"""

import re
import contextlib
import hashlib
import json
import os
import glob
import shutil
import click
import logging

logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"),
                    format='[%(asctime)s %(levelname)s] (%(name)s:%(lineno)d) - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

_log = logging.getLogger(__name__ + '_echo_through_log')
_log.propagate = False
h = logging.StreamHandler()
h.setLevel(os.environ.get("LOGLEVEL", "INFO"))
h.setFormatter(logging.Formatter('[%(asctime)s ECHO] - %(message)s', '%Y-%m-%d %H:%M:%S'))
_log.addHandler(h)
def echo_through_log(message=None, file=None, nl=True, err=False, color=None):
    _log.info(re.sub('\s+',' ', message).strip())
click.echo = echo_through_log

log = logging.getLogger(__name__)
log.propagate = False
h = logging.StreamHandler()
h.setLevel(os.environ.get("LOGLEVEL", "INFO"))
h.setFormatter(logging.Formatter('[%(asctime)s %(levelname)s] - %(message)s', '%Y-%m-%d %H:%M:%S'))
log.addHandler(h)

from datetime import datetime, date
from dateutil.relativedelta import relativedelta
from flask import Flask, jsonify, send_file
from threading import Thread, Timer
from pathlib import Path

from modape_helper import get_first_date_in_raw_modis_tiles, get_last_date_in_raw_modis_tiles,\
    curate_downloads, has_collected_dates
from modape.scripts.modis_download import cli as modis_download
from modape.scripts.modis_collect import cli as modis_collect
from modape.scripts.modis_smooth import cli as modis_smooth
from modape.scripts.modis_window import cli as modis_window

from modape_helper.timeslicing import Dekad, ModisInterleavedOctad

try:
    from types import SimpleNamespace as Namespace
except ImportError:
    from argparse import Namespace

app_state = None


def generate_file_md5(filepath, blocksize=2 ** 16):
    m = hashlib.md5()
    with open(filepath, "rb") as f:
        while True:
            buf = f.read(blocksize)
            if not buf:
                break
            m.update(buf)
    return m.hexdigest()


def exists_smooth_h5s(tiles, basedir):
    # Check all tiles for a corresponding H5 archive in VIM/SMOOTH:
    return all(
     [Path(os.path.join(basedir, 'VIM', "SMOOTH", "MXD13A2.{}.006.txd.VIM.h5".format(tile))).exists() for tile in tiles]
    )


def app_index():
    global app_state
    if app_state.fetcherThread.is_alive() and not getattr(app_state, 'suspended', False):
        return "Fetcher is running (or suspended), try again later\n", 503
    else:
        files = {}
        for f in sorted(glob.glob(os.path.join(app_state.basedir, 'VIM', 'SMOOTH', 'EXPORT', app_state.file_pattern))):
            if os.path.isfile(f + '.md5'):
                with open(f + '.md5') as mdf:
                    files[os.path.basename(f)] = re.sub('\\s+', '', mdf.readline())
        return jsonify(files)


def app_download(filename):
    global app_state
    if app_state.fetcherThread.is_alive() and not getattr(app_state, 'suspended', False):
        return "Fetcher is running (or suspended), try again later\n", 503
    else:
      try:
        return send_file(os.path.join(app_state.basedir, 'VIM', 'SMOOTH', 'EXPORT', filename), as_attachment=True, mimetype=app_state.mimetype)
      except FileNotFoundError:
        return ('', 404)


def app_fetch():
    global app_state
    if app_state.fetcherThread.is_alive() or getattr(app_state, 'suspended', False):
        return "[{}] Fetcher is already running (or suspended), try again later\n".format(
            datetime.now().strftime('%Y-%m-%d %H:%M:%S')), 503
    else:
        # Check all tiles for a corresponding H5 archive in VIM/SMOOTH:
        if not exists_smooth_h5s(app_state.tile_filter, app_state.basedir):
            app_state.fetcherThread = Timer(5, app_do_init, ())
            app_state.fetcherThread.start()
            return "[{}] Initialisation is scheduled to start (or resume) in 5 seconds...\n".format(
                datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        else:
            app_state.fetcherThread = Timer(5, app_do_processing, ())
            app_state.fetcherThread.start()
            return "[{}] Fetching and processing is scheduled to start in 5 seconds...\n".format(
                datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


def app_suspend():
    global app_state
    app_state.suspended = True
    if app_state.fetcherThread.is_alive():
        return "Fetcher is busy suspending; please check back later to see the suspended state confirmed...\n", 503
    else:
        s = "[{}] Fetcher suspended; restart service to resume production.\n".format(
            datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        log.info(s)
        return s


def app_log(filename):
    global app_state
    return send_file(os.path.join(app_state.basedir, 'log', filename), as_attachment=False, mimetype='text/plain')


def app_do_init():
    global app_state
    do_init(app_state)


def app_do_processing():
    global app_state
    do_processing(app_state)


def do_processing(args, only_one_inc=False):
    # download and ingest:
    while True:
        last_date = get_last_date_in_raw_modis_tiles(os.path.join(args.basedir, 'VIM'))
        next_date = last_date + relativedelta(days=8)
        if last_date.year < next_date.year:
            # handle turning of the year:
            next_date = datetime(next_date.year, 1, 1).date()

        if getattr(args, 'download_only', False) or (
                not getattr(args, 'collect_only', False) and
                not getattr(args, 'smooth_only', False) and
                not getattr(args, 'export_only', False)):

            if next_date > date.today():  # stop after today:
                break

            log.info('Downloading: {}...'.format(next_date))
            downloaded = modis_download.callback(
                products=['M?D13A2'],
                begin_date=datetime.combine(next_date, datetime.min.time()),
                end_date=datetime.combine(next_date, datetime.min.time()),
                targetdir=args.basedir,
                roi=None, target_empty=False, tile_filter=','.join(args.tile_filter),
                username=args.username,
                password=args.password, match_begin=True, print_results=False,
                download=True, overwrite=True, robust=True, max_retries=-1,
                multithread=True, nthreads=4, collection='006'
            )

            # anything downloaded? (or redo-smoothing: app_state.redo_smoothing)
            if len(downloaded) < 1 or getattr(args, 'download_only', False):
                break  # while True

        if getattr(args, 'collect_only', False) or (
                not getattr(args, 'smooth_only', False) and
                not getattr(args, 'export_only', False)):

            # check download completeness:
            if not curate_downloads(args.basedir, args.tile_filter, next_date, next_date):
                break

            # We're OK; now collect;
            modis_collect.callback(
                src_dir=args.basedir, targetdir=args.basedir,  # modape appends VIM to the targetdir
                compression='gzip', vam_code='VIM', interleave=True, parallel_tiles=1,
                cleanup=True, force=False, last_collected=None
            )

            if getattr(args, 'collect_only', False):
                break

        # smooth by N/n
        if getattr(args, 'smooth_only', False) or (
                not getattr(args, 'export_only', False)):

            modis_smooth.callback(
                src=os.path.join(args.basedir, 'VIM'),
                targetdir=os.path.join(args.basedir, 'VIM', 'SMOOTH'),
                svalue=None, srange=[], pvalue=None, tempint=10, tempint_start=None,
                nsmooth=args.nsmooth, nupdate=args.nupdate, soptimize=False,
                parallel_tiles=1, last_collected=None
            )

            if getattr(args, 'smooth_only', False):
                break

        # export dekads, from back to front (n = 6):
        nexports = 1
        export_octad = \
            ModisInterleavedOctad(get_last_date_in_raw_modis_tiles(os.path.join(args.basedir, 'VIM')))
        export_dekad = \
            Dekad(export_octad.getDateTimeEnd(), True)

        while Dekad(export_octad.prev().getDateTimeEnd(), True)\
                .Equals(export_dekad) and nexports <= args.nupdate:
            nexports = nexports + 1
            export_octad = export_octad.prev()

        first_date = get_first_date_in_raw_modis_tiles(os.path.join(args.basedir, 'VIM'))
        while (not export_dekad.startsBeforeDate(first_date)) and nexports <= args.nupdate:

            for region, roi in args.export.items():
                if getattr(args, 'region_only', region) != region:
                    continue
                exports = modis_window.callback(
                    src=os.path.join(args.basedir, 'VIM', 'SMOOTH'),
                    targetdir=os.path.join(args.basedir, 'VIM', 'SMOOTH', 'EXPORT'),
                    begin_date=export_dekad.getDateTimeMid(),
                    end_date=export_dekad.getDateTimeMid(),
                    # convert from LLX,LLY,URX,URY to ULX,ULY,LRX,LRY:
                    # roi=[roi[0], roi[3], roi[2], roi[1]],
                    roi=[roi[0], roi[1], roi[2], roi[3]],
                    region=region, sgrid=False, force_doy=False,
                    filter_product=None, filter_vampc=None, target_srs='EPSG:4326',
                    co=["COMPRESS=LZW", "PREDICTOR=2", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"],
                    clip_valid=True, round_int=2, gdal_kwarg={
                        'xRes': 0.01, 'yRes': 0.01,
                        'metadataOptions': ['CONSOLIDATION_STAGE={}'.format(nexports-1),
                                            'FINAL={}'.format('FALSE' if nexports < args.nupdate else 'TRUE')]
                    },
                    overwrite=True
                )

                for exp in exports:
                    md5 = generate_file_md5(exp)
                    with contextlib.suppress(FileNotFoundError):
                        os.remove(exp + '.md5')
                    with open(exp + '.md5', 'w') as f:
                        f.write(md5)

            nexports = nexports + 1
            export_octad = export_octad.prev()
            export_dekad = Dekad(export_octad.getDateTimeEnd(), True)
            while Dekad(export_octad.prev().getDateTimeEnd(), True)\
                    .Equals(export_dekad) and nexports <= args.nupdate:
                nexports = nexports + 1
                export_octad = export_octad.prev()

        if only_one_inc or getattr(args, 'export_only', False):
            break  # while True


def app_setup(config):
    global app_state
    with open(config) as f:
        app_state = json.load(f)
    app_state = Namespace(**app_state)

    flask_app = Flask(app_state.app_name)
    flask_app.config['JSONIFY_PRETTYPRINT_REGULAR'] = True
    flask_app.add_url_rule('/fetch', 'fetch', app_fetch)
    flask_app.add_url_rule('/suspend', 'suspend', app_suspend)
    flask_app.add_url_rule('/download/<filename>', 'download', app_download)
    flask_app.add_url_rule('/log/<filename>', 'log', app_log)
    flask_app.add_url_rule('/', 'index', app_index)
    app_state.fetcherThread = Thread()
    return flask_app


@click.group(invoke_without_command=True)
@click.option('--debug/--no-debug', default=False)
@click.option('--region')
@click.option('--config')
@click.pass_context
def cli(ctx, config, region, debug):
    ctx.ensure_object(dict)
    ctx.obj['CONFIG'] = config
    ctx.obj['REGION'] = region
    ctx.obj['DEBUG'] = debug
    if ctx.invoked_subcommand is None:
        ctx.invoke(serve)


@cli.command()
@click.pass_context
def serve(ctx) -> None:
    if ctx.obj['DEBUG']:
        with open(ctx.obj['CONFIG']) as f:
            args = json.load(f)
        args = Namespace(**args)
        if not exists_smooth_h5s(args.tile_filter, args.basedir):
            raise SystemExit(
                "Cannot run a full time step increment on an uninitialised archive! Run the init command first or run "
                "as a service. "
            )
        else:
            if ctx.obj['REGION']:
                args.region_only = ctx.obj['REGION']
            do_processing(args, only_one_inc=True)
    else:
        flask_app = app_setup(ctx.obj['CONFIG'])
        assert (ctx.obj['REGION'] is None), "Cannot serve only a specific region! Please run with the debug flag."
        flask_app.run(port=5001, threaded=False)  # Configure for single threaded request handling


@cli.command()
@click.pass_context
def export(ctx) -> None:
    with open(ctx.obj['CONFIG']) as f:
        args = json.load(f)
    args = Namespace(**args)

    # Check all tiles for a corresponding H5 archive in VIM/SMOOTH:
    assert exists_smooth_h5s(args.tile_filter, args.basedir)

    if ctx.obj['REGION']:
        args.region_only = ctx.obj['REGION']
    args.export_only = True
    do_processing(args)


@cli.command()
@click.pass_context
def smooth(ctx) -> None:
    with open(ctx.obj['CONFIG']) as f:
        args = json.load(f)
    args = Namespace(**args)

    assert (ctx.obj['REGION'] is None), "Cannot smooth for only a specific region!"
    args.smooth_only = True
    do_processing(args)


@cli.command()
@click.pass_context
def collect(ctx) -> None:
    with open(ctx.obj['CONFIG']) as f:
        args = json.load(f)
    args = Namespace(**args)
    assert (ctx.obj['REGION'] is None), "Cannot collect for only a specific region!"
    args.collect_only = True
    do_processing(args)


@cli.command()
@click.pass_context
def download(ctx) -> None:
    with open(ctx.obj['CONFIG']) as f:
        args = json.load(f)
    args = Namespace(**args)
    assert (ctx.obj['REGION'] is None), "Cannot download for only a specific region!"
    args.download_only = True
    do_processing(args)


@cli.command()
@click.pass_context
def reset(ctx) -> None:
    with open(ctx.obj['CONFIG']) as f:
        args = json.load(f)
    args = Namespace(**args)
    assert (ctx.obj['REGION'] is None), "Cannot reset for only a specific region!"
    while os.path.isdir(args.basedir):
        sure = input("Flushing the entire production environment. Are you sure? [y/n]: ").lower().strip()
        if sure == "y" or sure == "yes":
            shutil.rmtree(args.basedir)
            return
        elif sure == "n" or sure == "no":
            log.info("Aborted.")
            return
    os.makedirs(os.path.join(args.basedir, 'log'))
    log.info("Done.")


@cli.command()
@click.option('--download-only', is_flag=True, default=False)
@click.option('--smooth-only', is_flag=True, default=False)
@click.option('--export-only', is_flag=True, default=False)
@click.pass_context
def init(ctx, download_only, smooth_only, export_only) -> None:
    with open(ctx.obj['CONFIG']) as f:
        args = json.load(f)
    args = Namespace(**args)

    if ctx.obj['REGION'] is not None:
        assert (export_only and exists_smooth_h5s(args.tile_filter, args.basedir)),\
            "Can only do export for a specific region on a initialised archive!"
        args.region_only = ctx.obj['REGION']

    args.download_only = download_only
    args.smooth_only = smooth_only
    args.export_only = export_only
    do_init(args)


def do_init(args):
    if not getattr(args, 'smooth_only', False) and not getattr(args, 'export_only', False):
        # Download and Collect:
        # ---------------------

        begin_date = get_last_date_in_raw_modis_tiles(os.path.join(args.basedir, 'VIM'))
        if begin_date is None:
            begin_date = datetime.strptime(args.init_start_date, '%Y-%m-%d').date()
        else:
            begin_date = ModisInterleavedOctad(begin_date).next().getDateTimeStart().date()

        end_date = datetime.strptime(args.init_end_date, '%Y-%m-%d').date()
        if not getattr(args, 'download_only', False):
            # We can do incremental processing if we're not restricted to downloading only:
            end_date = min([end_date, begin_date.nextYear().prev().getDateTimeStart().date()])

        while begin_date < end_date:
            if getattr(args, 'suspended', False):
                return

            log.info('Downloading: {} - {}...'.format(begin_date, end_date))
            downloads = modis_download.callback(
                products=['M?D13A2'],
                begin_date=datetime.combine(begin_date, datetime.min.time()),
                end_date=datetime.combine(end_date, datetime.min.time()),
                targetdir=args.basedir,
                roi=None, target_empty=False, tile_filter=','.join(args.tile_filter),
                username=args.username,
                password=args.password, match_begin=True, print_results=False,
                download=True, overwrite=False, robust=True, max_retries=-1,
                multithread=True, nthreads=4, collection='006'
            )
            if len(downloads) == 0:
                break

            # Check: all downloads are found on disk?
            any_download_missing = False
            for filename in downloads:
                if not os.path.exists(os.path.join(args.basedir, filename)):
                    log.error('Download missing on disk: {}'.format(filename))
                    any_download_missing = True
            if any_download_missing:
                return

            if not getattr(args, 'download_only', False):
                # Check download: for all distinct dates: is there a download for EACH selected tile?
                if not curate_downloads(args.basedir, args.tile_filter, begin_date, end_date):
                    return

                # We're OK; now collect;
                modis_collect.callback(
                    src_dir=args.basedir, targetdir=args.basedir,  # modape appends VIM to the targetdir
                    compression='gzip', vam_code='VIM', interleave=True, parallel_tiles=1,
                    cleanup=True, force=False, last_collected=None
                )

                # move on:
                begin_date = get_last_date_in_raw_modis_tiles(
                    os.path.join(args.basedir, 'VIM'))
                begin_date = ModisInterleavedOctad(begin_date).next().getDateTimeStart().date()
                end_date = min([datetime.strptime(args.init_end_date, '%Y-%m-%d').date(),
                                begin_date.nextYear().prev().getDateTimeStart().date()])

        if getattr(args, 'download_only', False):
            return

    if not exists_smooth_h5s(args.tile_filter, args.basedir):
        # Smooth and interpolate the collected archive
        # --------------------------------------------

        # Check if the raw grid stacks (each tile) contain (and *only* contain) the configured date range
        # for initialisation: init_start_date -- init_end_date:
        begin_date = datetime.strptime(args.init_start_date, '%Y-%m-%d').date()
        end_date = datetime.strptime(args.init_end_date, '%Y-%m-%d').date()
        dates = []
        ts = ModisInterleavedOctad(begin_date)
        while ts.getDateTimeStart().date() < begin_date:
            ts = ts.next()
        while ts.getDateTimeStart().date() < end_date:
            dates.append(str(ts))
            ts = ts.next()
        for tile in args.tile_filter:
            assert has_collected_dates(os.path.join(args.basedir, 'VIM', "MXD13A2.{}.006.VIM.h5".format(tile)), dates)

        modis_smooth.callback(
            src=os.path.join(args.basedir, 'VIM'), targetdir=os.path.join(args.basedir, 'VIM', 'SMOOTH'),
            svalue=None, srange=[], pvalue=None, tempint=10, tempint_start=None,
            nsmooth=0, nupdate=0, soptimize=True, parallel_tiles=1, last_collected=None
        )

        if getattr(args, 'smooth_only', False):
            return

    # Export smoothened slices
    # ------------------------

    # Check all tiles for a corresponding H5 archive in VIM/SMOOTH:
    assert exists_smooth_h5s(args.tile_filter, args.basedir)

    first_date = max([
        get_first_date_in_raw_modis_tiles(os.path.join(args.basedir, 'VIM')),
        datetime.strptime(args.init_start_date, '%Y-%m-%d').date()
    ])
    last_date = datetime.combine(
        get_last_date_in_raw_modis_tiles(os.path.join(args.basedir, 'VIM')) + relativedelta(days=8),
        datetime.min.time()
    )

    export_slice = Dekad(first_date)
    if export_slice.startsBeforeDate(first_date):
        export_slice = export_slice.next()

    to_slice = export_slice
    cnt = 1
    while True:
        if cnt == 9 or to_slice.next().getDateTimeMid() > last_date:
            # Batch-wise export: 9 dekads at a time
            for region, roi in args.export.items():
                if hasattr(args, 'this_region_only') and args.this_region_only != region:
                    continue
                log.info('{} -- Exporting {} to {} ...'.format(region, str(export_slice), str(to_slice)))
                exports = modis_window.callback(
                    src=os.path.join(args.basedir, 'VIM', 'SMOOTH'),
                    targetdir=os.path.join(args.basedir, 'VIM', 'SMOOTH', 'EXPORT'),
                    begin_date=export_slice.getDateTimeMid(),
                    end_date=to_slice.getDateTimeMid(),
                    roi=[roi[0], roi[1], roi[2], roi[3]],
                    region=region, sgrid=False, force_doy=False,
                    filter_product=None, filter_vampc=None, target_srs='EPSG:4326',
                    co=["COMPRESS=LZW", "PREDICTOR=2", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"],
                    clip_valid=True, round_int=2, gdal_kwarg={
                        'xRes': 0.01, 'yRes': 0.01, 'metadataOptions': ['FINAL=TRUE']
                    },
                    overwrite=True
                )
                for exp in exports:
                    md5 = generate_file_md5(exp)
                    with contextlib.suppress(FileNotFoundError):
                        os.remove(exp + '.md5')
                    with open(exp + '.md5', 'w') as f:
                        f.write(md5)

        if to_slice.next().getDateTimeMid() > last_date:
            # every date represents (a) the *mid* of the one composite and the *start* of the other
            break

        if cnt == 9:
            export_slice = to_slice.next()
            to_slice = export_slice
            cnt = 1
        else:
            to_slice = to_slice.next()
            cnt = cnt + 1


if __name__ == '__main__':
    this_dir, _ = os.path.split(__file__)
    cli(default_map={
        'config': os.path.join(this_dir, 'production.json')
    })
