#!/usr/bin/env python
"""
  arc_modis_ndvi.py: Flask Service app for collecting, processing and disseminating filtered NDVI.
           Production the time series leverages the WFP VAM MODAPE toolkit: https://github.com/WFP-VAM/modape

  Dependencies: arc-modape (1.0), Numpy, ...

  Author: Rob Marjot, (c) ARC 2020

"""
import contextlib
import hashlib
import json
import os
import sys
import glob
import click
import re

from datetime import datetime, date
from dateutil.relativedelta import relativedelta
from flask import Flask, jsonify, send_file
from threading import Thread, Timer

from modape_helper import get_first_date_in_raw_modis_tiles, get_last_date_in_raw_modis_tiles, curate_downloads
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


def app_index():
    global app_state
    if app_state.fetcherThread.is_alive() and not getattr(app_state, 'suspended', False):
        return "Fetcher is running (or suspended), try again later\n", 404
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
        return "Fetcher is running (or suspended), try again later\n", 404
    else:
        return send_file(os.path.join(app_state.basedir, 'VIM', 'SMOOTH', 'EXPORT', filename), as_attachment=True, mimetype=app_state.mimetype)


def app_fetch():
    global app_state
    if app_state.fetcherThread.is_alive() or getattr(app_state, 'suspended', False):
        return "[{}] Fetcher is already running (or suspended), try again later\n".format(
            datetime.now().strftime('%Y-%m-%d %H:%M:%S')), 404
    else:
        app_state.fetcherThread = Timer(5, app_do_processing, ())
        app_state.fetcherThread.start()
        return "[{}] Fetching and processing is scheduled to start in 5 seconds...\n".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


def app_suspend():
    app_state.suspended = True
    if app_state.fetcherThread.is_alive():
        return "Fetcher is busy suspending; please check back later to see the suspended state confirmed...\n", 404
    else:
        return "[{}] Fetcher suspended; restart to resume.\n".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


def app_do_processing(debug=False):
    global app_state
    try:
        # download and ingest:
        while True:
            last_date = get_last_date_in_raw_modis_tiles(os.path.join(app_state.basedir, 'VIM'))
            next_date = last_date + relativedelta(days=8)
            if last_date.year < next_date.year:
                # handle turning of the year:
                next_date = datetime(next_date.year, 1, 1).date()

            if getattr(app_state, 'download_only', False) or (
                    not getattr(app_state, 'collect_only', False) and
                    not getattr(app_state, 'smooth_only', False) and
                    not getattr(app_state, 'export_only', False)):

                if next_date > date.today():  # stop after today:
                    break

                print('Downloading: {}...'.format(next_date))
                downloaded = modis_download.callback(
                    products=['M?D13A2'],
                    begin_date=datetime.combine(next_date, datetime.min.time()),
                    end_date=datetime.combine(next_date, datetime.min.time()),
                    targetdir=app_state.basedir,
                    roi=None, target_empty=False, tile_filter=','.join(app_state.tile_filter),
                    username=app_state.username,
                    password=app_state.password, strict_dates=True, return_results=False,
                    download=True, overwrite=True, multithread=False, nthreads=1, collection='006'
                )

                # anything downloaded? (or redo-smoothing: app_state.redo_smoothing)
                if len(downloaded) < 1 or getattr(app_state, 'download_only', False):
                    break

            if getattr(app_state, 'collect_only', False) or (
                    not getattr(app_state, 'smooth_only', False) and
                    not getattr(app_state, 'export_only', False)):

                # check download completeness:
                if not curate_downloads(app_state.basedir, app_state.tile_filter, next_date, next_date):
                    break

                # We're OK; now collect;
                modis_collect.callback(
                    src_dir=app_state.basedir, targetdir=app_state.basedir,  # modape appends VIM to the targetdir
                    compression='gzip', vam_code='VIM', interleave=True, parallel_tiles=1,
                    cleanup=True, last_collected=None
                )

                if getattr(app_state, 'collect_only', False):
                    break

            # smooth by N/n
            if getattr(app_state, 'smooth_only', False) or (
                    not getattr(app_state, 'export_only', False)):

                modis_smooth.callback(
                    src=os.path.join(app_state.basedir, 'VIM'),
                    targetdir=os.path.join(app_state.basedir, 'VIM', 'SMOOTH'),
                    svalue=None, srange=[], pvalue=None, tempint=10, tempint_start=None,
                    nsmooth=app_state.nsmooth, nupdate=app_state.nupdate, soptimize=False,
                    parallel_tiles=1, last_collected=None
                )

                if getattr(app_state, 'smooth_only', False):
                    break

            # export dekads, from back to front (n = 6):
            nexports = 1
            export_octad = \
                ModisInterleavedOctad(get_last_date_in_raw_modis_tiles(os.path.join(app_state.basedir, 'VIM')))
            export_dekad = \
                Dekad(export_octad.getDateTimeEnd(), True)
            if debug:
                print('')
                print('Octad-end for last ingested date: {}'.format(str(export_octad.getDateTimeEnd())))
                print(' > Corresponding dekad: {}'.format(str(export_dekad)))

            while Dekad(export_octad.prev().getDateTimeEnd(), True)\
                    .Equals(export_dekad) and nexports <= app_state.nupdate:
                nexports = nexports + 1
                export_octad = export_octad.prev()

            first_date = get_first_date_in_raw_modis_tiles(os.path.join(app_state.basedir, 'VIM'))
            while (not export_dekad.startsBeforeDate(first_date)) and nexports <= app_state.nupdate:
                if debug:
                    print('>>Export: {} [Update: {}]'.format(str(export_dekad), str(nexports)))
                for region, roi in app_state.export.items():
                    if getattr(app_state, 'region_only', region) != region:
                        continue
                    exports = modis_window.callback(
                        src=os.path.join(app_state.basedir, 'VIM', 'SMOOTH'),
                        targetdir=os.path.join(app_state.basedir, 'VIM', 'SMOOTH', 'EXPORT'),
                        begin_date=export_dekad.getDateTimeMid(),
                        end_date=export_dekad.getDateTimeMid(),
                        # convert from LLX,LLY,URX,URY to ULX,ULY,LRX,LRY:
                        roi=[roi[0], roi[3], roi[2], roi[1]],
                        region=region, sgrid=False, force_doy=False,
                        filter_product=None, filter_vampc=None, target_srs='EPSG:4326',
                        co=["COMPRESS=LZW", "PREDICTOR=2", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"],
                        clip_valid=True, round_int=2, gdal_kwarg={
                            'xRes': 0.01, 'yRes': 0.01,
                            'metadataOptions': ['CONSOLIDATION_STAGE={}'.format(nexports-1),
                                                'FINAL={}'.format('FALSE' if nexports < app_state.nupdate else 'TRUE')]
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
                        .Equals(export_dekad) and nexports <= app_state.nupdate:
                    nexports = nexports + 1
                    export_octad = export_octad.prev()

            if debug or getattr(app_state, 'export_only', False):
                break

    finally:
        # app_state.fetcherThread = None
        pass


def setup(config):
    global app_state
    with open(config) as f:
        app_state = json.load(f)
        app_state = Namespace(**app_state)
    flask_app = Flask(app_state.app_name)
    flask_app.config['JSONIFY_PRETTYPRINT_REGULAR'] = True
    flask_app.add_url_rule('/fetch', 'fetch', app_fetch)
    flask_app.add_url_rule('/suspend', 'suspend', app_suspend)
    flask_app.add_url_rule('/download/<filename>', 'download', app_download)
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
    global app_state
    if ctx.obj['DEBUG']:
        with open(ctx.obj['CONFIG']) as f:
            app_state = json.load(f)
            app_state = Namespace(**app_state)
        if ctx.obj['REGION']:
            app_state.region_only = ctx.obj['REGION']
        app_do_processing(debug=True)
    else:
        flask_app = setup(ctx.obj['CONFIG'])
        assert (ctx.obj['REGION'] is None), "Cannot serve only a specific region! Please run with the debug flag."
        flask_app.run(port=5001, threaded=False)  # Configure for single threaded request handling


@cli.command()
@click.pass_context
def export(ctx) -> None:
    global app_state
    with open(ctx.obj['CONFIG']) as f:
        app_state = json.load(f)
    app_state = Namespace(**app_state)
    if ctx.obj['REGION']:
        app_state.region_only = ctx.obj['REGION']
    app_state.export_only = True
    app_do_processing(debug=ctx.obj['DEBUG'])


@cli.command()
@click.pass_context
def smooth(ctx) -> None:
    global app_state
    with open(ctx.obj['CONFIG']) as f:
        app_state = json.load(f)
    app_state = Namespace(**app_state)
    assert (ctx.obj['REGION'] is None), "Cannot smooth for only a specific region!"
    app_state.smooth_only = True
    app_do_processing(debug=ctx.obj['DEBUG'])


@cli.command()
@click.pass_context
def collect(ctx) -> None:
    global app_state
    with open(ctx.obj['CONFIG']) as f:
        app_state = json.load(f)
    app_state = Namespace(**app_state)
    assert (ctx.obj['REGION'] is None), "Cannot collect for only a specific region!"
    app_state.collect_only = True
    app_do_processing(debug=ctx.obj['DEBUG'])


@cli.command()
@click.pass_context
def download(ctx) -> None:
    global app_state
    with open(ctx.obj['CONFIG']) as f:
        app_state = json.load(f)
    app_state = Namespace(**app_state)
    assert (ctx.obj['REGION'] is None), "Cannot download for only a specific region!"
    app_state.download_only = True
    app_do_processing(debug=ctx.obj['DEBUG'])


@cli.command()
@click.option('--download-only', is_flag=True, default=False)
@click.option('--smooth-only', is_flag=True, default=False)
@click.option('--export-only', is_flag=True, default=False)
@click.pass_context
def init(ctx, download_only, smooth_only, export_only) -> None:
    with open(ctx.obj['CONFIG']) as f:
        args = json.load(f)
        args = Namespace(**args)
    assert (ctx.obj['REGION'] is None), "Cannot initialize for only a specific region!"

    urls = []
    if not smooth_only and not export_only:
        end_date = None
        # download and ingest:
        begin_date = get_last_date_in_raw_modis_tiles(os.path.join(args.basedir, 'VIM'))
        if begin_date is None:
            begin_date = datetime.strptime(args.init_start_date, '%Y-%m-%d').date()
        else:
            archive_year = begin_date.year
            begin_date = begin_date + relativedelta(days=8)
            if archive_year < begin_date.year:
                # handle turning of the year:
                begin_date = datetime(begin_date.year, 1, 1).date()

        end_date = datetime.strptime(args.init_end_date, '%Y-%m-%d').date()
        if not download_only:
            end_date = min([end_date, begin_date + relativedelta(years=1) - relativedelta(days=1)])

        while begin_date < end_date:
            print('Downloading: {} - {}...'.format(begin_date, end_date))
            urls = modis_download.callback(
                products=['M?D13A2'],
                begin_date=datetime.combine(begin_date, datetime.min.time()),
                end_date=datetime.combine(end_date, datetime.min.time()),
                targetdir=args.basedir,
                roi=None, target_empty=False, tile_filter=','.join(args.tile_filter),
                username=args.username,
                password=args.password, strict_dates=True, return_results=False,
                download=True, overwrite=False, multithread=True, nthreads=4, collection='006'
            )
            if len(urls) == 0:
                break

            # Check download: all urls are found on disk?
            for url in urls:
                fname = url['file_id'][url['file_id'].rfind('/') + 1:]
                if not os.path.exists(os.path.join(args.basedir, fname)):
                    raise SystemExit('Download missing on disk: {}'.format(fname))

            if not download_only:
                # Check download: for all distinct dates: is there a download for EACH selected tile?
                if not curate_downloads(args.basedir, args.tile_filter, begin_date, end_date):
                    exit(1)

                # We're OK; now collect;
                modis_collect.callback(
                    src_dir=args.basedir, targetdir=args.basedir,  # modape appends VIM to the targetdir
                    compression='gzip', vam_code='VIM', interleave=True, parallel_tiles=1,
                    cleanup=True, last_collected=None
                )

                # move on:
                begin_date = get_last_date_in_raw_modis_tiles(
                    os.path.join(args.basedir, 'VIM')) + relativedelta(days=8)
                end_date = min([datetime.strptime(args.init_end_date, '%Y-%m-%d').date(),
                                begin_date + relativedelta(years=1) - relativedelta(days=1)])

        if download_only:
            sys.exit(int(not curate_downloads(args.basedir, args.tile_filter, begin_date, end_date)))

    if len(urls) > 0 or smooth_only:
        # smooth downloaded archive: setting the 'init_only' to True, this can be done only once per product tile:
        modis_smooth.callback(
            src=os.path.join(args.basedir, 'VIM'), targetdir=os.path.join(args.basedir, 'VIM', 'SMOOTH'),
            svalue=None, srange=[], pvalue=None, tempint=10, tempint_start=None,
            nsmooth=0, nupdate=0, soptimize=True, parallel_tiles=1, last_collected=None
        )

        if smooth_only:
            exit(0)

    # Export dekads:
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
            for region, roi in args.export.items():
                if (ctx.obj['REGION'] is not None) and ctx.obj['REGION'] != region:
                    continue
                print('\n{} -- Exporting {} to {} ...'.format(region, str(export_slice), str(to_slice)))
                exports = modis_window.callback(
                    src=os.path.join(args.basedir, 'VIM', 'SMOOTH'),
                    targetdir=os.path.join(args.basedir, 'VIM', 'SMOOTH', 'EXPORT'),
                    begin_date=export_slice.getDateTimeMid(),
                    end_date=to_slice.getDateTimeMid(),
                    roi=[roi[0], roi[3], roi[2], roi[1]],
                    # convert from LLX,LLY,URX,URY to ULX,ULY,LRX,LRY
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
