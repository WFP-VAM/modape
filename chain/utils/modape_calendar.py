#!/usr/bin/env python
"""
  modape_calendar.py: Utility script to produce the Release Calendar and MODAPE CLI parameters 

  Dependencies: arc-modape (1.0), Numpy, ...

  Author: Rob Marjot, (c) ARC 2021

"""

import datetime
from dateutil.relativedelta import relativedelta
from modape_helper.timeslicing import Dekad, ModisInterleavedOctad

def generate_params(begin_date, end_date, cs, filename, html):
    cs += 1
    curModisStep = ModisInterleavedOctad(begin_date)

    with open(filename, 'w') as w:

        if html:
            w.write('<table>\n')

        b_count = 0
        while curModisStep.getDateTimeStart() < end_date:
            b_count += 1
            if html:
                w.write('  <tr>\n')
            export_octad = ModisInterleavedOctad(curModisStep.getDateTimeStart())
            export_dekad = Dekad(export_octad.getDateTimeEnd(), True)
            collected_octad = export_octad.prev()

            if html:
                w.write('    <td>&lt;YYYY&gt;{}</td>\n'.format(export_octad.getDateTimeStart().strftime('%j')))

                w.write('    <td>download -b/e: {}<br>window -b/e:<ul>'.format(
                    curModisStep.getDateTimeStart().strftime('%Y-%m-%d'))
                )
            else:
                w.write('# B{0}:\npython modis_download.py --download --multithread \\\n'
                        '  --username=$CMR_USERNAME --password=$CMR_PASSWORD \\\n'
                        '  --robust --target-empty --match-begin --tile-filter $TILES \\\n'
                        '  -b {1} -e {1} M?D13A2\n'.format(
                    b_count,
                    curModisStep.getDateTimeStart().strftime('%Y-%m-%d')
                ))
                w.write('python modis_collect.py --interleave --cleanup --last-collected {} .\n'.format(
                    collected_octad.getDateTimeStart().strftime('%Y%j')
                ))
                collected_octad = collected_octad.next()
                w.write('python modis_smooth.py --nsmooth 64 --nupdate 6 --tempint 10 \\\n'
                        '  --last-collected {} -d ./VIM/SMOOTH ./VIM\n'.format(
                    collected_octad.getDateTimeStart().strftime('%Y%j')
                ))
            nexports = 1
            while Dekad(export_octad.prev().getDateTimeEnd(), True).Equals(export_dekad) and nexports <= cs:
                nexports = nexports + 1
                export_octad = export_octad.prev()

            while nexports <= cs:
                exp_date = export_dekad.getDateTimeMid()
                if html:
                    w.write('<li>&lt;YYYY{}&gt-{}</li>'.format(
                        '' if exp_date.year == curModisStep.Year else '-1', exp_date.strftime('%m-%d')
                    ))
                else:
                    w.write('python modis_window.py -b {0} -e {0} --clip-valid --round-int 2 \\\n'
                            '  --roi -26.0,-35.0,58.0,38.0 --gdal-kwarg xRes=0.01 --gdal-kwarg yRes=0.01 \\\n'
                            '  --overwrite --last-smoothed {1} -d ./VIM/SMOOTH/EXPORT ./VIM/SMOOTH\n'.format(
                        exp_date.strftime('%Y-%m-%d'),
                        collected_octad.getDateTimeStart().strftime('%Y%j')
                    ))
                nexports = nexports + 1
                export_octad = export_octad.prev()
                export_dekad = Dekad(export_octad.getDateTimeEnd(), True)
                while Dekad(export_octad.prev().getDateTimeEnd(), True).Equals(export_dekad) and nexports <= cs:
                    nexports = nexports + 1
                    export_octad = export_octad.prev()

            curModisStep = curModisStep.next()
            if html:
                w.write('</ul></td>\n  </tr>\n')
            else:
                w.write('\n')

        if html:
            w.write('  </table>\n')

def build_calendar(begin_date, end_date, filename):
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    with open(filename, 'w') as w:

        iterDekad = Dekad(begin_date)
        w.write('date;doy')
        while iterDekad.getDateTimeStart() < end_date:
            w.write(';{}-d{};'.format(months[iterDekad.Month - 1], iterDekad.SliceInMonth))
            iterDekad = iterDekad.next()
        w.write('\n')

        curModisStep = ModisInterleavedOctad(begin_date)
        while curModisStep.getDateTimeEnd() + relativedelta(days=2) > begin_date:
            curModisStep = curModisStep.prev()

        expectedRelease = curModisStep.next().getDateTimeEnd() + relativedelta(days=2)
        while expectedRelease < end_date:

            stages = {}
            nexports = 1
            export_octad = ModisInterleavedOctad(curModisStep.getDateTimeStart())
            export_dekad = Dekad(export_octad.getDateTimeEnd(), True)
            while Dekad(export_octad.prev().getDateTimeEnd(), True).Equals(export_dekad) and nexports <= 6:
                nexports = nexports + 1
                export_octad = export_octad.prev()

            while nexports <= 6:
                latency = int((expectedRelease - export_dekad.getDateTimeEnd()) / datetime.timedelta(days=1))
                stages['{}-d{}'.format(
                    months[export_dekad.Month - 1], export_dekad.SliceInMonth)] = '{};{}'.format(nexports - 1,
                                                                                                   latency)

                nexports = nexports + 1
                export_octad = export_octad.prev()
                export_dekad = Dekad(export_octad.getDateTimeEnd(), True)
                while Dekad(export_octad.prev().getDateTimeEnd(), True).Equals(export_dekad) and nexports <= 6:
                    nexports = nexports + 1
                    export_octad = export_octad.prev()

            iterDekad = Dekad(begin_date)
            w.write('"{}-{}";"{:03d}"'.format(
                months[expectedRelease.month - 1],
                expectedRelease.day,
                int(curModisStep.getDateTimeStart().strftime('%j'))
            ))
            while iterDekad.getDateTimeStart() < end_date:
                w.write(';{}'.format(stages.get('{}-d{}'.format(
                    months[iterDekad.Month - 1], iterDekad.SliceInMonth), ';')))
                iterDekad = iterDekad.next()
            w.write('\n')

            curModisStep = curModisStep.next()
            expectedRelease = curModisStep.next().getDateTimeEnd() + relativedelta(days=2)

if __name__ == '__main__':
    # generate_params(datetime.datetime(2019, 1, 1), datetime.datetime(2019, 12, 31), 2, 'params.html', True)
    # generate_params(datetime.datetime(2020, 1, 1), datetime.datetime(2020, 12, 31), 2, 'paramsLeap.html', True)
    # generate_params(datetime.datetime(2020, 1, 1), datetime.datetime(2020, 12, 31), 2, 'forward2020.sh', False)
    # generate_params(datetime.datetime(2021, 1, 1), datetime.datetime(2021, 12, 31), 2, 'forward2021.sh', False)
    # build_calendar(datetime.datetime(2019, 1, 1), datetime.datetime(2019, 12, 31), 'calendar.txt')
    # build_calendar(datetime.datetime(2020, 1, 1), datetime.datetime(2020, 12, 31), 'calendarLeap.txt')
