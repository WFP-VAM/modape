import sys
from itertools import chain
from pathlib import Path
from typing import Generator

import click
import h5py


@click.group()
def cli():
    pass


@cli.command("range")
@click.argument("pathglob", type=click.STRING)
@click.argument("min", type=click.INT)
@click.argument("max", type=click.INT)
def edit_range(pathglob: str, min: int, max: int):
    assert min <= max

    path = Path("./")
    glob = pathglob
    try:
        split_at = pathglob.rindex("/")
        path = Path(pathglob[:split_at])
        glob = pathglob[split_at + 1 :]
    except ValueError:
        pass

    for filepath in path.glob(glob):
        click.echo(f"Updating: {filepath.as_posix()}")
        with h5py.File(filepath, "r+", libver="latest") as h5f:
            ds_data = h5f.get("data")
            ds_data.attrs.update(dict(valid_range=(min, max)))


@cli.command("dates")
@click.argument("pathglob", type=click.STRING)
@click.argument("start", type=click.INT)
@click.argument("end", type=click.INT)
def check_dates(pathglob: str, start: int, end: int):

    class Oct46Generator:
        def __init__(self, year: int, within: tuple[int, int]):
            self.__year = year
            self.__within = within

        def __iter__(self) -> Generator[int, None, None]:
            for doy in range(1, 367, 8):
                step = (self.__year * 1000) + doy
                if step < self.__within[0] or step > self.__within[1]:
                    continue
                yield step

        @staticmethod
        def for_start_end(start, end):
            """
            Parameters:
              - start: integer YYYYDOY
              - end: integer YYYYDOY
            """
            return chain(
                *[
                    Oct46Generator(year, (start, end))
                    for year in range(start // 1000, (end // 1000) + 1)
                ]
            )

    path = Path("./")
    glob = pathglob
    try:
        split_at = pathglob.rindex("/")
        path = Path(pathglob[:split_at])
        glob = pathglob[split_at + 1 :]
    except ValueError:
        pass

    checked = 0
    steps = 0
    for filepath in path.glob(glob):
        with h5py.File(filepath, "r", libver="latest") as h5f:
            dates = h5f.get("dates")
            date_values = set([int(x.decode()) for x in dates[...]])
            steps = 0
            for year_doy in Oct46Generator.for_start_end(start, end):
                if year_doy not in date_values:
                    click.echo(f"Missing date {year_doy} in file {filepath.name}!")
                steps += 1
        checked += 1

    click.echo(f"Done. Checked: {checked} files for {steps} dates (each).")


def cli_wrap():
    """Wrapper for cli"""

    if len(sys.argv) == 1:
        cli.main(["--help"])
    else:
        cli()  # pylint: disable=E1120


if __name__ == "__main__":
    cli_wrap()
