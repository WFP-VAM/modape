#!/usr/bin/env python
""" Print modape version info"""
from . import __version__

def main():
    """main method"""
    msg = """
    +-+ +-+ +-+ +-+ +-+ +-+
    |m| |o| |d| |a| |p| |e|
    +-+ +-+ +-+ +-+ +-+ +-+
    Ahoy there! You are using modape version {}
    Godspeed.
    """.format(__version__)

    print(msg)
