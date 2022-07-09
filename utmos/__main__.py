#!/usr/bin/env python3
"""
Utmos main entrypoint
"""
import sys
import argparse

from utmos import __version__
from utmos.convert import cvt_main
from utmos.calculate import calc_main
# from utmos.plot import plot_main

def version(args):
    """Print the version"""
    print(f"Utmos v{__version__}")

TOOLS = {'convert': cvt_main,
         'select': calc_main,
         #'plot': plot_main
        }

USAGE = f"""\
Utmos v{__version__} - Reimplementation of SVCollector

    CMDs:
        convert  Extract genotypes from VCFs
        select   Select samples
        plot     Create qc plots (in-dev)
"""

def main():
    """
    Main entrypoint for utmos
    """
    parser = argparse.ArgumentParser(prog="utmos", description=USAGE,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="CMD", choices=TOOLS.keys(), type=str, default=None,
                        help="Command to execute")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()
    args = parser.parse_args()

    TOOLS[args.cmd](args.options)


if __name__ == '__main__':
    main()
