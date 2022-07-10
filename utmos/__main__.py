#!/usr/bin/env python3
"""
Utmos main entrypoint
"""
import sys
import argparse

from utmos import __version__
from utmos.convert import cvt_main
from utmos.select import select_main

def version(args):
    """Print the version"""
    print(f"Utmos v{__version__}")

TOOLS = {'convert': cvt_main,
         'select': select_main,
         'version': version
        }

USAGE = f"""\
Utmos v{__version__} - Maximum-coverage algorithm to select samples for validation and resequencing

    CMDs:
        convert  Extract genotypes from VCFs
        select   Select samples
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
