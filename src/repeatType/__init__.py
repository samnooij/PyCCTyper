#!/usr/bin/env python

import argparse

from cctyper.xgb import XGB
from cctyper.repeat import RepeatTyper


def parse_arguments():
    """
    Parse arguments from the command-line
    """
    ap = argparse.ArgumentParser()

    # Required
    ap.add_argument(
        "input",
        help="Input. A simple text file with one repeat sequence per line"
    )

    # Optional
    ap.add_argument(
        "--kmer",
        help="kmer size. Has to match training kmer size! [%(default)s].",
        default=4,
        type=int,
    )
    ap.add_argument("--db", help="Path to database.", default="", type=str)

    return ap


# Workflow starts here
def main():
    ap = parse_arguments()
    master = RepeatTyper(ap.parse_args())
    xgb = XGB(master)
    xgb.predict_repeats()
    xgb.print_xgb()

    return 0


if __name__ == "__main__":
    exit(main())
