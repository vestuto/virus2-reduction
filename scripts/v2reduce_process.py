#!/usr/bin/env python

import argparse as ap
import sys

from v2reduce.reduce import process


def get_args():

    parser = ap.ArgumentParser(add_help=True)

    parser.add_argument("infolder", help='''Input folder''', type=str)

    parser.add_argument('outfolder', type=str,
                        help='''name of the output file''')

    parser.add_argument('date', type=str,
                        help='''string name for date, ex: 20250613''')

    parser.add_argument('-n', '--name', type=str,
                        help='''Name of the science target''', default=None)

    parser.add_argument("-ra", "--reduce_all",
                        help='''Reduce all files in folder''',
                        action="count", default=0)

    parser.add_argument("-bl", "--bias_label",
                        help='''The objet name for bias files''',
                        type=str, default='bias')

    parser.add_argument("-al", "--arc_label",
                        help='''The objet name for arc files''',
                        type=str, default='arc')

    parser.add_argument("-dl", "--dark_label",
                        help='''The objet name for arc files''',
                        type=str, default='dark')

    parser.add_argument("-fl", "--flat_label",
                        help='''The objet name for dome flat files''',
                        type=str, default='flat')

    parser.add_argument("-tfl", "--twilight_flat_label",
                        help='''The objet name for twilight flat files''',
                        type=str, default='twi')
    argv = None
    args = parser.parse_args(args=argv)

    return args


def main():

    args = get_args()

    process(args.infolder, args.outfolder,
            args.date, args.target_name, args.reduce_all,
            args.bias_label, args.arc_label, args.dark_label, args.flat_label, args.twilight_flat_label
            )

    return None


if __name__ == '__main__':
    sys.exit(main())
