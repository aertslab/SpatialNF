#!/usr/bin/env python3

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(description='Format sample coord files.')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input coord csv file.'
)

parser.add_argument(
    "-x", "--x_col_name",
    type=str,
    action="store",
    dest="x_col_name",
    default="x",
    help="Column name for x coord."
)

parser.add_argument(
    "-y", "--y_col_name",
    type=str,
    action="store",
    dest="y_col_name",
    default="y",
    help="Column name for y coord."
)

parser.add_argument(
    "-z", "--z_col_name",
    type=str,
    action="store",
    dest="z_col_name",
    default="z",
    help="Column name for z coord."
)

parser.add_argument(
    "-g", "--gene_col_name",
    type=str,
    action="store",
    dest="gene_col_name",
    default="gene",
    help="Column name for gene."
)

parser.add_argument(
    "-s", "--sample_id",
    type=str,
    action="store",
    dest="sample_id",
    default=None,
    help="Sample id."
)

parser.add_argument(
    "-d", "--delimiter",
    type=str,
    action="store",
    dest="delimiter",
    default=",",
    help="Column separation char."
)

parser.add_argument(
    "--header",
    type=str,
    action="store",
    dest="header",
    default="infer",
    help="Header line."
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output formatted csv file.'
)

args = parser.parse_args()

FILE_PATH_IN = args.input
FILE_PATH_OUT = args.output.name

FILE_PATH_IN = FILE_PATH_IN.name
df = pd.read_csv(FILE_PATH_IN, sep=args.delimiter, header=args.header)

df = df.loc[:, [args.x_col_name, args.y_col_name, args.z_col_name, args.gene_col_name]].copy()
df.columns = ['x', 'y', 'z', 'gene']
df[['sample_id']] = args.sample_id
       
df.to_csv(FILE_PATH_OUT, sep=',', index=False)
