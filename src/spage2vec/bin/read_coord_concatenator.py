#!/usr/bin/env python3

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(description='Concatenate sample coord files.')

parser.add_argument(
    "input",
    nargs='+',
    type=argparse.FileType('r'),
    help='Input h5ad files.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output concatenated csv file.'
)

args = parser.parse_args()

FILE_PATH_OUT = args.output.name

files = []
for FILE_PATH_IN in args.input:
    print(FILE_PATH_IN)
    try:
        FILE_PATH_IN = FILE_PATH_IN.name
        df = pd.read_csv(FILE_PATH_IN, sep=',')
        files.append(df)
    except IOError:
        raise Exception("VSN ERROR: Wrong input format. Expects .csv files, got .{}".format(FILE_PATH_IN))
df = pd.concat(files)
df.reset_index(inplace=True, drop=True)

df.to_csv(FILE_PATH_OUT, sep=',', index=False)
