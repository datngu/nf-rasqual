#!/usr/bin/env python

import sys
import argparse
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description = "Merge lsgkm weights of many models ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--weight", nargs='+', help = "list of lsgkm weights", required = True)
parser.add_argument("--out", default = "merged_weights.txt", help = "output file name")

args = parser.parse_args()

path_list = args.weight
out = args.out

#path_list = ["10_nr10mer_scores.txt", "1_nr10mer_scores.txt", "9_nr10mer_scores.txt", "8_nr10mer_scores.txt"]

df = pd.DataFrame()
for i, path in enumerate(path_list):
    print(path)
    tem = pd.read_csv(path, sep = "\t", header=None)
    df[i] = tem[1]


res = pd.DataFrame()
res[0] = tem[0]
res[1] = df.mean(axis=1)

res.to_csv(out, sep = "\t", header = False, index = False)


