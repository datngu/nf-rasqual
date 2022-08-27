#!/usr/bin/env python

import pandas as pd
import os
import sys

meta_path = sys.argv[1]
meta_path = "brain_aqua-faang.csv"
meta = pd.read_csv(meta_path, header = 0)
#os.getcwd()

for i in meta.index:
    new_bam = meta['genotype_id'][i] + ".bam"
    old_bam = meta['attac_bam_id'][i]
    cmd = "cp " + old_bam + " " + new_bam
    print(cmd)
    os.system(cmd)
    new_bam2 = meta['genotype_id'][i] + ".bam.bai"
    old_bam2 = meta['attac_bam_id'][i] + ".bai"
    cmd2 = "cp " + old_bam2 + " " + new_bam2
    print(cmd2)
    os.system(cmd2)

