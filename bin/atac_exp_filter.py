#!/usr/bin/env python3

import os
import sys

exp_path = sys.argv[1]
out_path = sys.argv[2]
chr = sys.argv[3]


fi = open(exp_path,'rt')
fo = open(out_path,'w')
count = 0
for l in fi:
    count = count + 1
    if(count < 2):
        fo.writelines(l)
    else:
        row = l.split()
        if(row[1] == chr):
            fo.writelines(l)
fo.close()
fi.close() 