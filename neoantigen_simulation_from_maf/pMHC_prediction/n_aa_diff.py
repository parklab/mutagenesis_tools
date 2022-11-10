#!/usr/bin/env python


import os,sys,re

col1 = int(sys.argv[2])
col2 = int(sys.argv[3])
infile = sys.argv[1]

with open(infile,'r') as f:
    for i in f:
        iline = i.strip('\n').split('\t')
        str1 = iline[col1 - 1]
        str2 = iline[col2 - 1]
        ndiff = abs(len(str1) - len(str2))
        for j,k in zip(str1,str2):
            if j != k:
                ndiff += 1
        iline += [str(ndiff)]
        print('\t'.join(iline))
