#!/usr/local/bin/python

import matplotlib
import numpy as np
from scipy import stats
# matplotlib.use("macosx")
import matplotlib.pyplot as plt

f = open("../sc-diffusion.log",'r')
line = f.readline()
line=line.split()
nFilaments=len(line)-1
time=[]
dat=[]
for i in xrange(nFilaments):
    dat.append([])
for line in f:
    line=line.split()
    time.append(float(line[0]))
    for i in xrange(nFilaments):
        dat[i].append(float(line[i+1]))
f.close()

