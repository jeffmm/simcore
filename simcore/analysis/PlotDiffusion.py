#!/usr/local/bin/python

f = open("../sc-diffusion.log", "r")
line = f.readline()
line = line.split()
nFilaments = len(line) - 1
time = []
dat = []
for i in range(nFilaments):
    dat.append([])
for line in f:
    line = line.split()
    time.append(float(line[0]))
    for i in range(nFilaments):
        dat[i].append(float(line[i + 1]))
f.close()
