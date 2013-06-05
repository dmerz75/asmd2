#!/usr/bin/env python
import os,sys,time
import numpy as np
from glob import *

my_dir = os.path.abspath(os.path.dirname(__file__))

content_list = []

dir_list = [str(f).zfill(2) for f in range(1,12)]
d = dir_list[0]

winners = []
def stage_comparison(d):
    dct_text = {}
    for path in glob('%s/*/daOut.coor.*' % d):
        print path
        sd = path.split('/')[-1].split('.')[-1]
        print sd
        f = open(path)
        text = f.read()
        f.close()
        dct_text[sd]=text
        content_list.append(text)

    real = []
    nextdir = str(int(d)+1).zfill(2)
    for path in glob('%s/*.coor' % nextdir):
        print path
        f = open(path)
        text = f.read()
        f.close()
        real.append(text)

    keys = dct_text.keys()

    for k in keys:
        try:
            print dct_text[k].split('\n')[-3]
            print k
            print real[0].split('\n')[-3]
        except:
            print "possible error in file %s" % k
        if dct_text[k].split('\n')[-3] == real[0].split('\n')[-3]:
            print 'Winner is %s' % k
            winners.append(k)

[stage_comparison(d) for d in dir_list]
print winners
