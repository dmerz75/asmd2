#!/usr/bin/env python
import sys,os,time
import numpy as np
from glob import *

my_dir = os.path.abspath(os.path.dirname(__file__))

for path in glob(os.path.join(my_dir,'[0-9][0-9].txt')):
    x = np.loadtxt(path)
    print x.shape
    y = x.reshape(10,2)
    ddir = '/'.join(path.split('/')[0:-1])
    num = path.split('/')[-1].split('.')[-2]
    newpath = os.path.join(ddir,num)+'/'+'%s.txt' % num
    print newpath
    print y
    #np.savetxt(newpath,y,fmt='%5.f')
