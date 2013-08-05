#!/usr/bin/env python
import os,sys,time
my_dir=os.path.abspath(os.path.dirname(__file__))
import cPickle as pickle

sys.path.append(my_dir)
from asmd.asmdmethod import *

# get pylib
from pylib.cp import *
from pylib.regex import *

print "HELLO from the builder"

from glob import glob
asmdmethods = glob(os.path.join(my_dir,'*/*.pkl'))

def asmd_expansion(dir_loc,pkl):
    os.chdir(dir_loc)
    f = pickle.load(open(pkl,'r'))
    print dir(f)
    attrs = vars(f)
    print ' \n'.join("%s: %s" % item for item in attrs.items())
    f.populate_work_dir()

for i in asmdmethods:
    print ('/').join(i.split('/')[0:-1])
    print i
    pkl_dir = ('/').join(i.split('/')[0:-1])
    asmd_expansion(pkl_dir,i)
    

