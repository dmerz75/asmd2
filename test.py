#!/usr/bin/env python
import os,sys,time

my_dir = os.path.abspath(os.path.dirname(__file__))

def print_stuff(*args):
    for i,obj in enumerate(*args):
        print i,obj

print_stuff([76,'string',824.1])

def print_ab(*args):
    keys = ['a','b','c','d']
    # for i,obj in enumerate(*args):
    #     print i,obj
    dct = dict(zip(keys,*args))
    print dct
    return dct

alist = [0.1981,0.2182,0.3154]
blist = ['vac','imp']
clist = [100,10,1]
dlist = ['black']

dct_list = [print_ab([a,b,c,d]) for a in alist for b in blist for c in clist for \
 d in dlist]

# print len(dct_list)
# print dct_list

#print time.strftime("%H:%M:%S",60)

print my_dir
print os.getcwd()
base = os.path.basename(my_dir)
print base
print os.path.lexists(my_dir)
