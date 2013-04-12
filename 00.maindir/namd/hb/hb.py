#!/usr/bin/env python
import MDAnalysis
import MDAnalysis.analysis.hbonds
import MDAnalysis.analysis.distances
from sys import argv
import numpy as np
import os,sys,pickle

#______________universe________________________________________________________
u = MDAnalysis.Universe('../../../../00.struc/xxenvironxx/00.psf','daOut.dcd',\
                        permissive=True)

def analyze_bond(univ,seg1,seg2):
    try:
        name1=seg1.replace(' ','')
        name2=seg2.replace(' ','')
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(univ,seg1,seg2, \
                                                 distance=4.0, angle=140.0)
        results = h.run()
        pickle.dump(h.timeseries,open('%s-hb_%s_%s.pkl.%s' % (sys.argv[1], \
                                              name1,name2,sys.argv[2]),'w'))
    except:
        pass

#__analyze__bonds______________________________________________________________
analyze_bond(u,'protein','protein')
analyze_bond(u,'protein','segid WT1')
