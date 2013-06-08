#!/usr/bin/env python
import os,sys,time,shutil,subprocess
from glob import glob

my_dir = os.path.abspath(os.path.dirname(__file__))

def cp_file(f_dir,f,d_dir,d):
    shutil.copy(os.path.join(f_dir,f),os.path.join(d_dir,d))

def main():
    def copy():
        for path in glob(os.path.join(my_dir,'*/*/*/plotpkl.py')):
            #print path
            d_dir = '/'.join(path.split('/')[:-1])
            cp_file(my_dir,'allplotpkl3.py',d_dir,'plotpkl_rev1.py')
    def run():
        count = 0
        for path in glob(os.path.join(my_dir,'*/*/*/plotpkl_rev1.py')):
            count +=1
            print path
            d_dir = '/'.join(path.split('/')[:-1])
            #print d_dir
            #print 'plotting %d times' % count
            os.chdir(d_dir)
            command = ['/usr/bin/python',path]
            pipe=subprocess.Popen(command,stdin=subprocess.PIPE, \
                       stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            stout,stderr = pipe.communicate()
            print stout
            print 'stderr >>> ',stderr
            if count >=60: break

    copy()
    run()
    # end main________________________________________________________

if __name__ == "__main__":
    main()
