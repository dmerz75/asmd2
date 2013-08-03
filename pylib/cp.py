import os,sys,re,shutil

def cp_file(f_dir,f,d_dir,d):
    ''' copy a file
    '''
    shutil.copy(os.path.join(f_dir,f),os.path.join(d_dir,d))

def cp_tree(f_dir,f,d_dir,d):
    ''' copy a directory
        i.e. (file_dir,filename,destination_dir,as filename)
    '''
    shutil.copytree(os.path.join(f_dir,f),os.path.join(d_dir,d))
