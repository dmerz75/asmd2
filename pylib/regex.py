import os,sys,re,shutil

def reg_ex(s_dir,script,subout,subin):
    ''' perform regular expressions on a text file
        i.e. (script_directory,script,'xx_sub_xx','regular expression')
    '''
    full_path_script = os.path.join(s_dir,script)
    o = open(full_path_script,'r+')
    text=o.read()
    text=re.sub(subout,subin,text)
    o.close()
    o = open(full_path_script,'w+')
    o.write(text)
    o.close()
