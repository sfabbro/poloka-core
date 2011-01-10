# -*- mode: pyton; -*- 

import os
import os.path as op
import sys



def get_out_name():
    """
    Generate an (architecture dependant) 
    build directory name.
    """
    ret = os.uname()
    out_name = 'build.'+ ret[0]+'-'+ret[-1]
    return op.join('build', out_name)

    

def options(opt):
    """
    Default options for the frogs' packages.
    """
    opt.load('compiler_cc')
    opt.load('compiler_cxx')
    opt.load('compiler_fc')
    
    opt.add_option('--with-cfitsio', 
                   action='store', 
                   help='Path to the cfitsio root dir')
    opt.add_option('--with-sex', 
                   action='store', 
                   help='Path to the cfitsio root dir')    
    opt.add_option('--no-global-lapack', 
                   action='store_false',
                   default=True, 
                   dest='global_lapack',
                   help='use/no do use system lapack library')    
    opt.add_option("--no-cernlib", 
                   action='store_false', 
                   default=True, 
                   dest='cernlib', 
                   help='do not try to link with cernlib')

    
