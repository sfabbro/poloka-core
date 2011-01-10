# -*- mode: python; -*- 


import os
import os.path as op
import sys
import commands
import shlex

import Options
import Configure


APPNAME  = 'poloka'
VERSION  = '0.1.0'
top   = '.'
out   = 'build'
description = "This is poloka."
requirements = [
    ('sex', '2.4.4', True), 
    ('cfitsio', '3.006', True) ]


def options(opt):
    
    opt.load('frogs')
    opt.load('flex')
    opt.load('bison')
        

def configure(conf):
    
    conf.load('frogs')
    
    # flex and bison
    conf.check_tool( 'flex' )
    conf.env['FLEXFLAGS'] = '' 
    conf.check_tool( 'bison' )
    conf.env['BISONFLAGS'] = ['-y', '-l', '-d']
    
    # various headers & libraries 
    conf.check_cc( header_name='math.h' )
    conf.check_cc( lib='z', msg='checking for zlib' )    
    
    # cernlib
    conf.check_cernlib()
    
    # requirements 
    conf.check_packages(requirements)    
    conf.write_config_header("config.h")


def build(bld):
    
    bld.add_subdirs( [ "src_base",
                       "src_utils", 
                       "src", 
                       "psf", 
                       "flat", 
                       "lc",
                       "mc",
                       "simphot",
                       "utils", 
                       "cern_stuff"] )
    
    if not bld.env.global_lapack:
        bld.add_subdirs(["lapack_stuff"])

    if Options.options.cernlib and bld.env.HAVE_CERN:
        bld.add_subdirs("cern_utils")

    # pkg-config file 
    gen_pkgconfig(bld)
