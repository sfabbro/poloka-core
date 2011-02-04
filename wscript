# -*- mode: python; -*- 


import os
import os.path as op
import sys
import commands
import shlex

import Options
import Configure
import frogsutils 


APPNAME  = 'poloka'
VERSION  = '0.1.0'
top   = '.'
out   = frogsutils.get_out_name()
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
    conf.load( 'flex' )
    conf.env['FLEXFLAGS'] = '' 
    conf.load( 'bison' )
    conf.env['BISONFLAGS'] = ['-y', '-l', '-d']
    
    # doxygen ? 
    conf.find_program('doxygen', VAR='DOXYGEN', 
                      mandatory=False)
    
    # various headers & libraries 
    conf.check_cc(header_name='math.h')
    conf.check_cc(header_name='fenv.h')
    conf.check_cc(fragment="""#define _GNU_SOURCE
#include <fenv.h>
int main() {
void *p;
p=(void*)(feenableexcept);
return 0;}""", 
                  define_name="HAVE_FEENABLEEXCEPT", 
                  msg = "Checking for feenableexcept",
                  mandatory=False)
    conf.check_cc(lib='z', msg='Checking for zlib')    
    
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
    
    if bld.env['DOXYGEN']:
        bld.add_subdirs("doc");

    # pkg-config file 
    gen_pkgconfig(bld)

