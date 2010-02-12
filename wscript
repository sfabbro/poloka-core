# -*- mode: python; -*- 


import os
import sys
sys.path.append('./wtools')
import Options
import Configure


APPNAME  = 'poloka'
VERSION  = '0.1.0'
srcdir   = '.'
blddir   = 'build'



def set_options( ctx ):
    
    ctx.tool_options('compiler_cc')
    ctx.tool_options('compiler_cxx')
    ctx.tool_options('compiler_fortran', tooldir='./wtools')
    ctx.tool_options('flex')
    ctx.tool_options('bison')


def configure( conf ):
    
    # fortran compilers (not included by default in waf-1.5)
    conf.check_tool( 'compiler_fortran', tooldir='./wtools' )
    
    # c compiler 
    conf.check_tool( 'compiler_cc' )
    # TODO: why is -fPIC -DPIC not specified by default for gcc ? 
    #       is there a SharedObject() method like in scons ? 
    conf.env['CCFLAGS'] = [ '-fPIC', '-DPIC' ]
    
    # c++ compiler 
    conf.check_tool( 'compiler_cxx' )
    
    # flex and bison
    conf.check_tool( 'flex' )
    conf.env['FLEXFLAGS'] = '' 
    conf.check_tool( 'bison' )
    conf.env['BISONFLAGS'] = ['-y', '-l', '-d']
    
    # file substitutions
    conf.check_tool( 'misc' )
    
    # various headers & libraries 
    conf.check_cc( header_name='math.h' )
    conf.check_cc( lib='z', msg='checking for zlib' )    
    
    # cernlib 
    try:
        conf.find_program( 'cernlib', mandatory = True )
        conf.check_cfg( path='cernlib', args='', 
                        package='mathlib packlib',
                        uselib_store='cern' )
    except Configure.ConfigurationError:
        conf.fatal('CERNLIB not found.')
        
        
    # pkg config 
    conf.find_program('pkg-config')
    
#     # sextractor 
#     try:
#         conf.check_cfg( path='sex-config',
#                         args = '--cflags --libs', 
#                         package = 'sextractor', 
#                         mandatory = True )
#     except Configure.ConfigurationError:
#         conf.fatal('unable to locate sextractor.')
        
#     # cfitsio 
#     try:
#         conf.check_cfg( path='cfitsio-config', 
#                         args = '--cflags --libs', 
#                         package = 'cfitsio', 
#                         mandatory = True )
#     except Configure.ConfigurationError:
#         conf.fatal('unable to locate cfitsio')
        
    


def build( bld ):
    
    bld.add_subdirs( [ "lapack_stuff", 
                       "cern_stuff",
                       "src_base",
                       "src_utils", 
                       "cern_utils",
                       "src", 
                       "psf", 
                       "flat", 
                       "lc",
                       "utils" ] )
    
    # if the thing is not installed 
    
    obj = bld( 'subst', 
               target = 'poloka.pc', 
               source = 'poloka.pc.in' )
    
    obj.dict = { 'NAME': APPNAME, 
                 'VERSION': VERSION, 
                 'REQUIREMENTS': '',
                 'CONFLICTS': '',
                 'LIBS': bld.env.LIBS, 
                 'LIBS_PRIVATE': '', 
                 'CFLAGS': bld.env.CFLAGS }

#     dd = obj.env.get_merged_dict()
#     import pprint
#     pprint.pprint( dd )
#     print bld.path.__class__
