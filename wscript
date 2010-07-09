# -*- mode: python; -*- 


import os
import os.path as op
import sys
import Options
import Configure
import commands

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
    ctx.add_option('--with-cfitsio', action='store', help='Path to the cfitsio root dir')
    ctx.add_option('--with-sex', action='store', help='Path to the cfitsio root dir')
    
    ctx.add_option('--no-global-lapack', action='store_false',
                   default=True, dest='global_lapack',
                   help='use/no do use system lapack library')
    

def configure( conf ):
    
    if conf.find_program('fs'):
        ret = commands.getstatusoutput('fs sys')[1]
        targ = 'build.'+ ret.split("'")[-2]
    else:
        ret = os.uname()
        targ = 'build.'+ ret[0]+'-'+ret[-1]
    
    
    conf.env.NAME=targ
    conf.env.set_variant(targ)
        
    # fortran compilers (not included by default in waf-1.5)
    conf.check_tool( 'compiler_fortran', tooldir='./wtools' )
    
    # c compiler 
    conf.check_tool( 'compiler_cc' )
    # TODO: why is -fPIC -DPIC not specified by default for gcc ? 
    #       is there a SharedObject() method like in scons ? 
    conf.env['CCFLAGS'] = [ '-fPIC', '-DPIC', '-g' ]
    
    # c++ compiler 
    conf.check_tool( 'compiler_cxx' )
    conf.env.append_value( 'CXXFLAGS', ['-g'] );
    
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
    
    # lapack
    if Options.options.global_lapack:
        if conf.check_cc( lib='lapack', msg='checking for lapack' ):
            conf.env.global_lapack = True
        else:
            conf.env.global_lapack = False
    
    #    conf.env.global_lapack = False
    
    # cernlib 
    try:
        conf.find_program( 'cernlib', mandatory = True )
        conf.check_cfg( path='cernlib', args='', 
                        package='mathlib pawlib',
                        #                        package='mathlib packlib',
                        uselib_store='cern' )
        conf.env.HAVE_CERN = 1
    except Configure.ConfigurationError:
        print 'CERNLIB not found.'
        conf.env.HAVE_CERN = 0
        
    # pkg config 
    try:
        conf.find_program('pkg-config')
        pkgcpath = op.join( Options.options.prefix, 'lib', 'pkgconfig')
        if os.environ.has_key('PKG_CONFIG_PATH'):
            os.environ['PKG_CONFIG_PATH'] = pkgcpath + ":" + os.environ['PKG_CONFIG_PATH']
        else:
            os.environ['PKG_CONFIG_PATH'] = pkgcpath
    except Configure.ConfigurationError:
        conf.fatal('pkg-config not found')
    
    # sextractor 
    try:
        conf.check_cfg( args = '--cflags --libs', 
                        package = 'sex-2.2.2', 
                        mandatory = True, 
                        uselib_store = 'SEX' )
    except Configure.ConfigurationError:
        conf.fatal('unable to locate sextractor.')
        
    # cfitsio 
    try:
        conf.check_cfg( args = '--cflags --libs', 
                        package = 'cfitsio-3.0.0', 
                        mandatory = True, 
                        uselib_store = 'CFITSIO' )
    except Configure.ConfigurationError:
        conf.fatal('unable to locate cfitsio')
        
    
    conf.env['POLOKA_VERSION'] = VERSION
    conf.write_config_header("config.h")


def build( bld ):
    
    bld.add_subdirs( [ "cern_stuff",
                       "src_base",
                       "src_utils", 
                       #                       "cern_utils",
                       "src", 
                       "psf", 
                       "flat", 
                       "lc",
                       "mc",
                       "utils" ] )
    
    if not bld.env.global_lapack:
        bld.add_subdirs(["lapack_stuff"])


    obj = bld( 'subst', 
               target = 'poloka.pc', 
               source = 'poloka.pc.in' )
    
    obj.dict = { 'PREFIX': bld.env['PREFIX'], 
                 'VERSION': bld.env['POLOKA_VERSION'] }
    
    bld.install_as( '${PREFIX}/lib/pkgconfig/poloka-${POLOKA_VERSION}.pc', 
                    'poloka.pc' )
    
