# -*- mode: python; -*- 


import os
import os.path as op
import sys
import Options
import Configure
import commands
import shlex

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
    


def parse_cernlib_flags(line,uselib,env,at_sys_name):
    """
    Special workaround to parse cernlib flags
    """
    
    lst=shlex.split(line)
    while lst:
        x=lst.pop(0)
        st=x[:2]
        ot=x[2:]
        if st=='-I'or st=='/I':
            if not ot:ot=lst.pop(0)
            env.append_unique('CPPPATH_'+uselib,ot)
        elif st=='-D':
            if not ot:ot=lst.pop(0)
            env.append_unique('CXXDEFINES_'+uselib,ot)
            env.append_unique('CCDEFINES_'+uselib,ot)
        elif st=='-l':
            if not ot:ot=lst.pop(0)
            env.append_unique('LIB_'+uselib,ot)
        elif st=='-L':
            if not ot:ot=lst.pop(0)
            env.append_unique('LIBPATH_'+uselib,ot)
        elif x=='-pthread'or x.startswith('+'):
            env.append_unique('CCFLAGS_'+uselib,x)
            env.append_unique('CXXFLAGS_'+uselib,x)
            env.append_unique('LINKFLAGS_'+uselib,x)
        elif x=='-framework':
            env.append_unique('FRAMEWORK_'+uselib,lst.pop(0))
        elif x.startswith('-F'):
            env.append_unique('FRAMEWORKPATH_'+uselib,x[2:])
        elif x.startswith('-std'):
            env.append_unique('CCFLAGS_'+uselib,x)
            env.append_unique('LINKFLAGS_'+uselib,x)
        elif x.startswith('-Wl'):
            env.append_unique('LINKFLAGS_'+uselib,x)
        elif x.startswith('-m')or x.startswith('-f'):
            env.append_unique('CCFLAGS_'+uselib,x)
            env.append_unique('CXXFLAGS_'+uselib,x)
        elif x.startswith('/'):
            t = x.replace('@sys', at_sys_name)
            env.append_unique('LIBPATH_'+uselib, op.dirname(t))
            
            libname = op.basename(x).replace('.a', '')
            if libname[0:3] == 'lib': libname = libname[3:]
            env.append_unique('STATICLIB_'+uselib, libname)



def configure( conf ):

    at_sys_name = None
    
    if conf.find_program('fs') and os.system('fs sys') == 0:
        ret = commands.getstatusoutput('fs sys')[1]
        at_sys_name = ret.split("'")[-2]
        targ = 'build.'+ at_sys_name
    else:
        ret = os.uname()
        targ = 'build.'+ ret[0]+'-'+ret[-1]
        try:
            os.remove('core')
        except OSError :
            pass
    
    
    conf.env.NAME=targ
    conf.env.set_variant(targ)
        
    # fortran compilers (not included by default in waf-1.5)
    conf.check_tool( 'compiler_fortran', tooldir='./wtools' )
    
    # c compiler 
    conf.check_tool( 'compiler_cc' )
    # TODO: why is -fPIC -DPIC not specified by default for gcc ? 
    #       is there a SharedObject() method like in scons ? 
    conf.env['CCFLAGS'] = [ '-fPIC', '-DPIC', '-g', '-O3' ]
    
    # c++ compiler 
    conf.check_tool( 'compiler_cxx' )
    conf.env.append_value( 'CXXFLAGS', ['-g','-O3'] );
    
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
#        ret = commands.getstatusoutput('cernlib mathlib pawlib')
        ret = commands.getstatusoutput('cernlib packlib')
        
        if ret[0] != 0:
            raise Configure.ConfigurationError
        
        line = ret[1]
        try:
            parse_cernlib_flags(line, 'cern', conf.env, at_sys_name)
        except:
            print sys.exc_info()
        
        conf.env.HAVE_CERN = 1
    except Configure.ConfigurationError:
        #        print 'CERNLIB not found.'        
        conf.env.HAVE_CERN = 0
        # for the moment, this is still a fatal error
        conf.fatal('cernlib not found. please install it first.')
        
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
                        package = 'sex-2.4.4', 
                        mandatory = True, 
                        uselib_store = 'SEX' )
        conf.env.LINKFLAGS_SEX += ['-Wl,-rpath']+conf.env.LIBPATH_SEX
    except Configure.ConfigurationError:
        conf.fatal('unable to locate sextractor.')
        
    # cfitsio 
    try:
        conf.check_cfg( args = '--cflags --libs', 
                        package = 'cfitsio-3.006', 
                        mandatory = True, 
                        uselib_store = 'CFITSIO' )
        conf.env.LINKFLAGS_CFITSIO += ['-Wl,-rpath']+conf.env.LIBPATH_CFITSIO
    except Configure.ConfigurationError:
        conf.fatal('unable to locate cfitsio')

    
    conf.env['POLOKA_VERSION'] = VERSION
    conf.write_config_header("config.h")


def build( bld ):
    
    bld.add_subdirs( [ "src_base",
                       "src_utils", 
                       "src", 
                       "psf", 
                       "flat", 
                       "lc",
                       "mc",
                       "utils" ] )
    
    if not bld.env.global_lapack:
        bld.add_subdirs(["lapack_stuff"])

    if bld.env.HAVE_CERN:
        bld.add_subdirs("cern_utils")
        bld.add_subdirs("cern_stuff")


    obj = bld( 'subst', 
               target = 'poloka.pc', 
               source = 'poloka.pc.in' )
    
    obj.dict = { 'PREFIX': bld.env['PREFIX'], 
                 'VERSION': bld.env['POLOKA_VERSION'] }
    
    bld.install_as( '${PREFIX}/lib/pkgconfig/poloka-${POLOKA_VERSION}.pc', 
                    'poloka.pc' )
    
