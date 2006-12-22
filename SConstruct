#!/usr/bin/env python
# 
# main SConstruct file
# 
# MAIN TARGETS:
#   scons config             run the autoconf tests
#   scons build  [default]   build the distribution
#   scons install            install the files into prefix/
#   scons dist               build a dist tarball
#
# OPTIONS:
#    builddir=<buildir>      build the files into <builddir>
#    prefix=<prefix>         the files are to be installed into prefix
#                            NOTE: the RPATH is positioned accordingly
#    debug=0|1               if debug=1, compile with -g -Wall
#    relase=0|1              if release=1, compile with -O3
#    frogshome=<path>        search the external dependencies
#                            (cfitsio, sex...) in frogshome
#    
import os
import sys
import re
import string
import fnmatch
import commands

sys.path.append('.' + os.sep + 'scons-local')
sys.path.append('.' + os.sep + 'scons-local' + os.sep + 'scons-local-0.96.90')
sys.path.append('/usr/lib/scons/')


import SCons
from SCons.Script.Main import OptParser
import config.Package as package


# version
packageVersion = (0, 1, 0)
packageName = 'poloka'
fullPackageName = packageName + '-%d.%d.%d' % packageVersion


# options
opts = Options('options.cache')
opts.AddOptions(
              ('builddir','specify the build directory', 'build'),
              ('prefix', 'Installation Prefix', '.'),
    BoolOption('debug', 'Compile a debug version', 1),
    BoolOption('release', 'release build (optimize + install)', 0),
    PathOption('frogshome', 'frogs software installation prefix', '/usr/local'),

    
    PathOption('boost_prefix', 'boost library installation prefix', '/usr/local/lib'),
    BoolOption('enable_persistence', 'compile w/ the boost persistence', 1),
    BoolOption('fetch', 'fetch automatically missing software', 1),
    BoolOption('distcc', 'try to compile the software using distcc', 0),
    BoolOption('ccache', 'try to compile the software using ccache', 1),
    BoolOption('rpath', 'enable building libraries w/ the --rpath option (gcc only)', 1),
    
    BoolOption('cvs_cfitsio', 'always retrieve cfitsio from the frogs CVS server', 0),
    BoolOption('cvs_sex', 'always retrieve sextractor from the frogs CVS server', 0),
    BoolOption('cvs_lapack', 'always retrieve LAPACK from the frogs CVS server', 0)
    )


# platform
arch = os.uname()[4]
os_name = os.uname()[0]
platform = sys.platform
Export('arch os_name platform')


# environment
g_env = Environment(options=opts, ENV=os.environ)
g_env.Alias('config')
g_env.Alias('build')
g_env.Alias('install')
g_env.Alias('dist')



# targets ?
parser = OptParser()
options,targets = parser.parse_args(sys.argv[1:])


def dictgen(target, source, env):
    """
    Builder action.
    Builds a dictionary (calls dictgen2) from a dictfile...
    """
    dg = SCons.Node.FS.default_fs.File('#/%s/objio/dictgen' % (env['builddir'])).abspath
    tgpath = os.path.dirname(SCons.Node.FS.default_fs.File(target[0]).abspath)
    cmd = "%s -o %s %s" % (dg, target[0], source[0])
    print 'About to run: ', cmd
    res = commands.getstatusoutput(cmd)
    if res[0] != 0:
        print 'Unable to create ', str(target[0]), ' from ', str(source[0])
        return
    

bld = g_env.Builder(action=dictgen,
                    suffix='_dict.cc',
                    src_suffix='.h')
g_env.Append(BUILDERS = {'Dictgen': bld})


# distribution
distbase = 'poloka-0.1.0.tar.gz'
g_env['ENV']['distbase'] = '#'  + distbase
g_env.Alias('dist', g_env['ENV']['distbase'])
g_env.AddToDist(['SConstruct',
                 'AUTHORS', 'ChangeLog', 'COPYING',
                 'INSTALL', 'README', 'NEWS', 'TODO',
                 'configure.ac', 'Makefile.am',
                 'config.h'])

if g_env.WhereIs('ccache') != None:
    g_env.Replace(CC = 'ccache gcc')
    g_env.Replace(CXX = 'ccache g++')


# we should protect all of this ... 
g_env.Append(CPPFLAGS=['-Wall'])
if g_env['debug'] == 1:
    g_env.Append(CPPFLAGS=['-g'])
    g_env.Append(LDFLAGS=['-g'])
if g_env['release'] == 1:
    g_env.Append(CPPFLAGS='-O3')
    g_env.Append(CFLAGS='-03')
    g_env.Append(F77FLAGS='-O3')


# autoconf section
if (not os.path.isfile('config.log')) or \
    'config' in targets:
    
    config_env = g_env.Copy()
    
    cern = package.Cernlib()
    if not cern.locate():
        print 'Unable to find the CERN LIBS'
        print 'Calling cern.fetch()'
        print 'Calling cern.build()'
    g_env.RegisterExternalPackage(cern)
    config_env.TestExternalPackage(cern)
    
    sex = package.SExtractor(frogshome=g_env['frogshome'])
    if not sex.locate():
        print 'Unable to find the -lsex and -lwcs libraries'
        print 'Calling sex.fetch()'
        print 'Calling sex.build()'
    g_env.RegisterExternalPackage(sex)
    config_env.TestExternalPackage(sex)
    
    cfitsio = package.Cfitsio(frogshome='/home/nrl/software/cfitsio')
    if not cfitsio.locate():
        print 'Unable to find the -lcfitsio library'
        print 'Calling cfitsio.fetch()'
        print 'Calling cfitsio.build()'
    print cfitsio.includepath
    print cfitsio.libpath
    print cfitsio.libs
    g_env.RegisterExternalPackage(cfitsio)
    config_env.TestExternalPackage(cfitsio)
    
    config = config_env.Configure()
    
    # note: to detect commands: use the WhereIs method...
    # headers 
    config.CheckCXXHeader('math.h')
    config.CheckCXXHeader('vector')
    config.CheckCXXHeader('list')
    config.CheckCXXHeader('string')
    config.CheckCXXHeader('unistd.h')
    
    # libmath 
    if not config.CheckLib('m'):
        print '   *** Error not able to locate the C math library'
        print '   *** Check your C installation.\n'
        config_env.Exit(1)

    # libg2c
    if not config.CheckLib('g2c'):
        print '   *** Error: not able to locate the g2c library'
        print '   *** check your gcc installation.\n'
        config_env.Exit(1)
        
    # Now, less common packages
    if not config.CheckLib('wcs') or not config.CheckLib('sex') or g_env['cvs_sex']:
        print '   *** sextractor does not seem to be installed on your system'
        config_env.Exit(1)
        #        g_env.SConscript('sex' + os.sep + 'SConscript')
        #        g_env.SourceCode('sex', None)
        
    if not config.CheckLib('cfitsio') or g_env['cvs_cfitsio'] == 1:
        print '   *** libcfitsio does not seem to be installed on your system'
        config_env.Exit(1)
        #        g_env.SConscript('cfitsio' + os.sep + 'SConscript')
        #        g_env.SourceCode('cfitsio', None)
        
    # The CERN libraries -- check also the $CERN variable
    #    if not g_env.WhereIs('cernlib'):
    #        print g_env['ENV']['PATH']
    #        print '   *** the CERN library does not seem to be installed on your system'
    #        g_env.Exit(1)
    #    else:
    #        # here, I should do this properly (add a -L//// + the -l...)
    #        g_env.ParseConfig('cernlib')
    if not config.CheckLib(['packlib', 'g2c']):
        print '   *** The CERN library does not seem to be installed on your system'
        config_env.Exit(1)
    
    
    #    # LAPACK and BLAS
    #    if not config.CheckLib('blas') or not config.CheckLib('lapack') or g_env['cvs_lapack']:
    #        print '   *** LAPACK/BLAS does not seem to be installed on your system'
    #        print '   *** getting it from the frogs CVS server ?'
    
    config_env = config.Finish()
    g_env.SaveExternalPackageInfo('extpackages.pkl')
    
else:
    if os.path.isfile('extpackages.pkl'):
        g_env.LoadExternalPackageInfo('extpackages.pkl')

# get rid of this.
g_env.Append(LIBS = ['g2c'])
        

# Help and options
Help(opts.GenerateHelpText(g_env))
opts.Save('options.cache', g_env)


# files and subdirs
SUBDIRS = ['objio',
           'lapack_stuff',
           'cern_stuff',
           'src_base',
           'src_utils',
           'src',
           'dao',
           'cern_utils',
           'flat',
           'lc',
           'mc',
           'utils',
           'telinst',
           'm4',
           'doc',
           'datacards',
           'scons-local']

EXTRA_DIST = ['aclocal.m4',
              'AUTHORS',
              'autogen.sh',
              'ChangeLog',
              'configure.ac',
              'COPYING',
              'INSTALL',
              'Makefile.am',
              'NEWS',
              'README',
              'SConstruct',
              'TODO']


# Builds
build_dir = "build-%s-%s" % (arch, os_name)
libs = []
headers = []
binaries = []

Export('g_env')


# build dirs
for s in SUBDIRS:
    g_env.BuildDir(g_env['builddir'] + os.sep + s, s, duplicate=0)


# targets
g_env.Alias('build', [g_env['builddir'] + os.sep + s for s in SUBDIRS])
Default('build')


# build the submodules
for s in [g_env['builddir'] + os.sep + s + os.sep + 'SConscript' for s in SUBDIRS]:
    print 'Examining %s' % (s)
    g_env.SConscript(s)



#    config.CheckFunc('getcwd')
#    config.CheckFunc('memchr')
#    config.CheckFunc('memset')
#    config.CheckFunc('pow')
#    config.CheckFunc('rint')
#    config.CheckFunc('setenv')
#    config.CheckFunc('sqrt')
#    config.CheckFunc('strchr') 
#    config.CheckFunc('strrchr')
#    config.CheckFunc('strspn')
#    config.CheckFunc('strstr')



