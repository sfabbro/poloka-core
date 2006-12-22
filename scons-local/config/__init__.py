#!/usr/bin/env python

import os
import sys
import fnmatch
import string
import shutil
import distutils.archive_util as arch
import pickle


if os.path.isdir('/usr/lib/scons'):
    sys.path.append('/usr/lib/scons')
elif os.path.isdir('/usr/local/lib/scons'):
    sys.path.append('/usr/local/lib/scons')

import SCons
import SCons.Script
from SCons.Util import Split
from SCons.Defaults import Chmod


def Glob(includes=Split('*'), excludes=None, dir='.'):
    """
    Glob function. Similar to glob.glob
    But here, we glob SCons.Node objects instead of simple file names,
    in order to be able to safely build the software in a separate
    build dir...
    """
    
    def filename_filter(node):
        fn = os.path.basename(str(node))
        match = 0
        for p in includes:
            if fnmatch.fnmatchcase(fn, p):
                match=1
                break
        if match == 1 and excludes!=None:
            for p in excludes:
                if fnmatch.fnmatchcase(fn,p):
                    match = 0
                    break
        return match
    
    
    def filter_nodes(nodelist):
        """
        Filter all the nodes by name.
        Return a list of SCons.Node (Files or Dirs)
        """
        children = filter( filename_filter, nodelist )
        nodes = []
        for c in children:
            nodes.append( c )
        return nodes

    
    def gen_node(nm):
        """
        Convert the filename into a real Node.
        """
        if type(nm) in (type(''), type(u'')):
            path = nm
        else:
            path = n.abspath

        if os.path.isdir(path):
            return SCons.Node.FS.default_fs.Dir(path)
        else:
            return SCons.Node.FS.default_fs.File(path)
    
    
    from SCons.Scanner.Dir import DirScanner
    scanner = DirScanner()
    def scanDir(dir):
        nodes = []
        nodes.extend(dir.all_children())
        nodes.extend(scanner(dir,None))
        return nodes

    here = SCons.Node.FS.default_fs.Dir(dir)
    nodes = scanDir(here)
    nodes = filter_nodes(nodes)
    nodes_srcs = [n.srcnode() for n in nodes]
    
    src = here.srcnode()
    if not src is here:
        subdirs = scanDir(src)
        for s in filter_nodes(subdirs):
            if not s in nodes_srcs:
                nodes.append( gen_node(os.path.join(dir, os.path.basename(str(s)))) )


    for n in nodes_srcs:
        print '  --> ', n
        
    return nodes




def Distribute(target, source, env):
    """
    Builder action.
    Builds an archive from the list of distributed files.
    TODO:
      * archive format
    """
    
    tmpdir = str(target[0])
    tmpdir = string.replace(tmpdir, '.tar.gz', '')
    target = tmpdir
    
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)

    for s in source:
        s = str(s)
        s_name = os.path.basename(s)
        s_path = tmpdir + os.sep + os.path.dirname(s)
        if not os.path.exists(s_path):
            os.makedirs(s_path)
        shutil.copy2(s, os.path.join(s_path, s_name))
        
    archive_name = arch.make_archive(target, 'gztar', root_dir='.', base_dir=tmpdir)
    shutil.rmtree(tmpdir)




def AddInternalDep(env, dirname, libname=None):
    """
    Wrapper function.
    Declares another subpackage we depend from.

    Should we use the Dir() stuff ?
    """
    
    if env.has_key('prefix'):
        prefix = env['prefix']
    else:
        prefix = '#'
        
    if env.has_key('builddir'):
        builddir = env['builddir']
    else:
        builddir = ''

    if env.has_key('libdir'):
        libdir = env['libdir']
    else:
        libdir = 'lib'
        
    # the -I option. 
    if dirname != '.':
        env.Append(CPPPATH = ['#' + dirname])
    else:
        env.Append(CPPPATH = ['.'])

    if libname != None:
        # the -l stuff
        env.Append(LIBS    = [libname])
        
        # the -L options
        if builddir != None:
            env.Append(LIBPATH = ['#' + builddir + os.sep + dirname])
        else:
            env.Append(LIBPATH = ['#' + dirname])
        
        if prefix == '.':
            if builddir != None:
                env.Append(RPATH   = [SCons.Node.FS.default_fs.Dir('#' + builddir + os.sep + dirname).abspath])
            else:
                env.Append(RPATH   = [SCons.Node.FS.default_fs.Dir('#' + dirname).abspath])
        else:
            env.Append(RPATH   = [SCons.Node.FS.default_fs.Dir(prefix + os.sep + libdir).abspath])



def AddExternalDep(env, includes, libpath, libs):
    """
    REMOVE: not used anymore.
    Wrapper function.
    Declare an external package we depend from.
    """
    if includes != None:
        env.Append(CPPPATH = includes)
    if libpath != None:
        env.Append(LIBPATH = libpath)
        env.Append(RPATH   = libpath)
    if libs != None:
        env.Append(LIBS    = libs)

        
def RegisterExternalPackage(env, pkg):
    """
    Wrapper function.
    Register an external package in the Environment.
    """
    if not env['ENV'].has_key('EXT-PACKAGES'):
        env['ENV']['EXT-PACKAGES'] = {}
        
    name = pkg.name
    if not env['ENV']['EXT-PACKAGES'].has_key(name):
        d = {}
        env['ENV']['EXT-PACKAGES'][name] = d
        
    d['INCLUDEPATH'] = pkg.includepath
    d['LIBPATH'] = pkg.libpath
    d['LIBS'] = pkg.libs
    

def ActivateExternalPackage(env, pkgname):
    """
    Wrapper function.
    Add the package related flags to the compile and link flags.
    """
    
    if not env['ENV'].has_key('EXT-PACKAGES'):
        return
    
#    import pprint
#    pprint.pprint(env['ENV']['EXT-PACKAGES'])
    
    if not env['ENV']['EXT-PACKAGES'].has_key(pkgname):
        print 'Did not found ', pkgname
        return
    
    d = env['ENV']['EXT-PACKAGES'][pkgname]
    
    if len(d['INCLUDEPATH']) > 0:
        env.Append(CPPPATH = d['INCLUDEPATH'])
    if len(d['LIBPATH']) > 0:
        env.Append(LIBPATH = d['LIBPATH'])
        env.Append(RPATH = d['LIBPATH'])
    if len(d['LIBS']) > 0:
        env.Append(LIBS = d['LIBS'])

    
    #    print env['LIBS']
    
        
def TestExternalPackage(env, pkg):
    """
    Wrapper function.
    Register the package build flags, but not the LIBS.
    """
    
    if len(pkg.includepath) > 0:
        env.Append(CPPPATH = pkg.includepath)
    if len(pkg.libpath) > 0:
        env.Append(LIBPATH = pkg.libpath)
        env.Append(RPATH   = pkg.libpath)

#    if len(d['LIBS']) > 0:
#        env.Append(LIBS = d['LIBS'])

    

def SaveExternalPackageInfo(env, filename):
    """
    Wrapper function.
    Retrieve the external package information from the specified file.
    """
    f = open(filename, 'w')
    if env['ENV'].has_key('EXT-PACKAGES'):
        pickle.dump(env['ENV']['EXT-PACKAGES'], f)
    f.close()


def LoadExternalPackageInfo(env, filename):
    """
    Wrapper function.
    Save the external package information to the specified file.
    """
    f = open(filename)
    try:
        epi = pickle.load(f)
    except:
        f.close()
        return
    env['ENV']['EXT-PACKAGES'] = epi
    f.close()
    
        
def ExternalPackage(env, package, config=False):
    """
    REMOVE THIS.
    Wrapper function.
    Declare an external package.
    """
    if len(package.includepath) > 0:
        env.Append(CPPPATH = package.includepath)
    if len(package.libpath) > 0:
        env.Append(LIBPATH = package.libpath)
        env.Append(RPATH   = package.libpath)
    if not config and len(package.libs) > 0:
        env.Append(LIBS = package.libs)
    
        
def AddToDist(env, sources):
    """
    Wrapper function.
    Add files to the dist targets.
    """
    if not env['BUILDERS'].has_key('Distribute'):
        env['BUILDERS']['Distribute'] = env.Builder(action=Distribute,
                                                    multi=1)
        
    if not env['ENV'].has_key('distbase'):
        env['ENV']['distbase'] = 'archive'
        
    env.Distribute(env['ENV']['distbase'], sources)
    
    
def InstallProgs(env, progs):
    """
    Wrapper function.
    Defines an install target for the program
    and Chmod() it to a reasonnable value.
    """
    if not env.has_key('prefix')  or env['prefix'] == '.':
        return None
    prefix = env['prefix']

    if not env.has_key('bindir'):
        bindir = 'bin'
    else:
        bindir = env['bindir']

    #    env.Alias('install', prefix + os.sep + bindir)
    inst = env.Install(prefix + os.sep + bindir, progs)
    env.AddPostAction(inst, Chmod('$TARGETS', 0755))
    return inst


def InstallLibs(env, libs):
    """
    Wrapper function.
    Defines an install target for the library
    and Chmod() it to a reasonnable value.
    """
    if not env.has_key('prefix') or env['prefix'] == '.':
        return None
    prefix = env['prefix']
    
    if not env.has_key('libdir'):
        libdir = 'lib'
    else:
        libdir = env['libdir']

    #    env.Alias('install', prefix + os.sep + libdir)
    inst = env.Install(prefix + os.sep + libdir, libs)
    env.AddPostAction(inst, Chmod('$TARGETS', 0755))
    return inst



def InstallHeaders(env, headers):
    """
    Wrapper function.
    Defines an install target for the headers
    and Chmod() them to a reasonnable value.
    """
    if not env.has_key('prefix')  or env['prefix'] == '.':
        return None
    prefix = env['prefix']
    
    if not env.has_key('includedir'):
        includedir = 'include'
    else:
        includedir = env['includedir']

    #    env.Alias('install', prefix + os.sep + includedir)
    inst = env.Install(prefix + os.sep + includedir, headers)
    env.AddPostAction(inst, Chmod('$TARGETS', 0644))
    return inst

    
def InstallData(env, data):
    """
    Wrapper function.
    Defines an install target for the data
    and Chmod() them to a reasonnable value.
    """
    if not env.has_key('prefix')  or env['prefix'] == '.':
        return None
    prefix = env['prefix']
    
    if not env.has_key('datadir'):
        # TODO: poloka should be replace by something else.
        datadir = 'share' + os.sep + 'poloka'
    else:
        datadir = env['datadir']
    
    env.Alias('install', prefix + os.sep + datadir)
    inst = env.Install(prefix + os.sep + datadir, data)
    env.AddPostAction(inst, Chmod('$TARGETS', 0644))
    return inst



from SCons.Script.SConscript import SConsEnvironment
SConsEnvironment.AddToDist       = AddToDist
SConsEnvironment.AddInternalDep  = AddInternalDep
SConsEnvironment.AddExternalDep  = AddExternalDep
SConsEnvironment.InstallProgs    = InstallProgs
SConsEnvironment.InstallHeaders  = InstallHeaders
SConsEnvironment.InstallLibs     = InstallLibs
SConsEnvironment.InstallData     = InstallData
SConsEnvironment.ExternalPackage = ExternalPackage
SConsEnvironment.RegisterExternalPackage = RegisterExternalPackage
SConsEnvironment.ActivateExternalPackage = ActivateExternalPackage
SConsEnvironment.TestExternalPackage = TestExternalPackage
SConsEnvironment.SaveExternalPackageInfo = SaveExternalPackageInfo
SConsEnvironment.LoadExternalPackageInfo = LoadExternalPackageInfo
