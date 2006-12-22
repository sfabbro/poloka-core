#!/usr/bin/env python

import re
import os
import stat
import glob
import string
import fnmatch
import commands


class Package:
    """
    Base class for the external packages.
    """
    
    def __init__(self, name, frogshome=None):
        """
        Default Constructor.
        """
        self.name = name
        
        self.progpattern = []
        self.libpattern = []
        self.headerpattern = []
        
        self.progs = []
        self.progsfound = []
        self.libs = []
        self.headers = []
        
        self.prefixsearchpath = None
        self.binsearchpath = []
        self.libsearchpath = []
        self.includesearchpath = []
        
        self.binpath = []
        self.libpath = []
        self.includepath = []
        
        if frogshome != None:
            self.prefixsearchpath = frogshome
            
        self.libsearchpath.append('/usr/local/lib')
        self.libsearchpath.append('/usr/lib')
        self.libsearchpath.append('/lib')
        
        self.includesearchpath.append('/usr/local/include')
        self.includesearchpath.append('/usr/include')

        self.binsearchpath.append('/usr/local/bin')
        self.binsearchpath.append('/usr/bin')
        self.binsearchpath.append('/bin')
        
        
    def locate(self):
        """
        Attempts to locate the package features
        at the usual places
        """
        h = self.locateHeaders()
        l = self.locateLibs()
        p = self.locatePrograms()
        
        return h or l or p
    
    
    def locateFilesInPath(self, names, mode, path):
        """
        Try to locate the specified file with the
        specified (octal) mode in the specified path.
        """
        
        ret = []
        for n in names:
            for d in path:
                f = os.path.join(d, n)
                detect = glob.glob(f)
                for f in detect:
                    if os.path.isfile(f):
                        try:
                            st = os.stat(f)
                        except:
                            continue
                        if stat.S_IMODE(st[stat.ST_MODE]) & mode:
                            ret.append(f)
                            
        return ret
    
    
    def locateFilesUnderPrefix(self, patterns, mode, prefix):
        """
        Try to locate the specified file with the specified
        (octal) mode, under the specified prefix.
        """
        
        ret = []
        
        def match(arg, dir, files):
            patterns = arg[0]
            matches = arg[1]
            # mode ???
            for f in files:
                if not os.path.isfile(dir + os.sep + f): continue
                for p in patterns:
                    if not fnmatch.fnmatchcase(f,p): continue
                    ret.append(dir + os.sep + f)

        os.path.walk(prefix, match, (patterns, ret) )
        
        return ret
        
        # the os.walk function seems to be too news...
        #        for root,dirs,files in os.path.walk(prefix):
        #            for f in files:
        #                for p in patterns:
        #                    if fnmatch.fnmatchcase(f, p):
        #                        ret.append(os.path.join(root,f))
        #        return ret
    
    
    
    def locateHeaders(self):

        headers = []
        if self.prefixsearchpath != None and len(self.prefixsearchpath) > 0:
            headers = self.locateFilesUnderPrefix(self.headerpattern, 0444,
                                                  self.prefixsearchpath)

        if len(headers) == 0:
            headers = self.locateFilesInPath(self.headerpattern, 0444,
                                             self.includesearchpath)

        if len(headers) == 0:
            return False
        
        for h in headers:
            d = os.path.dirname(h)
            if not d in self.includepath:
                self.includepath.append(d)
                
        return True
            
    def locateLibs(self):
        
        libs = []
        if self.prefixsearchpath != None and len(self.prefixsearchpath) > 0:
            libs = self.locateFilesUnderPrefix(self.libpattern, 0111,
                                               self.prefixsearchpath)
            
        if len(libs) == 0:
            libs = self.locateFilesInPath(self.libpattern, 0111,
                                          self.libsearchpath)
            
        if len(libs) == 0:
            return False
        
        exp = re.compile('lib([\w\d-]+)\..*')
        for l in libs:
            d = os.path.dirname(l)
            if not d in self.libpath:
                self.libpath.append(d)
            m = exp.search(os.path.basename(l))
            if m == None:
                return False
            self.libs.append(m.group(1))
        
        return True

    
    def locatePrograms(self):
        
        progs = []
        if self.prefixsearchpath != None and len(self.prefixsearchpath) > 0:
            progs = self.locateFilesUnderPrefix(self.progpattern, 0111,
                                                self.prefixsearchpath)
            
        if len(progs) == 0:
            progs = self.locateFilesInPath(self.progpattern, 0111,
                                           self.binsearchpath)
            
        if len(progs) == 0:
            return False
        
        for l in progs:
            d = os.path.dirname(l)
            if not d in self.binpath:
                self.binpath.append(d)
            b = os.path.basename(l)
            if not b in self.progsfound:
                self.progsfound.append(b)
                self.progs.append(l)
                
        return True

    
    def fetch(self, url):
        """
        Attempts to retrieve the package.
        """
        pass
    
    
    def build(self):
        """
        Extract and build the package.
        """
        pass



class Cernlib(Package):
    """
    Attempt to find the CERN library location.
    Usually hidden behind the $CERN environment variable.
    """
    
    def __init__(self):
        Package.__init__(self, 'cern')



    def parseConfig(self, config):

        ret = {
            'CPPPATH': [],
            'LIBPATH': [],
            'CPPFLAGS': [],
            'LDFLAGS': [],
            'LIBS': []
            }
        
        fields = string.split(config)
        for f in fields:
            f = string.strip(f)
            if f[0:1] != '-':
                if f[0:1] != '/':
                    ret['LIBS'].append(f)
                else:
                    ret['LIBPATH'].append(os.path.dirname(f))
                    ret['LIBS'].append(os.path.basename(f))
                
            elif f[0:2] == '-L':
                ret['LIBPATH'].append(f[2:])
            elif f[0:2] == '-l':
                ret['LIBS'].append(f[2:])

        return ret
                
        
        
    def locate(self):

        print 'Cernlib.locate'
        
        env = os.environ
        if env.has_key('CERN'):
            self.prefixsearchpath = env['CERN']
            cernlib = self.locateFilesUnderPrefix(['cernlib'], 0111,
                                                  self.prefixsearchpath)
            print 'CERNLIB = ', cernlib
            
        if len(cernlib) == 0:
            cernlib = self.locateFilesInPath('cernlib', 0111,
                                             string.split(env['PATH'], os.pathsep))
            
        if len(cernlib) == 0:
            return False
            
        cernlib_cmd = cernlib[0]
        ret = commands.getstatusoutput(cernlib_cmd)
        if ret[0] == 0:
            cernflags = ret[1]
        else:
            return False

        print 'Cernlib.locate: cernlib output: ', cernflags

        ret = self.parseConfig(cernflags)

        self.libs.extend(ret['LIBS'])
        self.libpath.extend(ret['LIBPATH'])
        
        
#        cernflags = string.replace(cernflags, '-l', 'lib')
#        cernflags = string.split(cernflags)
#        libpath = []
#        libs = []
#        for c in cernflags:
#            if c[0] == os.sep:
#                p = os.path.dirname(c)
#                if not p in libpath:
#                    libpath.append(p)
#                l = os.path.basename(c)
#                libs.append(string.split(l,'.')[0])
#            else:
#                libs.append(string.split(c,'.')[0])
#                
#        self.libs.extend(libs)
#        self.libpath.extend(libpath)

        print self.libs
        print self.libpath
        
        return True
    
    
    
    def fetch(self):
        pass

    
    def build(self):
        pass

    
    
class Cfitsio(Package):
    """
    The CFITSio package.
    Usually hidden down below the frogshome package.
    """
    
    def __init__(self, frogshome=None):
        """
        Constructor. Just specify the frogshome.
        """
        
        Package.__init__(self, 'cfitsio', frogshome=frogshome)
        # get rid of this !
        if os.path.isdir(frogshome + os.sep + 'cfitsio' + os.sep + 'v2r440'):
            self.prefixsearchpath = frogshome + os.sep + 'cfitsio' + os.sep + 'v2r440'
            
        # CHECK: env methods --> library prefix
        self.libpattern = ['libcfitsio.so']
        self.headerpattern = ['fitsio.h', 'fitsio2.h']

    def fetch(self):
        pass
    

    def build(self):
        pass



class SExtractor(Package):
    """
    Attempt to find the sextractor and WCS libraries.
    Usually under the frogshome.
    """
    def __init__(self, frogshome=None):
        Package.__init__(self, 'sex', frogshome=frogshome)

        # get rid of this !
        if os.path.isdir(frogshome + os.sep + 'sex' + os.sep + 'v222'):
            self.prefixsearchpath = frogshome + os.sep + 'sex' + os.sep + 'v222'
            
        self.libpattern = ['libsex.so', 'libwcs.so']
        self.headerpattern = ['define.h']

        
    def fetch(self):
        pass


    def build(self):
        pass


class Boost(Package):
    """
    Boost is usually installed somewhere in the standard
    locations.
    """
    def __init__(self, frogshome=None):
        Package.__init__(self, 'boost', frogshome=frogshome)
        
        self.libpattern = ['libboost_serialization*']
        self.headerpattern = ['serialization/serialization.hpp']
        
        
    def fetch(self):
        pass
    
    
    def build(self):
        pass

    
    
def FindPackage(name):
    """
    Factory function. Return the requested package object.
    """
    pass









#    def locateLibs(self):
#        """
#        Try to locate the package libraries.
#        """
#        
#        if len(self.libs) == 0:
#            return False
#        
#        if len(self.prefixsearchpath) > 0:
#            path,libs = self._locateFilesUnderPrefix(self.libpattern, 0111,
#                                                     self.prefixsearchpath)
#            if path != None:
#                self.libpath = path
#                self.libs = libs
#                return True

#        if len(self.libsearchpath) > 0:
#            path,libs = self._locateFilesInPath(self.libpattern, 0111,
#                                                self.libsearchpath)
#            if path != None:
#                self.libpath = path
#                self.libs = libs
#                return True
#
#        return False
    
    
#    def locateHeaders(self):
#        """
#        Try to locate the package headers.
#        """
#        
#        if len(self.headers) == 0:
#            return
        
#        if len(self.prefixsearchpath) > 0:
#            path,headers = self._locateFilesUnderPrefix(self.headerpattern, 0444,
#                                                        self.prefixsearchpath)
#            if path != None:
#                self.headerpath = path
#                self.headers = headers
#                return True
#            
#        return False
    
    
#    def locatePrograms(self):
#        """
#        Try to locate the package programs (in the user PATH)...
#        """
#
#        if len(self.programs) == 0:
#            return False
#
#        if len(self.prefixsearchpath) > 0:
#            path,progs = self._locateFilesUnderPrefix(self.progpattern, 0111,
#                                                      self.prefixsearchpath)
#            if path != None:
#                self.progpath = path
#                self.programs = progs
#                return True
#            
#        return False




#    def locateHeaders(self):
#
#        headers = []
#        if len(self.prefixsearchpath) > 0:
#            headers = self.locateFilesUnderPrefix(self.headerpattern, 0444,
#                                                  self.prefixsearchpath)
#
#        if len(headers) == 0:
#            headers = self.locateFilesUnderPrefix(self.headerpattern, 0444,
#                                                  self.prefixsearchpath)
#
#        if len(headers) == 0:
#            return False
#        
#        for h in headers:
#            d = os.path.dirname(h)
#            if not d in self.includepath:
#                self.includepath.append(d)
            
            
#    def locateLibs(self):
#        
#        libs = []
#        if len(self.prefixsearchpath) > 0:
#            libs = self.locateFilesUnderPrefix(self.libpattern, 0111,
#                                               self.prefixsearchpath)
#            
#        if len(libs) == 0:
#            libs = self.locateFilesInPath(self.libpattern, 0111,
#                                          self.libsearchpath)
#            
#        if len(libs) == 0:
#            return False
#        
#        exp = re.compile('lib([\w\d-]+)\..*')
#        for l in libs:
#            self.libpath.append(os.path.dirname(l))
#            m = exp.search(os.path.basename(l))
#            if m == None:
#                return False
#            self.libs.append(m.group(1))
##
#
#      print self.libpath
#        print self.libs
#        
#        return True
