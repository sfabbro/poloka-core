EXTRA_DIST =

bin_PROGRAMS = fringefinder defringe medianstack trim correlation

AM_CPPFLAGS = -Wall -I$(top_srcdir)/src_base -I$(top_srcdir)/src_utils -I$(top_srcdir)/cern_stuff  @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@ 
AM_LDFLAGS = -L$(top_builddir)/src_base -L$(top_builddir)/src_utils  -L$(top_builddir)/cern_stuff -L$(top_builddir)/lapack_stuff -L$(top_builddir)/flat -R @libdir@ -L@libdir@   -ltoadsutils -llapack -lcernstuff -ltoadsbase  @CFITSIO_LDFLAGS@ @CMT_LDFLAGS@


AM_CPPFLAGS += $($<_cppflags)
optimize_it = -O3


fringeutils.cc_cppflags =    $(optimize_it)



defringe_LDFLAGS      =  -ltoadsflat
fringefinder_LDFLAGS  =  -ltoadsflat
medianstack_LDFLAGS   =  -ltoadsflat
trim_LDFLAGS          =  -ltoadsflat
correlation_LDFLAGS   =  -ltoadsflat

defringe_SOURCES = defringe.cc
fringefinder_SOURCES = fringefinder.cc
medianstack_SOURCES = medianstack.cc
trim_SOURCES = trim.cc
correlation_SOURCES = correlation.cc


lib_LTLIBRARIES = libtoadsflat.la
libtoadsflat_la_SOURCES = \
  fringeutils.cc              fringeutils.h    \
  fitsimagearray.cc           fitsimagearray.h \
  fitsset.cc                  fitsset.h \
  superflat.cc		      superflat.h



include_HEADERS = \
fringeutils.h    \
fitsimagearray.h \
fitsset.h \
superflat.h



######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package

if AUTOINSTALL

@prefix@/bin/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* $@

@prefix@/lib/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* $@

@prefix@/include/% : @abs_srcdir@/%
	ln -fs @abs_srcdir@/$* $@

all : $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/libtoadsflat.*))  $(patsubst %,@prefix@/include/%,${include_HEADERS})
all : $(patsubst %,@prefix@/bin/%,${bin_PROGRAMS})
endif
