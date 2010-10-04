
AM_CPPFLAGS = -Wall -I$(top_srcdir)/src_base  -I$(top_srcdir)/cern_stuff  @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@ 
AM_LDFLAGS = -L$(top_builddir)/src_base -L$(top_builddir)/cern_stuff -L$(top_builddir)/lapack_stuff -L@libdir@ -ltoadsbase -llapack -lcernstuff @CFITSIO_LDFLAGS@ @CMT_LDFLAGS@


AM_CPPFLAGS += $($<_cppflags)
optimize_it = -O3

histo1d.cc_cppflags =    $(optimize_it)
histo2d.cc_cppflags =    $(optimize_it)
histo4d.cc_cppflags =    $(optimize_it)
matvect.cc_cppflags =    $(optimize_it)
imageinterpolation.cc_cppflags =    $(optimize_it)
imageback.cc_cppflags = $(optimize_it)

lib_LTLIBRARIES = libtoadsutils.la
libtoadsutils_la_SOURCES = \
	imageback.cc		      imageback.h \
	matvect.cc		      matvect.h \
	histo1d.cc		      histo1d.h \
	histo2d.cc		      histo2d.h \
	histo4d.cc		      histo4d.h \
	imageinterpolation.cc	      imageinterpolation.h \
	fitsslice.cc		      fitsslice.h
include_HEADERS = \
	matvect.h \
	imageback.h \
	histo1d.h \
	histo2d.h \
	histo4d.h \
	imageinterpolation.h \
	fitsslice.h

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

all : $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/libtoadsutils.*))  $(patsubst %,@prefix@/include/%,${include_HEADERS})

endif
