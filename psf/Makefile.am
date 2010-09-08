
# MATVECT_CHECK_BOUNDS checks access to Mat and Vect for index over/underflow

AM_CPPFLAGS = -Wall -I$(top_srcdir)/src_base -I$(top_srcdir)/src -I$(top_srcdir)/psf -I$(top_srcdir)/src_utils @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@ -DMATVECT_CHECK_BOUNDS
#AM_FFLAGS   = -Wall -pedantic -W  -O3 -funroll-loops
AM_FFLAGS   = -Wall -pedantic -W -ffortran-bounds-check
AM_LDFLAGS  = -L$(top_builddir)/src_base -L$(top_builddir)/src -ltoadsbase -ltoads -L@libdir@ -lg2c @CFITSIO_LDFLAGS@ @CMT_LDFLAGS@

AM_CPPFLAGS += $($<_cppflags)
optimize_it = -O3


lib_LTLIBRARIES = libtoadspsf.la

libtoadspsf_la_SOURCES = \
	analyticpsf.cc   analyticpsf.h nonlinmodel.h \
	psfstar.cc       \
	imagepsf.cc      \
	poly2.cc         \
	poly1.cc        \
	nonlinmodel.cc  
include_HEADERS = poly1.h imagepsf.h psfstar.h poly2.h sparsevect.h analyticpsf.h



analyticpsf.cc_cppflags = $(optimize_it)
psfstar.cc_cppflags = $(optimize_it)
imagepsf.cc_cppflags = $(optimize_it)
poly2.cc_cppflags = $(optimize_it)


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

all :  $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/libtoadspsf.*)) $(patsubst %,@prefix@/include/%,${include_HEADERS})

endif

