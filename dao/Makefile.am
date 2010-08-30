EXTRA_DIST = daocommon.f

AM_CPPFLAGS = -Wall -I$(top_srcdir)/src_base -I$(top_srcdir)/src @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@
AM_FFLAGS   = -Wall -pedantic -W  -O -funroll-loops -ffortran-bounds-check
#AM_FFLAGS   = -Wall -pedantic -W -ffortran-bounds-check
AM_LDFLAGS  = -L$(top_builddir)/src_base -L$(top_builddir)/src -ltoadsbase -ltoads -L@libdir@ -lg2c @CFITSIO_LDFLAGS@ @CMT_LDFLAGS@ 

if USESEXTRACTOR

lib_LTLIBRARIES = libtoadsdao.la

libtoadsdao_la_SOURCES = \
	addstr.f         \
	allstr.f	 \
	allstsubs.f	 \
	daotools.f	 \
	find.f		 \
	fotometry.f	 \
	fudge.f		 \
	getpsf.f	 \
	group.f		 \
	iosubs.f	 \
	mathsubs.f	 \
	nstar.f		 \
	pckpsf.f	 \
	peak.f		 \
	sort.f		 \
	substar.f	 \
	unxsubs.f	 \
	getdaocommon.f   \
	fitssubs.c	 \
	daophot.cc	 \
	daophotio.cc	 \
	daophotoption.cc \
	daophotpsf.cc	 \
	daophotutils.cc 

include_HEADERS = daophot.h cdaophot.h daophotio.h daophotoption.h daophotpsf.h daophotutils.h fitssubs.h


######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package

if AUTOINSTALL

@prefix@/lib/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* @prefix@/lib/$*

all :  $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/libtoadsdao.*))

endif

endif