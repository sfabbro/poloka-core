
AM_CPPFLAGS = -Wall -I../src -I../objio @XML_CFLAGS@ @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@
AM_FFLAGS   = -Wall -pedantic -W -O3 -funroll-loops
AM_LDFLAGS  = -L../src -ltoads -L../objio -L@libdir@ -ltoadio @XML_LIBS@ @CFITSIO_LDFLAGS@ @CMT_LDFLAGS@ 


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
	fitssubs.c	 \
	daophot.cc	 \
	daophotio.cc	 \
	daophotoption.cc \
	daophotpsf.cc	 \
	daophotutils.cc
