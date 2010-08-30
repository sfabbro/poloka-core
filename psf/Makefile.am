CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in *.new stamp-*

EXTRA_DIST = README.ls_biases README_HOWTO

include_HEADERS = \
	analyticpsf.h \
	imagepsf.h \
	nonlinmodel.h \
	poly1.h \
	poly2.h \
	psfstar.h \
	sparsevect.h

lib_LTLIBRARIES = libtoadspsf.la

libtoadspsf_la_SOURCES = \
	$(include_HEADERS) \
	analyticpsf.cc \
	imagepsf.cc \
	nonlinmodel.cc \
	poly2.cc \
	poly1.cc \
	psfstar.cc

# MATVECT_CHECK_BOUNDS checks access to Mat and Vect for index over/underflow
AM_CPPFLAGS = \
	-I$(top_srcdir)/src_base \
	-I$(top_srcdir)/src_utils \
	-I$(top_srcdir)/src \
	@CFITSIO_CFLAGS@ -DMATVECT_CHECK_BOUNDS

libtoadspsf_la_LIBADD = $(top_builddir)/src/libtoads.la

LDADD = libtoadspsf.la

bin_PROGRAMS = makepsf
makepsf_SOURCES = makepsf.cc

