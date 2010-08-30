CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in \
	stamp-* *.new

lib_LTLIBRARIES = libtoadsflat.la

include_HEADERS = \
	fitsimagearray.h \
	fitsset.h \
	fringeutils.h \
	superflat.h

libtoadsflat_la_SOURCES = \
	$(include_HEADERS) \
	fitsimagearray.cc \
	fitsset.cc \
	fringeutils.cc \
	superflat.cc

AM_CPPFLAGS = \
	-I$(top_srcdir)/src_base \
	-I$(top_srcdir)/src_utils \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src_base \
	 @CFITSIO_CFLAGS@

libtoadsflat_la_LIBADD = \
	$(top_builddir)/src_base/libtoadsbase.la \
	$(top_builddir)/src_utils/libtoadsutils.la

bin_PROGRAMS = correlation defringe elixir_to_toads flatfield fringefinder medianstack selection trim

LDADD = libtoadsflat.la

correlation_SOURCES = correlation.cc
defringe_SOURCES = defringe.cc
elixir_to_toads_SOURCES  = elixir_to_toads.cc
elixir_to_toads_LDADD = $(LDADD) $(top_builddir)/src/libtoads.la
flatfield_SOURCES  = flatfield.cc
flatfield_LDADD = $(LDADD) $(top_builddir)/src/libtoads.la
fringefinder_SOURCES = fringefinder.cc
medianstack_SOURCES = medianstack.cc
selection_SOURCES  = selection.cc
trim_SOURCES = trim.cc
