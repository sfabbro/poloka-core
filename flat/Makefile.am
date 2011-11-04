CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in \
	stamp-* *.new

lib_LTLIBRARIES = libtoadsflat.la

flat_includedir = $(includedir)/poloka
flat_include_HEADERS = \
	fitsimagearray.h \
	fitsset.h \
	fringeutils.h \
	superflat.h

libtoadsflat_la_SOURCES = \
	$(flat_include_HEADERS) \
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

LDADD = libtoadsflat.la \
	$(top_builddir)/src_base/libtoadsbase.la \
	$(top_builddir)/src_utils/libtoadsutils.la

bin_PROGRAMS = correlation defringe fringefinder trim

correlation_SOURCES = correlation.cc
defringe_SOURCES = defringe.cc
fringefinder_SOURCES = fringefinder.cc
trim_SOURCES = trim.cc
