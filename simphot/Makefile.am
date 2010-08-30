CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in \
	stamp-* *.new

EXTRA_DIST = README

lib_LTLIBRARIES = libtoadssimphot.la

include_HEADERS = \
	array2d.h \
	array4d.h \
	calibratedstar.h \
	imagepsfserver.h \
	intframe.h \
	intpoint.h \
	lightcurvefile.h \
	lightcurvesyntaxerror.h \
	model.h \
	objecttofit.h \
	pixelblock.h \
	resampler.h \
	simphotfit.h \
	vignette.h

libtoadssimphot_la_SOURCES = \
	$(include_HEADERS) \
	array4d.cc \
	calibratedstar.cc \
	imagepsfserver.cc \
	intframe.cc \
	lightcurvefile.cc \
	model.cc \
	objecttofit.cc \
	pixelblock.cc \
	resampler.cc \
	simphotfit.cc \
	vignette.cc

# MATVECT_CHECK_BOUNDS checks access to Mat and Vect for index over/underflow
AM_CPPFLAGS = \
	-I$(top_srcdir)/src_base \
	-I$(top_srcdir)/src_utils \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/psf \
	@CFITSIO_CFLAGS@ -DMATVECT_CHECK_BOUNDS

libtoadssimphot_la_LIBADD = \
	$(top_builddir)/src/libtoads.la \
	$(top_builddir)/psf/libtoadspsf.la

LDADD = libtoadssimphot.la

check_PROGRAMS =  testlc testpsf
bin_PROGRAMS = mklc

mklc_SOURCES = mklc.cc
testlc_SOURCES = testlc.cc
testpsf_SOURCES = testpsf.cc
