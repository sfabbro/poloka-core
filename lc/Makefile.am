CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in *.new stamp-*
lib_LTLIBRARIES = libtoadslc.la

lc_includedir = $(includedir)/poloka
lc_include_HEADERS = \
	fiducial.h \
	gausspsf.h \
	lcio.h \
	lightcurve.h \
	lightcurvepoint.h \
	photstar.h \
	refstar.h \
	simfit.h \
	simfitphot.h \
	simfitvignet.h \
	vignet.h \
	vignetphot.h \
	vignetserver.h

libtoadslc_la_SOURCES =  \
	$(lc_include_HEADERS) \
	gausspsf.cc \
	lcio.cc \
	lightcurve.cc \
	lightcurvepoint.cc \
	photstar.cc \
	refstar.cc \
	simfit.cc \
	simfitphot.cc \
	simfitvignet.cc \
	vignet.cc \
	vignetphot.cc \
	vignetserver.cc

AM_CPPFLAGS = \
	-I$(top_srcdir)/src_base \
	-I$(top_srcdir)/src_utils \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/psf \
	@CFITSIO_CFLAGS@

libtoadslc_la_LIBADD = \
	$(top_builddir)/src/libtoads.la \
	$(top_builddir)/psf/libtoadspsf.la
