CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in \
	*.new stamp-*

src_utils_includedir = $(includedir)/poloka
src_utils_include_HEADERS = \
	fitsslice.h \
	histo1d.h \
	histo2d.h \
	histo4d.h \
	imageback.h \
	imageinterpolation.h \
	matvect.h

lib_LTLIBRARIES = libtoadsutils.la
libtoadsutils_la_SOURCES = \
	$(src_utils_include_HEADERS) \
	fitsslice.cc \
	histo1d.cc \
	histo2d.cc \
	histo4d.cc \
	imageback.cc \
	imageinterpolation.cc \
	matvect.cc

libtoadsutils_la_CPPFLAGS = \
	-I$(top_srcdir)/src_base \
	-I$(top_srcdir)/cern_stuff \
	@CFITSIO_CFLAGS@

libtoadsutils_la_LIBADD = \
	$(top_builddir)/src_base/libtoadsbase.la \
	$(top_builddir)/cern_stuff/libcernstuff.la \
	@LAPACK_LIBS@
