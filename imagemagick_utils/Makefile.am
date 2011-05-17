CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in \
	stamp-*

AM_CPPFLAGS = \
	-I$(top_srcdir)/src_base \
	 @CFITSIO_CFLAGS@ @IMAGEMAGICK_CFLAGS@
LDADD = \
	$(top_builddir)/src_base/libtoadsbase.la \
	@IMAGEMAGICK_LIBS@

bin_PROGRAMS = colorize makegif

colorize_SOURCES = colorize.cc
makegif_SOURCES  = makegif.cc
