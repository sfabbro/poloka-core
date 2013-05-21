CLEANFILES = *~
AINTAINERCLEANFILES = \
        Makefile.in \
        stamp-*

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_builddir) @CFITSIO_CFLAGS@

AM_DEFAULT_SOURCE_EXT = .cc

bin_PROGRAMS = \
	pka-fitshead pka-fitsmerge \
	pka-imselect pka-imls \
	pka-imarith pka-imstats pka-imsky \
	pka-imcosmic pka-imdeadpix pka-imtracks pka-imsatur \
	pka-imweight pka-imcatalog pka-imaper \
	pka-imcombine pka-imregister pka-imextract pka-imswarp \
	pka-listjoin pka-listmerge \
	pka-wcsconvert pka-wcsfit

LDADD = $(top_builddir)/poloka/libpoloka-core.la
