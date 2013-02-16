CLEANFILES = *~
AINTAINERCLEANFILES = \
        Makefile.in \
        stamp-*

AM_CPPFLAGS = @CFITSIO_CFLAGS@

AM_DEFAULT_SOURCE_EXT = .cc

bin_PROGRAMS = \
	pka-fitskey pka-fitshead pka-fitsmerge \
	pka-fitsselect pka-fitsstats \
	pka-imarith pka-imls pka-imstats pka-imsky \
	pka-imcosmic pka-imdeadpix pka-imtracks pka-imsatur \
	pka-imweight pka-imcatalog pka-imaper \
	pka-imcombine pka-imregister pka-imextract pka-imwarp \
	pka-listjoin pka-listmerge \
	pka-wcsconvert pka-wcsfit

LDADD = $(top_builddir)/poloka/libpoloka-core.la
