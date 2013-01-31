CLEANFILES = *~
AINTAINERCLEANFILES = \
        Makefile.in \
        stamp-*

AM_CPPFLAGS = \
        -I$(top_srcdir)/src \
        @CFITSIO_CFLAGS@

AM_DEFAULT_SOURCE_EXT = .cc
bin_PROGRAMS = \
	dbls dbimstat \
	astrom matchusno \
	changekey header skylev saturlev \
	image_combine image_match image_select image_stats \
	image_extract image_zp image_compute \
	joinlist merge_headers mergelists \
	make_catalog make_cosmic \
	make_dead make_satur make_weight \
	findsatellite mkcat2

LDADD = $(top_builddir)/src/libpoloka-core.la
