EXTRA_DIST =

bin_PROGRAMS = l2tup 2im2tup im2tup 
#mkskyhist sky_hist l2tup


AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/dbutils -I$(top_srcdir)/objio @XML_CFLAGS@ @SEX_CPPFLAGS@ @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@
AM_LDFLAGS = -L$(top_builddir)/src  -L$(top_builddir)/objio -L$(top_builddir)/cern_stuff  -L@libdir@ -ltoads -ltoadio   -lcernstuff  @CMT_LDFLAGS@ @CERNLIB_LDFLAGS@ @XML_LIBS@ @SEX_LDFLAGS@  @CFITSIO_LDFLAGS@ -lg2c


im2tup_SOURCES = im2tup.cc
l2tup_SOURCES = l2tup.cc
2im2tup_SOURCES = 2im2tup.cc
#mkskyhist_SOURCES = mkskyhist.cc
#sky_hist_SOURCES = sky_hist.cc
