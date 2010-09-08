EXTRA_DIST =

if IMAGEMAGICK

 bin_PROGRAMS = colorize makegif

 AM_CPPFLAGS = -I$(top_srcdir)/src_base  @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@ @IMAGEMAGICK_CPPFLAGS@
 AM_LDFLAGS = -L$(top_builddir)/src_base -L@libdir@ -ltoadsbase @CMT_LDFLAGS@  @CFITSIO_LDFLAGS@  @IMAGEMAGICK_LDFLAGS@ @IMAGEMAGICK_LIBS@ 

 colorize_SOURCES = colorize.cc
 makegif_SOURCES  = makegif.cc


if AUTOINSTALL

@prefix@/bin/% : %
	ln -fs @abs_builddir@/.libs/$* @prefix@/bin/$*

all   : $(patsubst %,@prefix@/bin/%,${bin_PROGRAMS})

endif


endif