## this library is not compiled if sextractor is disable ( --disable-sex in configure)

if USESEXTRACTOR

AM_CPPFLAGS = -Wall -I$(top_srcdir)/src_base -I$(top_srcdir)/src_utils -I$(top_srcdir)/src -I$(top_srcdir)/dao -I$(top_srcdir)/psf -I$(top_builddir)/src_base -I$(top_builddir)/src_utils -I$(top_builddir)/src -I$(top_builddir)/dao  -I$(top_builddir)/psf @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@ 

#### -DFNAME

AM_LDFLAGS  = -L$(top_builddir)/src_base -L$(top_builddir)/src_utils -L$(top_builddir)/src -L$(top_builddir)/dao  -L$(top_builddir)/psf -L@libdir@ -ltoadsbase -ltoadsutils -ltoads CFITSIO_LDFLAGS@ @CMT_LDFLAGS@


lib_LTLIBRARIES = libtoadslc.la

AM_CPPFLAGS += $($<_cppflags)
optimize_it = -O3

vignet.cc_cppflags = $(optimize_it)
vignetphot.cc_cppflags = $(optimize_it)
simfitvignet.cc_cppflags = $(optimize_it)
simfit.cc_cppflags = $(optimize_it)
simfitphot.cc_cppflags = $(optimize_it)
lightcurve.cc_cppflags = $(optimize_it)

libtoadslc_la_SOURCES =  \
	lightcurvepoint.cc \
	gausspsf.cc	 \
	photstar.cc      \
	vignet.cc	 \
	vignetphot.cc	 \
	refstar.cc       \
	lcio.cc	         \
	lightcurve.cc	 \
	simfitvignet.cc  \
	simfit.cc        \
	simfitphot.cc 	 \
	vignetserver.cc
#	lightcurvepoint_dict.cc 

include_HEADERS = lightcurvepoint.h 

######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package

if AUTOINSTALL

@prefix@/lib/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* $@

@prefix@/include/% : @abs_srcdir@/%
	ln -fs @abs_srcdir@/$* $@

all : $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/libtoadslc.*)) $(patsubst %,@prefix@/include/%,${include_HEADERS})

endif

endif
