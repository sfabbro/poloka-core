
AM_CPPFLAGS = -Wall -I../src -I../objio -I../dao @XML_CFLAGS@ @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@ -DFNAME
AM_LDFLAGS  = -L../src -L../dao -L../objio -L@libdir@ -ltoadio @XML_LIBS@ CFITSIO_LDFLAGS@ @CMT_LDFLAGS@

%_dict.cc : %.h
	../objio/dictgen -o $@ $<

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
	lightcurvepoint_dict.cc \
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
	simfitphot.cc 	 

include_HEADERS = lightcurvepoint.h

######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package

if AUTOINSTALL

@prefix@/lib/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* @prefix@/lib/$*

@prefix@/include/% : @abs_builddir@/%
	ln -fs @abs_builddir@/$* @prefix@/include/$*

all : $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/libtoadslc.*)) $(patsubst %,@prefix@/include/%,${include_HEADERS})

endif
