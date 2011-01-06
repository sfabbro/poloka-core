
# MATVECT_CHECK_BOUNDS checks access to Mat and Vect for index over/underflow



libtag=simphot

AM_CPPFLAGS = -Wall -I$(top_srcdir)/src_base -I$(top_srcdir)/src -I$(top_srcdir)/simphot -I$(top_srcdir)/psf -I$(top_srcdir)/src_utils @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@ -DMATVECT_CHECK_BOUNDS
#AM_FFLAGS   = -Wall -pedantic -W  -O3 -funroll-loops

AM_LDFLAGS  = -L$(top_builddir)/src_base -L$(top_builddir)/src -ltoadsbase -ltoads   -L$(top_builddir)/simphot  -R @libdir@ -L@libdir@ -lgfortran @CFITSIO_LDFLAGS@ @CMT_LDFLAGS@ 



AM_CPPFLAGS += $($<_cppflags)
optimize_it = -O3



lib_LTLIBRARIES = libtoadssimphot.la

libtoadssimphot_la_SOURCES = \
	lightcurvefile.cc lightcurvefile.h  lightcurvesyntaxerror.h \
	objecttofit.cc objecttofit.h \
	simphotfit.cc simphotfit.h \
	vignette.cc vignette.h \
	model.cc model.h \
	array4d.cc array2d.h array4d.h \
	pixelblock.cc pixelblock.h \
	resampler.cc resampler.h \
	intframe.cc intpoint.h intframe.h \
	imagepsfserver.cc imagepsfserver.h \
	calibratedstar.cc calibratedstar.h
include_HEADERS = lightcurvefile.h


resampler.cc_cppflags = $(optimize_it)
vignette.cc_cppflags = $(optimize_it)
array4d.cc_cppflags = $(optimize_it)

bin_PROGRAMS =  testlc testpsf mklc mklc_test

mklc_SOURCES = mklc.cc
mklc_LDFLAGS = -ltoadssimphot -R @libdir@ -L@libdir@ -ltoadspsf
mklc_test_SOURCES = mklc_test.cc
mklc_test_LDFLAGS = -ltoadssimphot -R @libdir@ -L@libdir@ -ltoadspsf
testlc_SOURCES = testlc.cc
testlc_LDFLAGS = -ltoadssimphot -ltoadspsf
testpsf_SOURCES = testpsf.cc
testpsf_LDFLAGS = -ltoadssimphot -ltoadspsf

######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package

if AUTOINSTALL

@prefix@/bin/% : %
	ln -fs @abs_builddir@/.libs/$* $@

all   : $(patsubst %,@prefix@/bin/%,${bin_PROGRAMS})

@prefix@/lib/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* $@

all :  $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/libtoadssimphot.*))

endif

