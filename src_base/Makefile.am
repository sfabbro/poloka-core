EXTRA_DIST = slasub.c alltelinst.cc alltelinst.h wcs_defines.h virtualinstrument.cc virtualinstrument.h


AM_CPPFLAGS = -Wall @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@ 
AM_LDFLAGS =  @CFITSIO_LDFLAGS@ @CMT_LDFLAGS@ -lz

AM_YFLAGS = -l -d

#%_dict.cc : %.h
#	$(top_builddir)/objio/dictgen -o $@ $<


AM_CPPFLAGS += $($<_cppflags)
optimize_it = -O3

alltelinst.cc alltelinst.h : $(top_srcdir)/telinst/*.cc $(top_srcdir)/telinst $(src)fitstoad.cc
	@ awk '{if (switch) print $$0}/FOR_MAKE/{switch=1}/^#endif/{if (switch) switch=0}' fitstoad.cc | $(MAKE) -f - `basename $@`


fitsimage.cc : alltelinst.h alltelinst.cc

vutils.cc_cppflags = $(optimize_it)
image.cc_cppflags  = $(optimize_it)
fitsimage.cc_cppflags = $(optimize_it)

lib_LTLIBRARIES = libtoadsbase.la
libtoadsbase_la_SOURCES = \
  image.cc                    image.h \
  fitsimage.cc                fitsimage.h \
  vutils.cc                   vutils.h \
  frame.cc                    frame.h \
  fileutils.cc                fileutils.h \
  stringlist.cc               stringlist.h \
  fitstoad.cc                 fitstoad.h \
  astroutils.cc               astroutils.h \
  wcscon.c                    wcscon.h \
  gzstream.cc                 gzstream.h \
  signalhandlers.cc fitsexception.h \
  trapfpe.c 
include_HEADERS = \
image.h \
fitsimage.h \
point.h \
vutils.h \
image.h \
frame.h \
fileutils.h \
stringlist.h \
fitstoad.h \
astroutils.h \
gzstream.h \
polokaexception.h \
alltelinst.h


noinst_HEADERS = wcscon.h


######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package

if AUTOINSTALL

@prefix@/bin/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* $@

@prefix@/lib/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* $@

@prefix@/include/% : @abs_srcdir@/%
	ln -fs @abs_srcdir@/$* $@

all : $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/libtoadsbase.*))  $(patsubst %,@prefix@/include/%,${include_HEADERS})

endif
