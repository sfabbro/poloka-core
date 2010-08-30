CLEANFILES = *~
MAINTAINERCLEANFILES = \
	aclocal.m4 \
	configure \
	Makefile.in \
	config.h \
	config.h.in \
	autobuild/* \
	stamp-* \
	poloka.pc

EXTRA_DIST = autogen.sh README NEWS ChangeLog INSTALL COPYING

ACLOCAL_AMFLAGS = -I m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA= poloka.pc

SUBDIRS = datacards doc \
	cern_stuff telinst \
	src_base src_utils src \
	flat psf lc simphot utils

if CERNLIB
SUBDIRS += cern_utils
endif

if MC
SUBDIRS += mc
endif

if IMAGEMAGICK
SUBDIRS += imagemagick_utils
endif
