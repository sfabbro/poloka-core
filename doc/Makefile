# Generated automatically from Makefile.in by configure.
# Makefile.in generated automatically by automake 1.4-p5 from Makefile.am

# Copyright (C) 1994, 1995-8, 1999, 2001 Free Software Foundation, Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.


SHELL = /bin/sh

srcdir = .
top_srcdir = ..
prefix = /tmp/toto
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var
libdir = ${exec_prefix}/lib
infodir = ${prefix}/info
mandir = ${prefix}/man
includedir = ${prefix}/include
oldincludedir = /usr/include

DESTDIR =

pkgdatadir = $(datadir)/toads
pkglibdir = $(libdir)/toads
pkgincludedir = $(includedir)/toads

top_builddir = ..

ACLOCAL = aclocal
AUTOCONF = autoconf
AUTOMAKE = automake
AUTOHEADER = autoheader

INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL} $(AM_INSTALL_PROGRAM_FLAGS)
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_SCRIPT = ${INSTALL_PROGRAM}
transform = s,x,x,

NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
host_alias = i686-pc-linux-gnu
host_triplet = i686-pc-linux-gnu
AS = @AS@
AWK = gawk
CC = gcc
CERN_LDFLAGS = -L/coda/software/cern/pro/lib -lpacklib -lmathlib -lnsl -lcrypt -ldl -lg2c
CPP = gcc -E
CXX = c++
DLLTOOL = @DLLTOOL@
DOXYGEN = /coda/software/bin/i386_linux/doxygen
ECHO = echo
EXEEXT = 
F77 = g77
LEX = flex
LIBCFITSIO_CPPFLAGS = -I/usr/local//include
LIBCFITSIO_LDFLAGS = -L/usr/local//lib -lcfitsio
LIBTOOL = $(SHELL) $(top_builddir)/libtool
LIBXML2_CPPFLAGS = -I/usr/local//include/libxml2
LIBXML2_LDFLAGS = -L/usr/local//lib -lxml2
LN_S = ln -s
MAKEINFO = makeinfo
OBJDUMP = @OBJDUMP@
OBJEXT = o
PACKAGE = toads
RANLIB = ranlib
SEX_CPPFLAGS = -I/home/nrl/Work/frogs/sex/v222/src
SEX_LDFLAGS = -L/home/nrl/Work/frogs/sex/v222/Linux-i686 -lsex -lwcs
STRIP = strip
VERSION = 0.1.0
YACC = bison -y

EXTRA_DIST = README.db refman.ps Doxyfile.in

docdir = $(datadir)/doc
doc_DATA = html latex man
DOX = /coda/software/bin/i386_linux/doxygen
HEADER_LIST = ../src/*.h
mkinstalldirs = $(SHELL) $(top_srcdir)/config/mkinstalldirs
CONFIG_HEADER = ../config.h
CONFIG_CLEAN_FILES =  Doxyfile
DATA =  $(doc_DATA)

DIST_COMMON =  Doxyfile.in Makefile.am Makefile.in


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = gtar
GZIP_ENV = --best
all: all-redirect
.SUFFIXES:
$(srcdir)/Makefile.in: Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu doc/Makefile

Makefile: $(srcdir)/Makefile.in  $(top_builddir)/config.status $(BUILT_SOURCES)
	cd $(top_builddir) \
	  && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status

Doxyfile: $(top_builddir)/config.status Doxyfile.in
	cd $(top_builddir) && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status
tags: TAGS
TAGS:


distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)

subdir = doc

distdir: $(DISTFILES)
	here=`cd $(top_builddir) && pwd`; \
	top_distdir=`cd $(top_distdir) && pwd`; \
	distdir=`cd $(distdir) && pwd`; \
	cd $(top_srcdir) \
	  && $(AUTOMAKE) --include-deps --build-dir=$$here --srcdir-name=$(top_srcdir) --output-dir=$$top_distdir --gnu doc/Makefile
	@for file in $(DISTFILES); do \
	  d=$(srcdir); \
	  if test -d $$d/$$file; then \
	    cp -pr $$d/$$file $(distdir)/$$file; \
	  else \
	    test -f $(distdir)/$$file \
	    || ln $$d/$$file $(distdir)/$$file 2> /dev/null \
	    || cp -p $$d/$$file $(distdir)/$$file || :; \
	  fi; \
	done
info-am:
info: info-am
dvi-am:
dvi: dvi-am
check-am: all-am
check: check-am
installcheck-am:
installcheck: installcheck-am
install-exec-am:
install-exec: install-exec-am

install-data-am: install-docDATA
install-data: install-data-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am
install: install-am
uninstall-am: uninstall-docDATA
uninstall: uninstall-am
all-am: Makefile $(DATA)
all-redirect: all-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) AM_INSTALL_PROGRAM_FLAGS=-s install
installdirs:
	$(mkinstalldirs)  $(DESTDIR)$(docdir)


mostlyclean-generic:

clean-generic:

distclean-generic:
	-rm -f Makefile $(CONFIG_CLEAN_FILES)
	-rm -f config.cache config.log stamp-h stamp-h[0-9]*

maintainer-clean-generic:
mostlyclean-am:  mostlyclean-generic

mostlyclean: mostlyclean-am

clean-am:  clean-generic mostlyclean-am clean-local

clean: clean-am

distclean-am:  distclean-generic clean-am
	-rm -f libtool

distclean: distclean-am

maintainer-clean-am:  maintainer-clean-generic distclean-am
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."

maintainer-clean: maintainer-clean-am

.PHONY: uninstall-docDATA install-docDATA tags distdir info-am info \
dvi-am dvi check check-am installcheck-am installcheck install-exec-am \
install-exec install-data-am install-data install-am install \
uninstall-am uninstall all-redirect all-am all installdirs \
mostlyclean-generic distclean-generic clean-generic \
maintainer-clean-generic clean mostlyclean distclean maintainer-clean


all: stamp-doc

html: stamp-doc

stamp-doc: Doxyfile $(HEADER_LIST)
	$(DOX) $<
	@echo > timestamp > $@

install-docDATA:
	$(mkinstalldirs) $(docdir);
	(cp -r html $(docdir); cd $(docdir)/html; rm -rf CVS)
	(cp -r man $(docdir);cd $(docdir)/man; rm -rf CVS)

uninstall-docDATA:
	cd $(docdir) ; rm -rf html
	cd $(docdir) ; rm -rf man
	rmdir $(docdir)

clean-local:
	@rm -rf html latex man rtf stamp-doc
	@mkdir html latex man rtf

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
