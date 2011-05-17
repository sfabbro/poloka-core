CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in stamp-*

bin_PROGRAMS = l2tup im2tup 2im2tup

AM_CPPFLAGS = -I$(top_srcdir)/src_base

LDADD = \
	$(top_builddir)/src_base/libtoadsbase.la \
	@CERNLIB_PACKLIBS@

im2tup_SOURCES = im2tup.cc
2im2tup_SOURCES = 2im2tup.cc
l2tup_SOURCES = l2tup.cc
