EXTRA_DIST = repack_python.swg


PYTHON_INCLUDES=-I/usr/include/python@PYTHON_VERSION@ -I../src
AM_CPPFLAGS=-I/usr/include/python@PYTHON_VERSION@ -I../src
AM_LDFLAGS=-L../src -ltoads -L../cern_stuff -lcernstuff @SEX_LDFLAGS@ -lg2c

lib_LTLIBRARIES = librepack_python.la
librepack_python_la_SOURCES = repack_python_wrap.cxx
#librepack_perl5_la_SOURCES = repack_perl5_wrap.cxx
#librepack_perl5_la_CPPFLAGS=$(PERL_INCLUDES)

repack_python_wrap.cxx: repack_python.swg
	$(SWIG) -interface librepack_python $(PYTHON_INCLUDES) -python $<

#repack_perl5_wrap.cxx: repack_perl5.swg
#	$(SWIG) $(PERL_INCLUDES) -I../toadbase -perl $<

CLEANFILES = *.cxx *.py *.pm *.pyc


# FIXME: should be PERL5_INCLUDES
# can be checked with perl -V:installarchlib
# PERL_INCLUDES=-I/usr/lib/perl/5.8.0/CORE -I../toadbase

