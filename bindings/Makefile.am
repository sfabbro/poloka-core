EXTRA_DIST = repack_python.swg

# FIXME: should be SWIG
SWIG=swig
# FIXME: should be PERL5_INCLUDES
# can be checked with perl -V:installarchlib
# PERL_INCLUDES=-I/usr/lib/perl/5.8.0/CORE -I../toadbase
PYTHON_INCLUDES=-I/usr/include/python2.2 -I../src
AM_CPPFLAGS=-I/usr/include/python2.2 -I../src
LDFLAGS=-L../src -ltoads -L../cern_stuff -lcernstuff -L/home/nrl/Work/frogs/sex/v222/Linux-i686 -lsex -lwcs -lg2c

lib_LTLIBRARIES = librepack_python.la
librepack_python_la_SOURCES = repack_python_wrap.cxx
#librepack_perl5_la_SOURCES = repack_perl5_wrap.cxx
#librepack_perl5_la_CPPFLAGS=$(PERL_INCLUDES)

repack_python_wrap.cxx: repack_python.swg
	$(SWIG) -c++ -interface librepack_python $(PYTHON_INCLUDES) -python $<

#repack_perl5_wrap.cxx: repack_perl5.swg
#	$(SWIG) -c++ $(PERL_INCLUDES) -I../toadbase -perl $<

CLEANFILES = *.cxx *.py *.pm *.pyc

