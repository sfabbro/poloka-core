EXTRA_DIST = snls.xsd

AM_CPPFLAGS = -I../src @XML_CFLAGS@

lib_LTLIBRARIES = libtoadio.la
libtoadio_la_SOURCES = \
 xmlstream.cc dict.cc swig_dict_reader.cc xmlschema_dict_writer.cc \
 typemgr.cc cpptype.cc cppclassmember.cc cppclass.cc cppswigclassreader.cc
libtoadio_la_LDFLAGS = @XML_LIBS@


%_dict.cc : %.h
	./dictgen -o $@ $<

include_HEADERS = \
 persister.h dict.h objio.h dictio.h reader.h \
 xmlstream.h toadtypes.h persistence.h typemgr.h \
 cpptype.h cppclassmember.h cppclass.h

noinst_HEADERS = \
 codegen.h

bin_PROGRAMS = dictgen dictgen2 test_swig test_dict_rw test_cppstuff

test_cppstuff_SOURCES = test_cppstuff.cc 
test_cppstuff_CPPFLAGS = -I../src @XML_CFLAGS@
test_cppstuff_LDFLAGS = -ltoadio -lxml2 @XML_LIBS@

dictgen_SOURCES = dictgen.cc codegen.cc codegen.h
dictgen_CPPFLAGS = -I../src @XML_CFLAGS@
dictgen_LDFLAGS = -ltoadio -lxml2 @XML_LIBS@

dictgen2_SOURCES = dictgen2.cc codegen.cc codegen.h
dictgen2_CPPFLAGS = -I../src @XML_CFLAGS@
dictgen2_LDFLAGS = -ltoadio -lxml2 @XML_LIBS@

test_swig_SOURCES = test_swig_dict.cc test_swig_dict.h test_swig.cc 
test_swig_CPPFLAGS = -I../src @XML_CFLAGS@
test_swig_LDFLAGS = -ltoadio -lxml2 @XML_LIBS@

test_dict_rw_SOURCES = test_dict_rw.cc 
test_dict_rw_CPPFLAGS = -I../src @XML_CFLAGS@
test_dict_rw_LDFLAGS = -ltoadio -lxml2 @XML_LIBS@

CLEANFILES = *_dict.h *_dict.cc *__persister.h *__persister.cc
