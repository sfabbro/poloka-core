EXTRA_DIST = snls.xsd

lib_LTLIBRARIES = libtoadio.la
libtoadio_la_SOURCES = xmlstream.cc dict.cc swig_dict_reader.cc xmlschema_dict_writer.cc typemgr.cc
libtoadio_la_CPPFLAGS = @XML_CFLAGS@
libtoadio_la_LDFLAGS = @XML_LIBS@


#%_dict.cc : %.h
#	./dictgen2 -o $@ $<

include_HEADERS = \
 persister.h dict.h objio.h dictio.h reader.h \
 xmlstream.h toadtypes.h persistence.h typemgr.h

noinst_HEADERS = \
 codegen.h

bin_PROGRAMS = dictgen2 test_swig test_dict_rw

dictgen2_SOURCES = dictgen2.cc codegen.cc codegen.h
dictgen2_CPPFLAGS = @XML_CFLAGS@
dictgen2_LDFLAGS = -ltoadio -lxml2 @XML_LIBS@

test_swig_SOURCES = test_swig_dict.cc test_swig_dict.h test_swig.cc 
test_swig_CPPFLAGS = @XML_CFLAGS@
test_swig_LDFLAGS = -ltoadio -lxml2 @XML_LIBS@

test_dict_rw_SOURCES = test_dict_rw.cc 
test_dict_rw_CPPFLAGS = @XML_CFLAGS@
test_dict_rw_LDFLAGS = -ltoadio -lxml2 @XML_LIBS@


CLEANFILES = *_dict.h *_dict.cc *__persister.h *__persister.cc
