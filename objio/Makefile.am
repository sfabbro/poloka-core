EXTRA_DIST = snls.xsd


lib_LTLIBRARIES = libtoadio.la
libtoadio_la_SOURCES = xmlstream.cc dict.cc swig_dict_reader.cc xmlschema_dict_writer.cc
libtoadio_la_CPPFLAGS = -I../toadbase -I/usr/include/libxml2 @LIBXML2_CPPFLAGS@
libtoadio_la_LDFLAGS = -lxml2 @LIBXML2_LDFLAGS@


include_HEADERS = \
 persister.h dict.h objio.h dictio.h reader.h \
 xmlstream.h

noinst_HEADERS = \
 codegen.h

bin_PROGRAMS = dictgen2 test_swig test_dict_rw
#dictgen_SOURCES = dictgen.cc
#dictgen_CPPFLAGS = -I/usr/include/libxml2 @LIBXML2_CPPFLAGS@
#dictgen_LDFLAGS = -lxml2 @LIBXML2_LDFLAGS@

dictgen2_SOURCES = dictgen2.cc codegen.cc codegen.h
dictgen2_CPPFLAGS = -I../toadbase -I/usr/include/libxml2 @LIBXML2_CPPFLAGS@
dictgen2_LDFLAGS = -ltoadio -lxml2 @LIBXML2_LDFLAGS@

test_swig_SOURCES = test_swig.cc test_swig_dict.cc test_swig_dict.h
test_swig_CPPFLAGS = -I../toadbase -I/usr/include/libxml2 @LIBXML2_CPPFLAGS@
test_swig_LDFLAGS = -ltoadio -lxml2 @LIBXML2_LDFLAGS@

test_dict_rw_SOURCES = test_dict_rw.cc 
test_dict_rw_CPPFLAGS = -I../toadbase -I/usr/include/libxml2 @LIBXML2_CPPFLAGS@
test_dict_rw_LDFLAGS = -ltoadio -lxml2 @LIBXML2_LDFLAGS@

%_dict.cc : %.h
	./dictgen2 -o $@ $<

