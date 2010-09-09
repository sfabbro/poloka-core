
## this library is not compiled if sextractor is disable ( --disable-sex in configure)



AM_CPPFLAGS = -Wall -I$(top_srcdir)/src_base -I$(top_srcdir)/src_utils  -I$(top_srcdir)/cern_stuff @SEX_CPPFLAGS@ @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@ 
AM_LDFLAGS = -L$(top_builddir)/src_base  -L$(top_builddir)/src_utils -L$(top_builddir)/cern_stuff  -L$(top_builddir)/lapack_stuff  -L@libdir@ -ltoadsbase -ltoadsutils -lcernstuff -llapackstuff @SEX_LDFLAGS@ @CFITSIO_LDFLAGS@ @CMT_LDFLAGS@

AM_YFLAGS = -l -d

if USESEXTRACTOR

#%_dict.cc : %.h
#	$(top_builddir)/objio/dictgen -o $@ $<

#point_dict.cc : $(top_srcdir)/src_base/point.h
#	$(top_builddir)/objio/dictgen -o $@ $<


AM_CPPFLAGS += $($<_cppflags)
optimize_it = -O3

cluster.cc_cppflags = $(optimize_it)
kernelfit.cc_cppflags = $(optimize_it)
dimage.cc_cppflags = $(optimize_it)
listmatch.cc_cppflags = $(optimize_it)
imagesubtraction.cc_cppflags =  $(optimize_it)
convolution.cc_cppflags =  $(optimize_it)
imagesum.cc_cppflags =  $(optimize_it)
fastfinder.cc_cppflags =  $(optimize_it)
gtransfo.cc_cppflags = $(optimize_it)
detection.cc_cppflags = $(optimize_it)
apersestar.cc_cppflags = $(optimize_it)
imageutils.cc_cppflags = $(optimize_it)


lib_LTLIBRARIES = libtoads.la
libtoads_la_SOURCES =  \
  wcsutils.cc                 wcsutils.h       \
  toadscards.cc               toadscards.h     \
  datacards.cc                datacards.h      \
  datdetec.cc                 datdetec.h       \
                              countedref.h     \
  basestar.cc                 basestar.h       \
  nstarmatch.cc               nstarmatch.h     \
  sestar.cc                   sestar.h         \
  starlist.cc                 starlist.h  starlistexception.h     \
  starmatch.cc                starmatch.h      \
  listmatch.cc                listmatch.h matchexception.h     \
  fastfinder.cc               fastfinder.h     \
  usnoutils.cc                usnoutils.h      \
  gtransfo.cc                 gtransfo.h       \
  align_box.cc                align_box.h      \
  quali_box.cc                quali_box.h      \
  seeing_box.cc               seeing_box.h     \
  sextractor_box.cc           sextractor_box.h \
  reducedutils.cc             reducedutils.h   \
  dimage.cc                   dimage.h           \
  detection.cc                detection.h        \
  allreducedimage.cc          allreducedimage.h  \
  dbimage.cc                  dbimage.h  dbconfigexception.h \
  dbfileParse.y                          \
  dbfileLex.l                            \
  imagelist.h        \
  imagematch.cc               imagematch.h       \
  imagesubtraction.cc         imagesubtraction.h \
  sub.cc                      sub.h              \
  imagesum.cc                 imagesum.h         \
  reducedimage.cc             reducedimage.h     \
  subimage.cc                 subimage.h         \
  transformedimage.cc         transformedimage.h \
  cluster.cc                  cluster.h  \
  convolution.cc              convolution.h \
                              convtemplate.h \
  kernelfit.cc                kernelfit.h iohelpers.h \
  kernelfitter.cc             kernelfitter.h  \
  imagepair.cc	              imagepair.h \
  scorescollection.cc         scorescollection.h \
  swarpstack.cc		      swarpstack.h \
  dicstar.cc                  dicstar.h \
  histopeakfinder.cc          histopeakfinder.h \
  dictfile.cc		      dictfile.h \
  apersestar.cc       apersestar.h \
  globalval.cc    globalval.h \
  hprofile1d.cc   hprofile1d.h \
  imageutils.cc imageutils.h \
  gaussianfit.cc gaussianfit.h\
                 ccd.h \
 fastifstream.cc fastifstream.h \
 lambert_dust.cc lambert_dust.h \
               lambert_tag.h \
  ellipticaper.cc  ellipticaper.h\
  multimagstar.h multimagstar.cc\
  pegase_utils.h pegase_utils.cc
#  imagesum_dict.cc \
#  imagesubtraction_dict.cc \
#  subimage_dict.cc \
#  swarpstack_dict.cc \
#  transformedimage_dict.cc 

include_HEADERS = \
dictfile.h \
dicstar.h \
wcsutils.h       \
toadscards.h     \
datacards.h      \
datdetec.h       \
countedref.h     \
basestar.h       \
fatpoint.h       \
nstarmatch.h     \
sestar.h         \
starlist.h       \
starmatch.h      \
listmatch.h      \
fastfinder.h     \
usnoutils.h      \
gtransfo.h       \
align_box.h      \
quali_box.h      \
seeing_box.h     \
sextractor_box.h \
reducedutils.h   \
dimage.h           \
detection.h        \
allreducedimage.h  \
dbimage.h  dbconfigexception.h        \
imagelist.h        \
imagematch.h       \
imagesubtraction.h \
sub.h              \
imagesum.h         \
reducedimage.h     \
subimage.h         \
transformedimage.h \
cluster.h  \
convolution.h exceptions.h\
convtemplate.h \
kernelfit.h \
swarpstack.h    \
nightelement.h  \
histopeakfinder.h \
apersestar.h \
ellipticaper.h \
pegase_utils.h \
multimagstar.h \
hprofile1d.h \
imageutils.h \
gaussianfit.h \
ccd.h \
lambert_dust.h \
fastifstream.h \
lambert_tag.h


wcsutils.o :: $(top_srcdir)/telinst/*.cc

noinst_HEADERS =  defs.h simps2d.h


CLEANFILES = *_dict.h *_dict.cc *__persister.h *__persister.cc

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

all : $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/libtoads.*))  $(patsubst %,@prefix@/include/%,${include_HEADERS})

endif

endif
