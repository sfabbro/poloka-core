EXTRA_DIST = slasub.c alltelinst.cc alltelinst.h wcs_defines.h


AM_CPPFLAGS = -Wall -I$(top_srcdir)/objio -I$(top_srcdir)/cern_stuff @XML_CFLAGS@ @SEX_CPPFLAGS@ @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@ 
AM_LDFLAGS = -L$(top_srcdir)/objio -L$(top_srcdir)/cern_stuff -L@libdir@ -ltoadio -lcernstuff @XML_LIBS@ @SEX_LDFLAGS@ @CFITSIO_LDFLAGS@ @CMT_LDFLAGS@

AM_YFLAGS = -l -d

%_dict.cc : %.h
	../objio/dictgen -o $@ $<


AM_CPPFLAGS += $($<_cppflags)
optimize_it = -O3



kernelfit.cc_cppflags = $(optimize_it)
image.cc_cppflags = $(optimize_it)
dimage.cc_cppflags = $(optimize_it)
listmatch.cc_cppflags = $(optimize_it)
imagesubtraction.cc_cppflags =  $(optimize_it)
convolution.cc_cppflags =  $(optimize_it)
vutils.cc_cppflags =      $(optimize_it)
imagesum.cc_cppflags =  $(optimize_it)
fastfinder.cc_cppflags =  $(optimize_it)
fringeutils.cc_cppflags =    $(optimize_it)
histo2d.cc_cppflags = $(optimize_it)
histo4d.cc_cppflags = $(optimize_it)
gtransfo.cc_cppflags = $(optimize_it)
detection.cc_cppflags = $(optimize_it)
matvec.cc_cppflags = $(optimize_it)

lib_LTLIBRARIES = libtoads.la
libtoads_la_SOURCES = \
  dbimage_dict.cc \
  frame_dict.cc \
  reducedimage_dict.cc \
  gtransfo_dict.cc \
  transformedimage_dict.cc \
  kernelfit_dict.cc \
  psfmatch_dict.cc \
  imagesubtraction_dict.cc \
  point_dict.cc \
  starlist_dict.cc \
  basestar_dict.cc \
  sestar_dict.cc \
  jimstar_dict.cc \
  standardstar_dict.cc \
  starmatch_dict.cc \
  detection_dict.cc \
 \
  matvect.cc                  matvect.h \
  vutils.cc                   vutils.h         \
  wcscon.c                    wcscon.h         \
  wcsutils.cc                 wcsutils.h       \
  toadscards.cc               toadscards.h     \
  datacards.cc                datacards.h      \
  datdetec.cc                 datdetec.h       \
  fileutils.cc                fileutils.h      \
  rootstuff.cc                rootstuff.h      \
  histo1d.cc                  histo1d.h        \
  histo2d.cc                  histo2d.h        \
  histo4d.cc                  histo4d.h        \
  stringlist.cc               stringlist.h     \
                              countedref.h     \
  \
  \
                              point.h          \
  basestar.cc                 basestar.h       \
  jimstar.cc                  jimstar.h        \
  nstarmatch.cc               nstarmatch.h     \
  sestar.cc                   sestar.h         \
  standardstar.cc             standardstar.h   \
  starlist.cc                 starlist.h       \
  starmatch.cc                starmatch.h      \
  listmatch.cc                listmatch.h      \
  fastfinder.cc               fastfinder.h     \
  usnoutils.cc                usnoutils.h      \
  \
  gtransfo.cc                 gtransfo.h       \
                              gtransfocomposition.h \
  \
  align_box.cc                align_box.h      \
  quali_box.cc                quali_box.h      \
  seeing_box.cc               seeing_box.h     \
  sextractor_box.cc           sextractor_box.h \
  astroutils.cc               astroutils.h     \
  fringeutils.cc              fringeutils.h    \
  reducedutils.cc             reducedutils.h   \
  superflat.cc                superflat.h      \
  \
  frame.cc                    frame.h            \
  image.cc                    image.h            \
  dimage.cc                   dimage.h           \
  detection.cc                detection.h        \
  allreducedimage.cc          allreducedimage.h  \
  dbimage.cc                  dbimage.h          \
  dbfileParse.y                          \
  dbfileLex.l                            \
  fitsimage.cc                fitsimage.h        \
  fitsset.cc                  fitsset.h          \
  fitsslice.cc                fitsslice.h        \
  fitsimagearray.cc           fitsimagearray.h   \
  imageback.cc                imageback.h        \
  imagelist.h        \
  imagematch.cc               imagematch.h       \
  imagesubtraction.cc         imagesubtraction.h \
  sub.cc                      sub.h              \
  imagesum.cc                 imagesum.h         \
  reducedimage.cc             reducedimage.h     \
  subimage.cc                 subimage.h         \
  transformedimage.cc         transformedimage.h \
  fitstoad.cc                 fitstoad.h \
  cluster.cc                  cluster.h  \
  convolution.cc              convolution.h \
                              convtemplate.h \
  kernelfit.cc                kernelfit.h \
  psfmatch.cc                 psfmatch.h  \
  swarpstack.cc		      swarpstack.h \
  dicstar.cc                  dicstar.h \
  dictfile.cc		      dictfile.h



#  imagesum_dict.cc \
#  imagesubtraction_dict.cc \
#  subimage_dict.cc \
#  swarpstack_dict.cc \
#  transformedimage_dict.cc 

include_HEADERS = \
dictfile.h \
dicstar.h \
vutils.h         \
wcsutils.h       \
toadscards.h     \
datacards.h      \
dynccd.h         \
datdetec.h       \
fileutils.h      \
rootstuff.h      \
histo1d.h        \
histo2d.h        \
histo4d.h        \
stringlist.h     \
countedref.h     \
point.h          \
basestar.h       \
jimstar.h        \
nstarmatch.h     \
sestar.h         \
standardstar.h   \
starlist.h       \
starmatch.h      \
listmatch.h      \
fastfinder.h     \
usnoutils.h      \
gtransfo.h       \
gtransfocomposition.h \
align_box.h      \
quali_box.h      \
seeing_box.h     \
sextractor_box.h \
astroutils.h     \
fringeutils.h    \
reducedutils.h   \
superflat.h      \
frame.h            \
image.h            \
dimage.h           \
detection.h        \
dodetection.h      \
allreducedimage.h  \
dbimage.h          \
fitsimage.h        \
fitsset.h          \
fitsslice.h        \
fitsimagearray.h   \
imageback.h        \
imagelist.h        \
imagematch.h       \
imagesubtraction.h \
sub.h              \
imagesum.h         \
reducedimage.h     \
subimage.h         \
transformedimage.h \
fitstoad.h \
cluster.h  \
convolution.h \
convtemplate.h \
kernelfit.h \
psfmatch.h  \
swarpstack.h    \
night.h     \
nightelement.h  \
matvect.h



noinst_HEADERS = wcscon.h defs.h simps2d.h


CLEANFILES = *_dict.h *_dict.cc *__persister.h *__persister.cc

######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package

if AUTOINSTALL

@prefix@/bin/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* @prefix@/bin/$*

@prefix@/lib/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* @prefix@/lib/$*

@prefix@/include/% : @abs_builddir@/%
	ln -fs @abs_builddir@/$* @prefix@/include/$*

all : $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/libtoads.*))  $(patsubst %,@prefix@/include/%,${include_HEADERS})

endif
