EXTRA_DIST = slasub.c alltelinst.cc alltelinst.h wcs_defines.h

# with CMT
#AM_CPPFLAGS = @CMT_INCLUDES@ 
#AM_LDFLAGS = @CMT_LDFLAGS@ @CERNLIB_LDFLAGS@

# without CMT
AM_CPPFLAGS = -I../objio @XML_CFLAGS@ @SEX_CPPFLAGS@ @CFITSIO_CPPFLAGS@ @LAPACK_CPPFLAGS@
AM_LDFLAGS = -L../objio -ltoadio @XML_LIBS@ @SEX_LDFLAGS@ @CFITSIO_LDFLAGS@ @LAPACK_LDFLAGS@ @CERNLIB_LDFLAGS@

AM_YFLAGS = -l -d

%_dict.cc : %.h
	../objio/dictgen2 -o $@ $<


#  psf.cc                      psf.h            



lib_LTLIBRARIES = libtoads.la
libtoads_la_SOURCES = \
 \
  frame_dict.cc \
  point_dict.cc \
  basestar_dict.cc \
  sestar_dict.cc \
  candidatestar_dict.cc \
  fidstar_dict.cc \
  jimstar_dict.cc \
  photstar_dict.cc \
  rollingstar_dict.cc \
  senearstar_dict.cc \
  standardstar_dict.cc \
  starmatch_dict.cc \
  yquemstar_dict.cc \
  histo1d_dict.cc \
  histo2d_dict.cc \
  histo4d_dict.cc \
  image_dict.cc \
  dimage_dict.cc \
 \
  vutils.cc                   vutils.h         \
  wcscon.c                    wcscon.h         \
  wcsutils.cc                 wcsutils.h       \
  toadscards.cc               toadscards.h     \
  datacards.cc                datacards.h      \
  dynccd.cc                   dynccd.h         \
  datdetec.cc                 datdetec.h       \
  fileutils.cc                fileutils.h      \
  rootstuff.cc                rootstuff.h      \
  histo1d.cc                  histo1d.h        \
  histo2d.cc                  histo2d.h        \
  histo4d.cc                  histo4d.h        \
  stringlist.cc               stringlist.h     \
  nbrandom.c                  nbrandom.h       \
                              nbconst.h        \
                              countedref.h     \
                              exceptions.h     \
  \
  \
                              point.h          \
  basestar.cc                 basestar.h       \
  candidatestar.cc            candidatestar.h  \
  candstar.cc                 candstar.h       \
  fidstar.cc                  fidstar.h        \
  jimstar.cc                  jimstar.h        \
  nstarmatch.cc               nstarmatch.h     \
  photstar.cc                 photstar.h       \
  rollingstar.cc              rollingstar.h    \
  senearstar.cc               senearstar.h     \
  sestar.cc                   sestar.h         \
  standardstar.cc             standardstar.h   \
  starlist.cc                 starlist.h       \
  starmatch.cc                starmatch.h      \
  yquemstar.cc                yquemstar.h      \
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
  dodetection.cc              dodetection.h      \
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
  transformedimageaddfakes.cc transformedimageaddfakes.h \
  transformedimage.cc         transformedimage.h \
  fitstoad.cc                 fitstoad.h \
  addfakes.cc                 addfakes.h \
  cluster.cc                  cluster.h  \
  convolution.cc              convolution.h \
                              convtemplate.h \
  kernelfit.cc                kernelfit.h \
  psfmatch.cc                 psfmatch.h  \
  \
  newsub.cc                   newsub.h    \
  night.cc                    night.h     \
                              nightelement.h

#  simulation.cc               simulation.h                       




# psf.h 
include_HEADERS = \
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
nbrandom.h       \
nbconst.h        \
countedref.h     \
exceptions.h     \
 \
 \
point.h          \
basestar.h       \
candidatestar.h  \
candstar.h       \
fidstar.h        \
jimstar.h        \
nstarmatch.h     \
photstar.h       \
rollingstar.h    \
senearstar.h     \
sestar.h         \
standardstar.h   \
starlist.h       \
starmatch.h      \
yquemstar.h      \
listmatch.h      \
fastfinder.h     \
usnoutils.h      \
 \
gtransfo.h       \
gtransfocomposition.h \
 \
align_box.h      \
quali_box.h      \
seeing_box.h     \
sextractor_box.h \
astroutils.h     \
fringeutils.h    \
reducedutils.h   \
superflat.h      \
		 \
frame.h            \
image.h            \
dimage.h           \
detection.h        \
dodetection.h      \
allreducedimage.h  \
dbimage.h          \
           \
           \
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
transformedimageaddfakes.h \
transformedimage.h \
fitstoad.h \
addfakes.h \
cluster.h  \
convolution.h \
convtemplate.h \
kernelfit.h \
psfmatch.h  \
 \
newsub.h    \
night.h     \
nightelement.h 

#simulation.h                       

noinst_HEADERS = wcscon.h defs.h simps2d.h


CLEANFILES = *_dict.h *_dict.cc *__persister.h *__persister.cc

