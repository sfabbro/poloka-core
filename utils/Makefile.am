EXTRA_DIST =

bin_PROGRAMS = \
 header astrom changekey checker create_precise_std_file \
 dbimstat dead_correction defringe defringe2  \
 elixir_to_toads findsatellite flatfield fringefinder \
 image_align image_combine image_register image_select \
 image_stats image_sub image_sum lbl2toads maglimeval \
 make_catalog make_cosmic make_phot \
 make_dead makefringe \
 joinlist \
 makesub make_weight matchNimages matchradec matchstandard matchusno \
 preptoads removeback selection \
 skylev sum_images trim valid_fits make_daophot make_lightcurve \
 split_fits2 rustine_for_fringe assign dbls split_fits testdao \
 fitnight invertdeadmap install_image vignet4gif modelvignet wrapswarp \
 compute

### CPPFLAGS for all programs
AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/objio @XML_CFLAGS@ @SEX_CPPFLAGS@ @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@
AM_LDFLAGS = -L$(top_builddir)/src  -L$(top_builddir)/objio  -L$(top_builddir)/lc -L$(top_builddir)/dao -L$(top_builddir)/cern_stuff -L@libdir@  -ltoads -ltoadio  -lcernstuff -llapack -lblas @CMT_LDFLAGS@  @XML_LIBS@ @SEX_LDFLAGS@  @CFITSIO_LDFLAGS@ -lg2c



### Programs depending on dao and lc
make_lightcurve_SOURCES = make_lightcurve.cc
make_lightcurve_CXXFLAGS = -I$(top_srcdir)/dao -I$(top_srcdir)/lc $(AM_CPPFLAGS)
make_lightcurve_LDFLAGS =  -ltoadslc -ltoadsdao

make_phot_SOURCES = make_phot.cc
make_phot_CXXFLAGS = -I$(top_srcdir)/dao -I$(top_srcdir)/lc $(AM_CPPFLAGS)
make_phot_LDFLAGS =  -ltoadslc -ltoadsdao

fitnight_SOURCES = fitnight.cc
fitnight_CXXFLAGS = -I$(top_srcdir)/dao -I$(top_srcdir)/lc $(AM_CPPFLAGS)
fitnight_LDFLAGS =  -ltoadslc -ltoadsdao

modelvignet_SOURCES = modelvignet.cc
modelvignet_CXXFLAGS = -I$(top_srcdir)/dao -I$(top_srcdir)/lc $(AM_CPPFLAGS)
modelvignet_LDFLAGS =  -ltoadslc -ltoadsdao


### Programs depending on dao
make_daophot_SOURCES = make_daophot.cc
make_daophot_CXXFLAGS = -I$(top_srcdir)/dao $(AM_CPPFLAGS)
make_daophot_LDFLAGS =  -ltoadsdao

testdao_SOURCES = testdao.cc
testdao_CXXFLAGS = -I$(top_srcdir)/dao $(AM_CPPFLAGS)
testdao_LDFLAGS =  -ltoadsdao


compute_SOURCES = compute.cc
wrapswarp_SOURCES = wrapswarp.cc
joinlist_SOURCES = joinlist.cc
vignet4gif_SOURCES = vignet4gif.cc
install_image_SOURCES = install_image.cc
invertdeadmap_SOURCES = invertdeadmap.cc
assign_SOURCES = assign.cc
dbls_SOURCES = dbls.cc
split_fits_SOURCES = split_fits.cc
rustine_for_fringe_SOURCES = rustine_for_fringe.cc
split_fits2_SOURCES = split_fits2.cc
header_SOURCES = header.cc
astrom_SOURCES = astrom.cc 
changekey_SOURCES = changekey.cc
checker_SOURCES = checker.cc
create_precise_std_file_SOURCES = create_precise_std_file.cc
dbimstat_SOURCES = dbimstat.cc
dead_correction_SOURCES = dead_correction.cc
defringe_SOURCES = defringe.cc
defringe2_SOURCES = defringe2.cc
elixir_to_toads_SOURCES = elixir_to_toads.cc
findsatellite_SOURCES = findsatellite.cc
flatfield_SOURCES = flatfield.cc
fringefinder_SOURCES = fringefinder.cc
image_align_SOURCES = image_align.cc
image_combine_SOURCES = image_combine.cc
# image_match_SOURCES = image_match.cc
image_register_SOURCES = image_register.cc
image_select_SOURCES = image_select.cc
image_stats_SOURCES = image_stats.cc
image_sub_SOURCES = image_sub.cc
image_sum_SOURCES = image_sum.cc
lbl2toads_SOURCES = lbl2toads.cc
#lmatch_SOURCES = lmatch.cc
maglimeval_SOURCES = maglimeval.cc 
make_catalog_SOURCES = make_catalog.cc
make_cosmic_SOURCES = make_cosmic.cc
make_dead_SOURCES = make_dead.cc
makefringe_SOURCES = makefringe.cc
#make_mosaique_SOURCES = make_mosaique.cc
#make_sub_SOURCES = make_sub.cc
makesub_SOURCES = makesub.cc
make_weight_SOURCES = make_weight.cc
matchNimages_SOURCES = matchNimages.cc
matchradec_SOURCES = matchradec.cc
matchstandard_SOURCES = matchstandard.cc
matchusno_SOURCES = matchusno.cc
preptoads_SOURCES = preptoads.cc
removeback_SOURCES = removeback.cc
selection_SOURCES = selection.cc
skylev_SOURCES = skylev.cc
sum_images_SOURCES = sum_images.cc
trim_SOURCES = trim.cc
valid_fits_SOURCES = valid_fits.cc




######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package

if AUTOINSTALL

@prefix@/bin/% : %
	ln -fs @abs_builddir@/.libs/$* @prefix@/bin/$*

all   : $(patsubst %,@prefix@/bin/%,${bin_PROGRAMS})

endif
