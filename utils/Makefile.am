EXTRA_DIST =

bin_PROGRAMS = \
 header astrom changekey checker create_precise_std_file \
 dbimstat dead_correction defringe defringe2 detectonedb \
 elixir_to_toads findsatellite flatfield fringefinder \
 image_align image_combine image_register image_select \
 image_stats image_sub image_sum lbl2toads maglimeval \
 make_catalog make_cosmic \
 make_dead makefringe \
 makesub make_weight matchNimages matchradec matchstandard matchusno \
 newmake_sub preptoads removeback selection \
 skylev sum_images trim valid_fits wiynpad


AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/objio @XML_CFLAGS@ @SEX_CPPFLAGS@ @CFITSIO_CPPFLAGS@ @LAPACK_CPPFLAGS@ @CMT_CFLAGS@
AM_LDFLAGS = -L$(top_builddir)/src  -L$(top_builddir)/objio -L$(top_builddir)/cern_stuff -L@libdir@ -ltoads -ltoadio   -lcernstuff  @CMT_LDFLAGS@  @XML_LIBS@ @SEX_LDFLAGS@  @CFITSIO_LDFLAGS@ @LAPACK_LDFLAGS@ -lg2c

header_SOURCES = header.cc
astrom_SOURCES = astrom.cc 
changekey_SOURCES = changekey.cc
checker_SOURCES = checker.cc
create_precise_std_file_SOURCES = create_precise_std_file.cc
dbimstat_SOURCES = dbimstat.cc
dead_correction_SOURCES = dead_correction.cc
defringe_SOURCES = defringe.cc
defringe2_SOURCES = defringe2.cc
detectonedb_SOURCES = detectonedb.cc
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
#make_daophot_SOURCES = make_daophot.cc
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
newmake_sub_SOURCES = newmake_sub.cc
preptoads_SOURCES = preptoads.cc
removeback_SOURCES = removeback.cc
selection_SOURCES = selection.cc
skylev_SOURCES = skylev.cc
sum_images_SOURCES = sum_images.cc
trim_SOURCES = trim.cc
valid_fits_SOURCES = valid_fits.cc
wiynpad_SOURCES = wiynpad.cc
#write_star_SOURCES = write_star.cc
