CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in \
	stamp-*

AM_CPPFLAGS = \
	-I$(top_srcdir)/src_base \
	-I$(top_srcdir)/src_utils \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src_base \
	@CFITSIO_CFLAGS@

LDADD = $(top_builddir)/src/libtoads.la

bin_PROGRAMS = \
	apercorr assign astrom \
	basicaperture \
	calibratedbim calibrationhacker catalog_images_match changekey checkelixirflats \
	checker check_fits checkmask compute compute_zp coordfromimage \
	dbimstat dbls dead_correction deadprop \
	enlargedead enlargeimage \
	findsatellite getvalueatpoint gridcorrection \
	header \
	image_align image_register image_select image_stats image_sub image_sum \
	install_image invertdeadmap \
	joinlist \
	lbl2toads \
	maglimeval make_catalog make_cosmic make_dead make_weight makesub \
	makvar mapfcangles mask2dead \
	matchusno mergelists mkcat2 mosaique_int \
	preptoads \
	removeback rustine_for_fringe \
	seeing simple_zpfitter skylev split_fits splitsubs \
	sum_images sumpixels \
	valid_fits vignet4gif \
	wcs_for_vignets \
	wrapswarp \
	zpfitter

check_PROGRAMS = testdao test_dust test_matvect testphot

TESTS = $(check_PROGRAMS)

apercorr_SOURCES = apercorr.cc
assign_SOURCES = assign.cc
astrom_SOURCES = astrom.cc
basicaperture_SOURCES = basicaperture.cc
calibratedbim_SOURCES = calibratedbim.cc
calibrationhacker_SOURCES = calibrationhacker.cc
catalog_images_match_SOURCES = catalog_images_match.cc
changekey_SOURCES = changekey.cc
checkelixirflats_SOURCES = checkelixirflats.cc
checker_SOURCES = checker.cc
check_fits_SOURCES = check_fits.cc
checkmask_SOURCES = checkmask.cc
compute_SOURCES = compute.cc
compute_zp_SOURCES = compute_zp.cc
coordfromimage_SOURCES = coordfromimage.cc
dbimstat_SOURCES = dbimstat.cc
dbls_SOURCES = dbls.cc
dead_correction_SOURCES = dead_correction.cc
deadprop_SOURCES = deadprop.cc
enlargedead_SOURCES = enlargedead.cc
enlargeimage_SOURCES = enlargeimage.cc
findsatellite_SOURCES = findsatellite.cc
getvalueatpoint_SOURCES = getvalueatpoint.cc
gridcorrection_SOURCES = gridcorrection.cc
header_SOURCES = header.cc
image_align_SOURCES = image_align.cc
image_register_SOURCES = image_register.cc
image_select_SOURCES = image_select.cc
image_stats_SOURCES = image_stats.cc
image_sub_SOURCES = image_sub.cc
image_sum_SOURCES = image_sum.cc
install_image_SOURCES = install_image.cc
invertdeadmap_SOURCES = invertdeadmap.cc
joinlist_SOURCES = joinlist.cc
lbl2toads_SOURCES = lbl2toads.cc
maglimeval_SOURCES = maglimeval.cc
make_catalog_SOURCES = make_catalog.cc
make_cosmic_SOURCES = make_cosmic.cc
make_dead_SOURCES = make_dead.cc
makesub_SOURCES = makesub.cc
make_weight_SOURCES = make_weight.cc
makvar_SOURCES = makvar.cc
mapfcangles_SOURCES = mapfcangles.cc
mask2dead_SOURCES = mask2dead.cc
matchusno_SOURCES = matchusno.cc
mergelists_SOURCES = mergelists.cc
mkcat2_SOURCES = mkcat2.cc
mosaique_int_SOURCES = mosaique_int.cc
preptoads_SOURCES = preptoads.cc
removeback_SOURCES = removeback.cc
rustine_for_fringe_SOURCES = rustine_for_fringe.cc
seeing_SOURCES = seeing.cc
simple_zpfitter_SOURCES = simple_zpfitter.cc
skylev_SOURCES = skylev.cc
split_fits_SOURCES = split_fits.cc
splitsubs_SOURCES = splitsubs.cc
sum_images_SOURCES = sum_images.cc
sumpixels_SOURCES = sumpixels.cc
testdao_SOURCES = testdao.cc
test_dust_SOURCES = test_dust.cc
test_matvect_SOURCES = test_matvect.cc
testphot_SOURCES = testphot.cc
valid_fits_SOURCES = valid_fits.cc
vignet4gif_SOURCES = vignet4gif.cc
wcs_for_vignets_SOURCES = wcs_for_vignets.cc
wrapswarp_SOURCES = wrapswarp.cc
zpfitter_SOURCES = zpfitter.cc
