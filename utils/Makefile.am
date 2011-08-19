CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in \
	stamp-*
DISTCLEANFILES = mat.fits

AM_CPPFLAGS = \
	-I$(top_srcdir)/src_base \
	-I$(top_srcdir)/src_utils \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/src_base \
	-I$(top_srcdir)/flat \
	-I$(top_srcdir)/psf \
	-I$(top_srcdir)/lc \
	-I$(top_srcdir)/simphot \
	@CFITSIO_CFLAGS@

bin_PROGRAMS = \
	header astrom changekey checker dbimstat \
	dead_correction elixir_to_toads \
        findsatellite image_align image_register \
        image_select image_stats image_sub image_sum \
        lbl2toads maglimeval make_catalog \
        make_cosmic make_dead joinlist make_weight \
        matchusno preptoads removeback selection \
        skylev sum_images valid_fits make_lightcurve \
        rustine_for_fringe assign dbls split_fits \
        fitnight invertdeadmap install_image \
        vignet4gif modelvignet wrapswarp compute \
        enlargeimage calibrate coordfromimage \
        enlargedead mask2dead mkcat2 seeing \
        mergelists modelpsf zpfitter \
        calibrationhacker basicaperture \
        calibratefixedpos flatfield makepsf \
        compute_zp splitsubs test_dust makvar \
        getvalueatpoint catalog_images_match \
        deadprop checkmask sumpixels wcs_for_vignets \
        apercorr simple_zpfitter satellite_trajectory \
        merge_headers image_match image_combine image_subtract \
	daosetup daoconvert makesub make_satur addstars magellan_to_poloka

LDADD = $(top_builddir)/src/libtoads.la

header_SOURCES = header.cc
astrom_SOURCES = astrom.cc
changekey_SOURCES = changekey.cc
checker_SOURCES = checker.cc
dbimstat_SOURCES = dbimstat.cc
dead_correction_SOURCES = dead_correction.cc
elixir_to_toads_SOURCES = elixir_to_toads.cc
elixir_to_toads_LDADD = $(LDADD) $(top_builddir)/flat/libtoadsflat.la
magellan_to_poloka_SOURCES = magellan_to_poloka.cc
findsatellite_SOURCES = findsatellite.cc
image_align_SOURCES = image_align.cc
image_register_SOURCES = image_register.cc
image_select_SOURCES = image_select.cc
image_stats_SOURCES = image_stats.cc
image_sub_SOURCES = image_sub.cc
image_sum_SOURCES = image_sum.cc
lbl2toads_SOURCES = lbl2toads.cc
maglimeval_SOURCES = maglimeval.cc
make_catalog_SOURCES = make_catalog.cc
make_cosmic_SOURCES = make_cosmic.cc
make_dead_SOURCES = make_dead.cc
joinlist_SOURCES = joinlist.cc
make_weight_SOURCES = make_weight.cc
matchusno_SOURCES = matchusno.cc
preptoads_SOURCES = preptoads.cc
removeback_SOURCES = removeback.cc
selection_SOURCES = selection.cc
selection_LDADD = $(LDADD) $(top_builddir)/flat/libtoadsflat.la
skylev_SOURCES = skylev.cc
sum_images_SOURCES = sum_images.cc
valid_fits_SOURCES = valid_fits.cc
make_lightcurve_SOURCES = make_lightcurve.cc
make_lightcurve_LDADD = $(LDADD) $(top_builddir)/lc/libtoadslc.la
rustine_for_fringe_SOURCES = rustine_for_fringe.cc
assign_SOURCES = assign.cc
dbls_SOURCES = dbls.cc
split_fits_SOURCES = split_fits.cc
fitnight_SOURCES = fitnight.cc
fitnight_LDADD = $(LDADD) $(top_builddir)/lc/libtoadslc.la
invertdeadmap_SOURCES = invertdeadmap.cc
install_image_SOURCES = install_image.cc
vignet4gif_SOURCES = vignet4gif.cc
modelvignet_SOURCES = modelvignet.cc
wrapswarp_SOURCES = wrapswarp.cc
compute_SOURCES = compute.cc
enlargeimage_SOURCES = enlargeimage.cc
calibrate_SOURCES = calibrate.cc
calibrate_LDADD = $(LDADD) $(top_builddir)/lc/libtoadslc.la
coordfromimage_SOURCES = coordfromimage.cc
enlargedead_SOURCES = enlargedead.cc
mask2dead_SOURCES = mask2dead.cc
mkcat2_SOURCES = mkcat2.cc
seeing_SOURCES = seeing.cc
mergelists_SOURCES = mergelists.cc
modelpsf_SOURCES = modelpsf.cc
modelpsf_LDADD = $(LDADD) $(top_builddir)/lc/libtoadslc.la
zpfitter_SOURCES = zpfitter.cc
calibrationhacker_SOURCES = calibrationhacker.cc
basicaperture_SOURCES = basicaperture.cc
calibratefixedpos_SOURCES = calibratefixedpos.cc
calibratefixedpos_LDADD = $(LDADD) $(top_builddir)/lc/libtoadslc.la
flatfield_SOURCES = flatfield.cc
flatfield_LDADD = $(LDADD) $(top_builddir)/flat/libtoadsflat.la
makepsf_SOURCES = makepsf.cc
makepsf_LDADD = $(LDADD) $(top_builddir)/psf/libtoadspsf.la
compute_zp_SOURCES = compute_zp.cc
splitsubs_SOURCES = splitsubs.cc
test_dust_SOURCES = test_dust.cc
makvar_SOURCES = makvar.cc
getvalueatpoint_SOURCES = getvalueatpoint.cc
catalog_images_match_SOURCES = catalog_images_match.cc
deadprop_SOURCES = deadprop.cc
checkmask_SOURCES = checkmask.cc
sumpixels_SOURCES = sumpixels.cc
wcs_for_vignets_SOURCES = wcs_for_vignets.cc
apercorr_SOURCES = apercorr.cc
simple_zpfitter_SOURCES = simple_zpfitter.cc
satellite_trajectory_SOURCES = satellite_trajectory.cc
merge_headers_SOURCES = merge_headers.cc
image_combine_SOURCES = image_combine.cc
image_subtract_SOURCES = image_subtract.cc
image_match_SOURCES = image_match.cc
daoconvert_SOURCES = daoconvert.cc
daosetup_SOURCES = daosetup.cc
makesub_SOURCES = makesub.cc
make_satur_SOURCES = make_satur.cc
addstars_SOURCES = addstars.cc
addstars_LDADD = $(LDADD) $(top_builddir)/psf/libtoadspsf.la
