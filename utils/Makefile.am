## these programs are not compiled if sextractor is disable ( --disable-sex in configure)
if ENABLE_MC
SNFIT_INC = -I/afs/in2p3.fr/throng/snovae/softsnls/snfit/v1r0/src
SNFIT_LIB = -L/afs/in2p3.fr/throng/snovae/softsnls/snfit/v1r0/$(CMTCONFIG)/ -lsnfitlib `cernlib packlib mathlib` -llapack -lblas
endif


if USESEXTRACTOR

EXTRA_DIST =


bin_PROGRAMS = \
 header astrom changekey checker \
 dbimstat dead_correction  \
 elixir_to_toads findsatellite \
 image_align image_register image_select \
 image_stats image_sub image_sum lbl2toads maglimeval \
 make_catalog make_cosmic \
 make_dead  \
 joinlist \
 make_weight matchusno \
 preptoads removeback selection \
 skylev sum_images valid_fits make_lightcurve \
 rustine_for_fringe assign dbls split_fits  \
 fitnight invertdeadmap install_image vignet4gif modelvignet wrapswarp \
 compute enlargeimage calibrate coordfromimage enlargedead mask2dead \
 mkcat2 seeing mergelists \
 modelpsf zpfitter gridcorrection calibrationhacker basicaperture calibratefixedpos flatfield\
 makepsf compute_zp  splitsubs test_dust makvar getvalueatpoint\
 catalog_images_match \
 deadprop checkmask sumpixels \
 wcs_for_vignets \
 apercorr simple_zpfitter compute_zp_phot \
 phot_cmp phot_cmp_snlc phot_cmp_snlc_fitnight
 #list_sne_cmp 


if ENABLE_MC
bin_PROGRAMS += make_mc_cat makesub  make_dualsub  mc_addsn make_lcsim
endif

wcs_for_vignets_SOURCES = wcs_for_vignets.cc

#testback_SOURCES = testback.cc
#testphot_SOURCES = testphot.cc
#testphot_LDFLAGS = `cernlib packlib`
#fitlin_SOURCES  = fitlin.cc
#fitlin_LDFLAGS = -ltoadspsf

# split_fits2

### CPPFLAGS for all programs
AM_CPPFLAGS = -I$(top_srcdir)/src_base  -I$(top_srcdir)/src_utils  -I$(top_srcdir)/src  -I$(top_srcdir)/src_base  -I$(top_srcdir)/src_utils  -I$(top_srcdir)/src -I$(top_srcdir)/dao -I$(top_srcdir)/psf @SEX_CPPFLAGS@ @CFITSIO_CPPFLAGS@ @CMT_CFLAGS@
AM_LDFLAGS = -L$(top_builddir)/mc -L$(top_builddir)/lc -L$(top_builddir)/dao -L$(top_builddir)/src -L$(top_builddir)/flat -L$(top_builddir)/src_utils -L$(top_builddir)/src_base  -L$(top_builddir)/psf -L$(top_builddir)/cern_stuff  -L$(top_builddir)/lapack_stuff  -R @libdir@ -L@libdir@  -ltoads -ltoadsutils -ltoadsbase  -lcernstuff -llapack @CMT_LDFLAGS@  @SEX_LDFLAGS@  @CFITSIO_LDFLAGS@ -lgfortran

apercorr_SOURCES = apercorr.cc
makepsf_SOURCES = makepsf.cc
makepsf_LDFLAGS = -ltoadspsf 

catalog_images_match_SOURCES = catalog_images_match.cc

#reduce_lmatch_SOURCES = reduce_lmatch.cc
### Programs additionnaly depending on flat
elixir_to_toads_SOURCES  = elixir_to_toads.cc
elixir_to_toads_CXXFLAGS = -I$(top_srcdir)/flat
elixir_to_toads_LDFLAGS  = -ltoadsflat

selection_SOURCES  = selection.cc
selection_CXXFLAGS = -I$(top_srcdir)/flat
selection_LDFLAGS  = -ltoadsflat

flatfield_SOURCES  = flatfield.cc
flatfield_CXXFLAGS = -I$(top_srcdir)/flat
flatfield_LDFLAGS  = -ltoadsflat



if ENABLE_MC 
### Programs additionnaly depending on mc and dao and lc

#genlcfakes_SOURCES = genlcfakes.cc
#genlcfakes_CXXFLAGS = -I$(top_srcdir)/dao -I$(top_srcdir)/mc -I$(top_srcdir)/lc $(SNFIT_INC)
#genlcfakes_LDFLAGS =  -ltoadsmc  -ltoadslc -ltoadsdao $(SNFIT_LIB)

make_lcsim_SOURCES = make_lcsim.cc
make_lcsim_CXXFLAGS =  -I$(top_srcdir)/mc -I$(top_srcdir)/dao -I$(top_srcdir)/psf $(AM_CPPFLAGS) $(SNFIT_INC)
make_lcsim_LDFLAGS =  -ltoadsmc  -ltoadsdao  -ltoadspsf $(SNFIT_LIB)

mc_addsn_SOURCES = mc_addsn.cc
mc_addsn_CXXFLAGS =  -I$(top_srcdir)/mc -I$(top_srcdir)/dao -I$(top_srcdir)/psf $(AM_CPPFLAGS) $(SNFIT_INC)
mc_addsn_LDFLAGS =  -ltoadsmc  -ltoadsdao  -ltoadspsf $(SNFIT_LIB)

make_dualsub_SOURCES = make_dualsub.cc
make_dualsub_CXXFLAGS =  -I$(top_srcdir)/mc -I$(top_srcdir)/dao -I$(top_srcdir)/psf $(AM_CPPFLAGS) $(SNFIT_INC)
make_dualsub_LDFLAGS =  -ltoadsmc  -ltoadsdao  -ltoadspsf $(SNFIT_LIB)

make_mc_cat_SOURCES = make_mc_cat.cc
make_mc_cat_CXXFLAGS =  -I$(top_srcdir)/mc -I$(top_srcdir)/dao $(AM_CPPFLAGS) $(SNFIT_INC)
make_mc_cat_LDFLAGS =  -ltoadsmc  -ltoadsdao  -ltoadspsf $(SNFIT_LIB)

### Programs additionnaly depending on mc and dao
makesub_SOURCES = makesub.cc
makesub_CXXFLAGS = -I$(top_srcdir)/dao -I$(top_srcdir)/mc $(AM_CPPFLAGS) $(SNFIT_INC)
makesub_LDFLAGS =  -ltoadsmc -ltoadsdao  -ltoadspsf $(SNFIT_LIB)

endif



### Programs additionnaly depending on psf and lc
make_lightcurve_SOURCES = make_lightcurve.cc
make_lightcurve_CXXFLAGS = -I$(top_srcdir)/psf -I$(top_srcdir)/lc $(AM_CPPFLAGS)
make_lightcurve_LDFLAGS =  -ltoadslc -ltoadspsf


fitnight_SOURCES = fitnight.cc
fitnight_CXXFLAGS = -I$(top_srcdir)/psf -I$(top_srcdir)/lc $(AM_CPPFLAGS)
fitnight_LDFLAGS =  -ltoadslc -ltoadspsf

modelvignet_SOURCES = modelvignet.cc
modelvignet_CXXFLAGS = -I$(top_srcdir)/psf -I$(top_srcdir)/lc $(AM_CPPFLAGS)
modelvignet_LDFLAGS =  -ltoadslc -ltoadspsf

modelpsf_SOURCES = modelpsf.cc
modelpsf_CXXFLAGS = -I$(top_srcdir)/psf -I$(top_srcdir)/lc $(AM_CPPFLAGS)
modelpsf_LDFLAGS =  -ltoadslc -ltoadspsf

calibrate_SOURCES = calibrate.cc
calibrate_CXXFLAGS = -I$(top_srcdir)/psf -I$(top_srcdir)/lc $(AM_CPPFLAGS)
calibrate_LDFLAGS =  -ltoadslc -ltoadspsf

calibratefixedpos_SOURCES = calibratefixedpos.cc
calibratefixedpos_CXXFLAGS = -I$(top_srcdir)/psf -I$(top_srcdir)/lc $(AM_CPPFLAGS)
calibratefixedpos_LDFLAGS =  -ltoadslc -ltoadspsf

### Programs additionnaly depending on dao
#make_daophot_SOURCES = make_daophot.cc
#make_daophot_CXXFLAGS = -I$(top_srcdir)/dao $(AM_CPPFLAGS)
#make_daophot_LDFLAGS =  -ltoadsdao

getvalueatpoint_SOURCES = getvalueatpoint.cc
deadprop_SOURCES = deadprop.cc
makvar_SOURCES = makvar.cc
test_dust_SOURCES = test_dust.cc
checkmask_SOURCES = checkmask.cc
compute_zp_SOURCES = compute_zp.cc
splitsubs_SOURCES = splitsubs.cc


### All other programs
#tst_kernels_SOURCES = tst_kernels.cc
coordfromimage_SOURCES = coordfromimage.cc
enlargeimage_SOURCES = enlargeimage.cc
#readagain_extinction_SOURCES = readagain_extinction.cc
#extinction_SOURCES = extinction.cc
compute_SOURCES = compute.cc
#mergevignets_SOURCES = mergevignets.cc
wrapswarp_SOURCES = wrapswarp.cc
joinlist_SOURCES = joinlist.cc
vignet4gif_SOURCES = vignet4gif.cc
install_image_SOURCES = install_image.cc
invertdeadmap_SOURCES = invertdeadmap.cc
assign_SOURCES = assign.cc
dbls_SOURCES = dbls.cc
split_fits_SOURCES = split_fits.cc
rustine_for_fringe_SOURCES = rustine_for_fringe.cc
#split_fits2_SOURCES = split_fits2.cc
header_SOURCES = header.cc
astrom_SOURCES = astrom.cc 
changekey_SOURCES = changekey.cc
checker_SOURCES = checker.cc
dbimstat_SOURCES = dbimstat.cc
dead_correction_SOURCES = dead_correction.cc
findsatellite_SOURCES = findsatellite.cc
#flatfield_SOURCES = flatfield.cc
image_align_SOURCES = image_align.cc
#image_combine_SOURCES = image_combine.cc
#image_match_SOURCES = image_match.cc
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
#makefringe_SOURCES = makefringe.cc
#make_mosaique_SOURCES = make_mosaique.cc
#make_sub_SOURCES = make_sub.cc
make_weight_SOURCES = make_weight.cc
#matchNimages_SOURCES = matchNimages.cc
#matchradec_SOURCES = matchradec.cc
matchusno_SOURCES = matchusno.cc
preptoads_SOURCES = preptoads.cc
removeback_SOURCES = removeback.cc
skylev_SOURCES = skylev.cc
sum_images_SOURCES = sum_images.cc
valid_fits_SOURCES = valid_fits.cc
#safesex_SOURCES = safesex.cc
#essai_SOURCES = essai.cc
mkcat2_SOURCES = mkcat2.cc
enlargedead_SOURCES = enlargedead.cc
mask2dead_SOURCES = mask2dead.cc
seeing_SOURCES = seeing.cc
mergelists_SOURCES = mergelists.cc
zpfitter_SOURCES = zpfitter.cc
gridcorrection_SOURCES = gridcorrection.cc
simple_zpfitter_SOURCES = simple_zpfitter.cc
compute_zp_phot_SOURCES = compute_zp_phot.cc
phot_cmp_SOURCES = phot_cmp.cc
phot_cmp_snlc_SOURCES = phot_cmp_snlc.cc
phot_cmp_snlc_fitnight_SOURCES = phot_cmp_snlc_fitnight.cc
#list_sne_cmp_SOURCES = list_sne_cmp.cc

#calibsummary_SOURCES = calibsummary.cc
calibrationhacker_SOURCES = calibrationhacker.cc
basicaperture_SOURCES = basicaperture.cc
sumpixels_SOURCES = sumpixels.cc
testsub_SOURCES = testsub.cc
update_largeref_satur_SOURCES = update_largeref_satur.cc

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

endif
