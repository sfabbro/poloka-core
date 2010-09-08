
#AM_CPPFLAGS = -pg
#AM_LDFLAGS = -pg

## 
## EXTRA_DIST =   datacards/default.conv  datacards/default.nnw      datacards/default.param  \
##                datacards/default.sex   datacards/seeing.datacard  datacards/sub.datacard \
##                config/config.guess     config/config.sub          config/install-sh \
##                config/ltmain.sh        config/missing             config/mkinstalldirs \
##                repack-config.in        repack-config \
##                repack.spec.in          repack.spec
## 

EXTRA_DIST =   autogen.sh datacards/default.conv datacards/default.nnw \
	datacards/default.param  datacards/default.sex   \
	datacards/seeing.datacard  datacards/sub.datacard \
	cmt/cmt-config.csh



SUBDIRS = cern_stuff lapack_stuff telinst src_base src_utils flat src psf simphot lc utils datacards cern_utils imagemagick_utils
if ENABLE_MC
SUBDIRS += mc
endif


## bindings doc  

#bin_SCRIPTS = repack-config 
#BUILT_SOURCES = repack-config
#repack-config: repack-config.in  my-installdirs

all : my-installdirs

include_HEADERS = config.h 

##
## certainly exists in automake but I just do it by hand
##
install_dirs = @prefix@/bin @prefix@/lib @prefix@/include
my-installdirs : 	
	@for dir in $(install_dirs); do \
	  test -z "$$dir" || $(mkdir_p) "$$dir"; \
	done
	cd @prefix@/lib; test -e ".libs" || ln -s . .libs

