EXTRA_DIST =   datacards/default.conv  datacards/default.nnw      datacards/default.param  \
               datacards/default.sex   datacards/seeing.datacard  datacards/sub.datacard \
               config/config.guess     config/config.sub          config/install-sh \
               config/ltmain.sh        config/missing             config/mkinstalldirs \
               repack-config.in        repack-config \
               repack.spec.in          repack.spec

SUBDIRS = objio cern_stuff lapack_stuff telinst src dao lc mc utils datacards cern_utils


#bindings doc  

bin_SCRIPTS = repack-config 
BUILT_SOURCES = repack-config
repack-config: repack-config.in  my-installdirs

include_HEADERS = config.h 

## certainly exists in automake but I just do it by hand

install_dirs = @prefix@/bin @prefix@/lib @prefix@/include
my-installdirs :
	@for dir in $(install_dirs); do \
	  test -z "$$dir" || $(mkdir_p) "$$dir"; \
	done
