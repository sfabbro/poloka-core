EXTRA_DIST =   datacards/default.conv  datacards/default.nnw      datacards/default.param  \
               datacards/default.sex   datacards/seeing.datacard  datacards/sub.datacard \
               config/config.guess     config/config.sub          config/install-sh \
               config/ltmain.sh        config/missing             config/mkinstalldirs \
               repack-config.in        repack-config \
               repack.spec.in          repack.spec

SUBDIRS = src cern_stuff dao_stuff utils datacards objio bindings 

bin_SCRIPTS = repack-config
BUILT_SOURCES = repack-config
repack-config: repack-config.in
