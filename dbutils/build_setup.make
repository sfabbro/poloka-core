#
# File generated Thu Feb 25 16:20:28 MET 1999
#

all :: csh config
#	@echo csh and config ok"

csh : $(tag).csh
#	@echo $(tag).csh ok
  
$(tag).csh : $(METHODSROOT)/mgr/requirements  requirements
	@echo Rebuilding setup and config. One requirements file has changed.
	@pck build_setup -tag=$(tag)

config : Make.$(tag)
#	@echo Make.$(tag) ok

Make.$(tag) : $(METHODSROOT)/mgr/requirements  requirements
	@touch Make.$(tag); ${METHODSROOT}/mgr/pck -quiet show macros -tag=$(tag) | 	sed -e "s#[\']##g" >t$$; mv t$$ Make.$(tag)

include Make.$(tag)

configclean :
	rm -f Make.$(tag)
	rm -f build_setup.make
	rm -f build_config.make


