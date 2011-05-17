datacardsdir = $(datadir)/@PACKAGE_NAME@/datacards
datacards_DATA = \
	default.conv \
	default.nnw \
	default.param \
	default.sex \
	match.datacard \
	seeing.datacard \
	sub.datacard

dbconfigdir = $(datadir)/@PACKAGE_NAME@
dbconfig_DATA = dbconfig

EXTRA_DIST = $(datacards_DATA) $(dbconfig_DATA)

MAINTAINERCLEANFILES = Makefile.in
CLEANFILES = *~
