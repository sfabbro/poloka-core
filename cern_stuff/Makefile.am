EXTRA_DIST = kernnum/pilot.h

lib_LTLIBRARIES = libcernstuff.la
libcernstuff_la_SOURCES = \
 dfact.F  dfinv.F  eisrs1d.F  feqn.inc  sfact.inc  tmprnt.F  tred2d.F \
 dfeqn.F  dsinv.F  fact.inc  finv.inc  sfactd.inc   sfinv.inc  tql2d.F

######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package


if AUTOINSTALL

@prefix@/lib/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* @prefix@/lib/$*

all :  $(patsubst %,@prefix@/bin/%,${bin_PROGRAMS}) $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/libcernstuff.*)) $(patsubst %,@prefix@/include/%,${include_HEADERS})

endif

