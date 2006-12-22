AM_FFLAGS   = -Wall -pedantic -W  -fno-f2c -O3 -funroll-loops

lib_LTLIBRARIES = liblapackstuff.la
liblapackstuff_la_SOURCES = \
		  ddot.f    \
		  dgemm.f   \
		  dgemv.f   \
		  dlauum.f  \
		  dlauu2.f  \
		  dposv.f   \
		  dpotf2.f  \
		  dpotrf.f  \
		  dpotri.f  \
		  dpotrs.f  \
		  dscal.f   \
		  dsyrk.f   \
		  dtrti2.f  \
		  dtrmm.f   \
		  dtrmv.f   \
		  dtrsm.f   \
		  dtrtri.f  \
		  ieeeck.f  \
		  ilaenv.f  \
		  lsame.f   \
		  dgesv.f \
		  dgetf2.f \
		  dgetrf.f \
		  dgetrs.f \
		  dlaswp.f \
	dsysv.f \
	dsytrf.f \
	dsytrs.f \
	dswap.f \
	dlasyf.f \
	dcopy.f \
	idamax.f \
	dger.f \
	dsytf2.f \
	dsyr.f \
	xerbla.f  

######## auto symbolic link #########
# if configure is set with option --enable-autoinstall
# symbolic links to all binaries are set in the install directory
# this prevents during developpement to install
# systematically the package

if AUTOINSTALL

@prefix@/lib/% : .libs/%
	ln -fs @abs_builddir@/.libs/$* @prefix@/lib/$*

all :  $(patsubst %,@prefix@/bin/%,${bin_PROGRAMS}) $(patsubst .libs/%,@prefix@/lib/%,$(wildcard .libs/liblapackstuff.*)) $(patsubst %,@prefix@/include/%,${include_HEADERS})

endif