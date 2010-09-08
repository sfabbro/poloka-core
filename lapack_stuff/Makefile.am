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
	dsysv.f 	\
	dsytrf.f 	\
	dsytrs.f 	\
	dswap.f 	\
	dlasyf.f 	\
	dcopy.f 	\
	idamax.f 	\
	dger.f 		\
	dsytf2.f 	\
	dsyr.f 		\
	xerbla.f  	\
	dsytri.f 	\
	dsyevd.f 	\
	dormtr.f 	\
	dsytrd.f 	\
	dlascl.f 	\
	dlansy.f 	\
	dsymv.f 	\
	dsterf.f 	\
	dstedc.f 	\
	dlacpy.f 	\
	dlamch.f 	\
	dlapy2.f  	\
	dormql.f  	\
	dlatrd.f  	\
	dlasrt.f  	\
	dsteqr.f  	\
	dlaed0.f  	\
	dormqr.f  	\
	dsyr2k.f  	\
	dlaset.f  	\
	dlassq.f  	\
	dsytd2.f  	\
	dlae2.f  	\
	dlanst.f  	\
	dlarfb.f   	\
	dlasr.f   	\
	dsyr2.f   	\
	dlaed1.f   	\
	dlaev2.f   	\
	dlarft.f   	\
	dorm2l.f   	\
	dlarfg.f   	\
	dorm2r.f   	\
	dlartg.f   	\
	dlaed7.f   	\
	daxpy.f   	\
	dlaed9.f   	\
	dlamrg.f   	\
	dlaed3.f   	\
	dlarf.f   	\
	dlaed8.f   	\
	dlaed2.f   	\
	dlaeda.f   	\
	dnrm2.f   	\
	drot.f 		\
	dlaed4.f	\
	dlaed5.f	\
	dlaed6.f


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