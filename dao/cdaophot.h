#ifndef CDAOPHOT__H
#define CDAOPHOT__H

#include "cfortran.h"

#ifdef __cplusplus
extern "C" {
#endif

  
  /*

  MODIFICATIONS IN DAO (release>1.2 of cdaophot.h)

  1) The parameters MAXBOX,MAXPSF,MAXPAR,MAXEXP,MAXSTR,MAXN,MAXSKY,MAXMAX,NOPT
  describe the maximum size of arrays used in daophot.
  They MUST be fixed parameters, this is done in daocommon.f
  This latter file in included in all subroutines that need at least one
  of these parameters.
  These parameters were removed from the argument list of routines 
  (except few cases when the arrays were in the argument list)
  
  2) The values of some parameters are needed in the C++ code.
  For these purpose, routines in getdaocommon.f were added.
  
  3) C interfaces to fortran routines that include strings in the parameter list
  are build with the help of cfortran  (../src/cfortran.h) 
  (since this is quite tricky to do manually).
  In the fortran code, string arguments are defined as CHARACTER*(*) instead of 
  CHARACTER*(256), otherwise you get errors. Doing this, a STUPID2 function was 
  added with two CHARACTER*(*) arguments cause you cannot concatenate CHARACTER*(*)
  as required to use STUPID routine (absolutly thrilling).
  
  4) a -ffortran-bounds-check g77 option is added

  5) Description of fixed parameters (see daocommon.f):
  
  
  MAXBOX is the square subarray that will hold the largest final PSF.
  If the maximum PSF radius permitted is R, then MAXBOX is the
  odd integer 2*INT(R)+1.  However, because we will be dealing
  with two levels of interpolation:  (1) interpolating the raw
  picture data to arrive at a PSF whose centroid coincides with
  the central pixel of the look-up table; and (2) interpolating
  within the PSF itself to evaluate it for comparison with the
  raw picture data for the program stars, PSF will have to
  operate on a square array which is larger by 7 pixels in X
  and Y than MAXBOX.  Hence, the dimensions below are all
  MAXBOX + 7.
  
  MAXPSF is the dimension of the largest lookup table that will ever 
  need to be generated.  Recall that the corrections from the
  Gaussian approximation of the PSF to the true PSF will be
  stored in a table with a half-pixel grid size.
  MAXPSF must then equal 2*[ 2*INT(PSFRADMAX) + 1 ] + 7.
  
  MAXPAR is the maximum number of parameters which may be used to 
  describe the analytic part of the PSF, *** IN ADDITION TO
  the central intensity, and the x and y centroids.
  
  MAXEXP is the maximum number of terms in the expansion of the
  look-up table (1 for constant, 3 for linear, 6 for quadratic).
  
  MAXSTR is the largest number of stars permitted in a data file.
  
  MAXN is the maximum permitted number of PSF stars.
  
  MAXSKY is the maximum number of sky pixels we can deal with,
  given the limited amount of working space.
  
  MAXMAX is the largest group for which a solution will ever be
  attempted = maximum permissible value of MAXGRP.
  
  NOPT  is the number of options

  */
  
/*! Get sky mean, median, mode of an image */
#define GETSKY getsky_
void GETSKY(float *d, float *s, int *index, const float *readns, const float *hibad, 
	    float *skymn, float *skymed, float *skymod, float *skysig, int *n);

/*! Mean Median Mode DAOPHOT routine, for sky level and such */
#define MMM mmm_
void MMM(const float *sky, const int *nsky, const float *hibad, const float *readns, 
	 float *skymn, float *skymed, float *skymod, float *sigma, float *skew);

/*! read the PSF DAOPHOT file */
PROTOCCALLSFFUN12(INT,RDPSF,rdpsf,STRING,INT,FLOATV,INT,FLOATV,INT,INT,INT,FLOAT,FLOAT,FLOAT,FLOAT) 
#define RDPSF(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12) CCALLSFFUN12(RDPSF,rdpsf,STRING,INT,FLOATV,INT,FLOATV,INT,INT,INT,FLOAT,FLOAT,FLOAT,FLOAT,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12)

/*! returns PSF value on a point distant of (dx,dy) of a star located on (deltax,deltay) */
#define USEPSF usepsf_
float USEPSF(const int *ipstyp, const float *dx, const float *dy,  const float *bright, 
	     const float *par, const float *psf,
	     const int *npsf, const int *npar, const int *nexp, const int *nfrac, 
	     const float *deltax, const float *deltay, float *dvdxc, float *dvdyc);

/*! This is the subroutine which does the actual one-star least-squares profile fit for DAOPK */
#define PKFIT pkfit_
void PKFIT(const float *f, const int *nx, const int *ny, 
	   float *x, float *y, float *scale, const float *sky, const float *radius, 
	   const float *lobad, const float *hibad, const float *phpadu, 
	   const float *ronois, const float *perr, const float *pkerr, 
	   const float *bright, const int *ipstyp, const float *par, 
	   const int *npar, const float *psf, 
	   const int *npsf, const int *nexp, 
	   const int *nfrac, float *deltax, float *deltay, float *errmag, float *chi, 
	   float *sharp, int *niter, const float *global_sky);

/*! set the options */
PROTOCCALLSFSUB8(OPTION,option,STRING,INT,STRING,FLOATV,FLOAT,FLOAT,STRING,INT)
#define OPTION(A1,A2,A3,A4,A5,A6,A7,A8) CCALLSFSUB8(OPTION,option,STRING,INT,STRING,FLOATV,FLOAT,FLOAT,STRING,INT,A1,A2,A3,A4,A5,A6,A7,A8)

/*! find stars with gaussian filter */
PROTOCCALLSFSUB10(FIND,find,FLOATV,INT,INT,FLOATV, INTV, FLOATV, SHORTV, FLOATV,STRING,FLOAT)
#define FIND(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10) CCALLSFSUB10(FIND,find,FLOATV,INT,INT,FLOATV, INTV, FLOATV, SHORTV, FLOATV,STRING,FLOAT,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10)

/*! performs aperture photometry with daophot */
PROTOCCALLSFSUB12(PHOTSB,photsb,FLOATV,FLOATV,INTV,INT,INT,FLOAT,FLOAT,FLOAT,STRING,STRING,STRING,FLOAT)
#define PHOTSB(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12) CCALLSFSUB12(PHOTSB,photsb,FLOATV,FLOATV,INTV,INT,INT,FLOAT,FLOAT,FLOAT,STRING,STRING,STRING,FLOAT,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12)

/*! pick stars to make a psf */
PROTOCCALLSFSUB14(PCKPSF,pckpsf,INTV,FLOATV,FLOATV,FLOATV,FLOATV,INTV,INT,FLOAT,FLOAT,FLOAT,STRING,STRING,INT,FLOAT)
#define PCKPSF(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14) CCALLSFSUB14(PCKPSF,pckpsf,INTV,FLOATV,FLOATV,FLOATV,FLOATV,INTV,INT,FLOAT,FLOAT,FLOAT,STRING,STRING,INT,FLOAT,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14)

/*! get a psf from an image */
PROTOCCALLSFSUB9(GETPSF,getpsf,FLOATV,INT,INT,FLOATV,STRING,STRING,STRING,STRING,FLOAT)
#define GETPSF(A1,A2,A3,A4,A5,A6,A7,A8,A9) CCALLSFSUB9(GETPSF,getpsf,FLOATV,INT,INT,FLOATV,STRING,STRING,STRING,STRING,FLOAT,A1,A2,A3,A4,A5,A6,A7,A8,A9)

/*!  Single-star profile-fitting routine on every stars of an image */
PROTOCCALLSFSUB9(DAOPK,daopk,FLOATV,FLOAT,FLOAT,FLOAT,FLOAT,STRING,STRING,STRING,FLOAT)
#define DAOPK(A1,A2,A3,A4,A5,A6,A7,A8,A9) CCALLSFSUB9(DAOPK,daopk,FLOATV,FLOAT,FLOAT,FLOAT,FLOAT,STRING,STRING,STRING,FLOAT,A1,A2,A3,A4,A5,A6,A7,A8,A9)

/*! group stars for multi-star NSTAR */
PROTOCCALLSFSUB7(GROUP,group,FLOAT,FLOAT,STRING,STRING,STRING,FLOAT,FLOAT)
#define GROUP(A1,A2,A3,A4,A5,A6,A7) CCALLSFSUB7(GROUP,group,FLOAT,FLOAT,STRING,STRING,STRING,FLOAT,FLOAT,A1,A2,A3,A4,A5,A6,A7)

/*!  Multi-star profile-fitting routine on every stars of an image */
PROTOCCALLSFSUB11(NSTAR,nstar,FLOATV,INT,INT,FLOAT,FLOAT,FLOAT,FLOAT,STRING,STRING,STRING,FLOAT)
#define NSTAR(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11) CCALLSFSUB11(NSTAR,nstar,FLOATV,INT,INT,FLOAT,FLOAT,FLOAT,FLOAT,STRING,STRING,STRING,FLOAT,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11)

/*! subtract stars from an image using a scaled PSF */
PROTOCCALLSFSUB8(SUBSTR,substr,FLOATV,INT,INT,FLOAT,STRING,STRING,STRING,STRING)
#define SUBSTR(A1,A2,A3,A4,A5,A6,A7,A8) CCALLSFSUB8(SUBSTR,substr,FLOATV,INT,INT,FLOAT,STRING,STRING,STRING,STRING,A1,A2,A3,A4,A5,A6,A7,A8)

/*! generates fake stars from PSF and add poisson noise at locations given in a file */
PROTOCCALLSFSUB13(ADDSTR,addstr,FLOATV,INT,INT,FLOAT,STRING,STRING,STRING,STRING,FLOAT,FLOATV,INT,INT,INT)
#define ADDSTR(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13) CCALLSFSUB13(ADDSTR,addstr,FLOATV,INT,INT,FLOAT,STRING,STRING,STRING,STRING,FLOAT,FLOATV,INT,INT,INT,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13)

/*! fit all stars group by group with a PSF on an image */
PROTOCCALLSFSUB21(ALLSTR,allstr,FLOATV,INT,INT,FLOATV,FLOATV,FLOAT,FLOAT,FLOAT,INT,INT,INT,FLOAT,FLOAT,FLOAT,FLOAT,INT,STRING,STRING,STRING,STRING,FLOAT)
#define ALLSTR(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21) CCALLSFSUB21(ALLSTR,allstr,FLOATV,INT,INT,FLOATV,FLOATV,FLOAT,FLOAT,FLOAT,INT,INT,INT,FLOAT,FLOAT,FLOAT,FLOAT,INT,STRING,STRING,STRING,STRING,FLOAT,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21)

/*! dump intensities values on screen around a point */
#define DUMP dump_
void DUMP(const float *f, const int *ncol, const int *nrow, const float *coords, const float *size);

#define GETMAXBOX getmaxbox_
int GETMAXBOX();

#define GETMAXSKY getmaxsky_
int GETMAXSKY();

#define GETNOPT getnopt_
int GETNOPT();

#define GETMAXPAR getmaxpar_
int GETMAXPAR();

#define GETMAXPSF getmaxpsf_
int GETMAXPSF();



#ifdef __cplusplus
}
#endif


#endif /* CDAOPHOT__H */
