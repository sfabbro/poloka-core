#ifndef CDAOPHOT_H
#define CDAOPHOT_H

#define FILNAM filnam_
#define SIZE size_
#define OPTION option_
#define MMM mmm_
#define GETSKY getsky_
#define PHOTSB photsb_
#define PCKPSF pckpsf_
#define GETPSF getpsf_
#define RDPSF rdpsf_
#define USEPSF usepsf_
#define DAOPK daopk_
#define PKFIT pkfit_
#define NSTAR nstar_
#define SUBSTAR substar_
#define ADDSTR addstr_
#define ALLSTR allstr_
#define ATTACH attach_
#define RDARAY rdaray_
#define WRARAY wraray_
#define COPPIC coppic_
#define CLPIC clpic_
#define DELPIC delpic_
#define LIST list_
#define SPLASH splash_

const int MAXPSF = 207; // maximum PSF array size allowed by DAOPHOT
const int MAXPAR = 6;   // maximum number of PSF parameters
const int MAXEXP = 10;  // maximum coeff for interpolating spatial variable table
const int MAXSKY = 10000;
#define NCHARFILE (int) 256

#ifdef __cplusplus
extern "C" {
#endif
  
struct Common_Filename
{
  char COOFIL[NCHARFILE];  /* Name of the object coordinate file. */
  char MAGFIL[NCHARFILE];  /* Name of the photometry file. */
  char PSFFIL[NCHARFILE];  /* Name of the psf description file. */
  char PROFIL[NCHARFILE];  /* Name of the psf fit file. */
  char GRPFIL[NCHARFILE];  /* Name of the object groups file. */
};
  
struct Common_Size
{
  int ncol;  /* The number of columns (ie. x axis range) in the image file. */
  int nrow;  /* The number of rows (ie. y axis type) in the image file. */
};

/* Define the instances of the FORTRAN COMMON blocks that DAOphot II and 
   Allstar use. */

extern struct Common_Filename FILNAM;
extern struct Common_Size SIZE;
  
//! Get sky mean, median, mode of an image
void GETSKY(float *d, float *s, int *index, int *max, float *readns, float *hibad, 
	    float *skymn, float *skymed, float *skymod, float *skysig, int *n);

//! set the options
void OPTION(const char *optfile, const int *nopt, const char *names,  float *opt, 
	    const float *omin, const float *omax, const char *prompt, int *istat);

//! Mean Median Mode DAOPHOT routine, for sky level and such
void MMM(const float *sky, const int *nsky, const float *hibad, const float *readns, 
	 float *skymn, float *skymed, float *skymod, float *sigma, float *skew);

//! performs aperture photometry with daophot
void PHOTSB (const float *f, float *sky, int *index, const int *maxsky, const int *ncol, 
	     const int *nrow, const float *inbad, const float *watch, float *estmtr,
	     const float *global_sky);

     

//! pick stars to make a psf
void PCKPSF  (int *id, float *x, float *y, float *m , float *s, int *index, const int *maxn, 
	      const float *fitrad, const float *psfrad, const float *vary);
  
//! get a psf from an image
void GETPSF (const float *array, const int *nx, const int *ny, float *par, float *psf, 
	     const float *opt, const int *nopt, const float *global_sky);

//! read the PSF DAOPHOT file in calibrated.psf
int RDPSF(const char *psffil, int *ipstyp, float *par, const int *maxpar, int *npar, 
	  float *psf, const int *maxpsf, const int *maxexp, int *npsf, int *nexp, int *nfrac, 
	  float *psfmag, float *bright, float *xpsf, float *ypsf);
    
//! returns PSF value on a point distant of (dx,dy) of a star located on (deltax,deltay) 
float USEPSF(const int *ipstyp, const float *dx, const float *dy,  const float *bright, 
	     const float *par, const float *psf, const int *npsf, const int *npar, 
	     const int *nexp, const int *nfrac, const float *deltax, const float *deltay, 
	     float *dvdxc, float *dvdyc);

//!  Single-star profile-fitting routine on every stars of an image
void DAOPK(const float *par, const int *maxpar, const float *psf, const int *maxpsf, 
	   const int *maxexp, const float *f, const float *watch, const float *fitrad, 
	   const float *pererr, const float *proerr, const float *global_sky);

//! This is the subroutine which does the actual one-star least-squares profile fit for DAOPK
void PKFIT(const float *f, const int *nx, const int *ny, const int *maxbox, 
	   float *x, float *y, float *scale, const float *sky, const float *radius, 
	   const float *lobad, const float *hibad, const float *phpadu, 
	   const float *ronois, const float *perr, const float *pkerr, 
	   const float *bright, const int *ipstyp, const float *par, 
	   const int *maxpar, const int *npar, const float *psf, 
	   const int *maxpsf, const int *maxexp, const int *npsf, const int *nexp, 
	   const int *nfrac, float *deltax, float *deltay, float *errmag, float *chi, 
	   float *sharp, int *niter, const float *global_sky);

//! generates fake stars from PSF and add poisson noise at locations given in a file
void ADDSTR(const float *par, const int *maxpar, const float *psf, 
	    const int *maxpsf, const int *maxexp, float *f, 
	    const int *ncol, const int *nrow, const float *watch);

//! fit all stars group by group with a PSF on an image
void ALLSTR(const float *array, const int *nx, const int *ny, float *subt, float *sigma, 
	    const float *fitrad, const float *watch, const float *half, const int *iexp, 
	    const int *center, const int *maxgrp, const float *pererr, const float *proerr, 
	    const float *skyin, const float *skyout, int *istat, const float *global_sky);


/*! Opens the image file without reading it.*/
void ATTACH(const char* image, int* open);

/*! Reads pixels from the attached image (or the working copy of the image).
  If text is "DATA" then the pixels will be copied from the image, if
  text is "COPY" then the pixels will be copied from the working image copy. */
void RDARAY(const char* text, int *lx, int *ly, int *mx, 
	    int* my, const int* maxx, float* func, int* ier);

/*! Writes pixels to the attached image (or the working copy of the image).
  If text is "DATA" then the pixels will be copied to the image, or if
  text is "COPY" then the pixels will be copied to the working image copy. */
void WRARAY(const char* text, int *lx, int *ly, int *mx, 
            int* my, const int* maxx, float* func, int* ier);

void COPPIC(const char* picture, int* ier);

/*! Closes an image file. */
void CLPIC(const char* text);

/*! Deletes the image whose name is specified as image. */
void DELPIC(const char* image, int* ier);

/*! Give the name of the current file */
void LIST(const char* file);

/*! Flush the working copy image */
void SPLASH();

#ifdef __cplusplus
}
#endif


#endif // CDAOPHOT__H
