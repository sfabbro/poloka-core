#ifndef CDAOPHOT_H
#define CDAOPHOT_H

#ifdef __cplusplus
extern "C" {
#endif

/*! Get sky mean, median, mode of an image */
#define GETSKY getsky_
void GETSKY(float *d, float *s, int *index, int *max, const float *readns, const float *hibad, 
	    float *skymn, float *skymed, float *skymod, float *skysig, int *n);

/*! Mean Median Mode DAOPHOT routine, for sky level and such */
#define MMM mmm_
void MMM(const float *sky, const int *nsky, const float *hibad, const float *readns, 
	 float *skymn, float *skymed, float *skymod, float *sigma, float *skew);

/*! read the PSF DAOPHOT file */
#define RDPSF rdpsf_
int RDPSF(const char *psffil, int *ipstyp, float *par, const int *maxpar, int *npar, 
	  float *psf, const int *maxpsf, const int *maxexp, int *npsf, int *nexp, int *nfrac, 
	  float *psfmag, float *bright, float *xpsf, float *ypsf);
    
/*! returns PSF value on a point distant of (dx,dy) of a star located on (deltax,deltay) */
#define USEPSF usepsf_
float USEPSF(const int *ipstyp, const float *dx, const float *dy,  const float *bright, 
	     const float *par, const float *psf, const int *npsf, const int *npar, 
	     const int *nexp, const int *nfrac, const float *deltax, const float *deltay, 
	     float *dvdxc, float *dvdyc);

/*! This is the subroutine which does the actual one-star least-squares profile fit for DAOPK */
#define PKFIT pkfit_
void PKFIT(const float *f, const int *nx, const int *ny, const int *maxbox, 
	   float *x, float *y, float *scale, const float *sky, const float *radius, 
	   const float *lobad, const float *hibad, const float *phpadu, 
	   const float *ronois, const float *perr, const float *pkerr, 
	   const float *bright, const int *ipstyp, const float *par, 
	   const int *maxpar, const int *npar, const float *psf, 
	   const int *maxpsf, const int *maxexp, const int *npsf, const int *nexp, 
	   const int *nfrac, float *deltax, float *deltay, float *errmag, float *chi, 
	   float *sharp, int *niter, const float *global_sky);

/*! set the options */
#define OPTION option_
void OPTION(const char *optfile, const int *nopt, const char *names,  float *opt, 
	    const float *omin, const float *omax, const char *prompt, int *istat);


/*! find stars with gaussian filter */
#define FIND find_
void FIND(const float *d, float *h, int *jcyln, float *g, short int *skip, const int *max, 
	  const int *maxbox, const int *maxcol, const int *maxsky, const float* opt, 
	  const int *nopt, const char* coofil, const float *global_sky);


/*! performs aperture photometry with daophot */
#define PHOTSB photsb_
void PHOTSB (const float *f, float *sky, int *index, const int *maxsky, const int *ncol, 
	     const int *nrow, const float *inbad, const float *watch, const float *estmtr,
	     const char* coofil, const char* magfil, const char* table, const float *global_sky);

     

/*! pick stars to make a psf */
#define PCKPSF pckpsf_
void PCKPSF  (int *id, float *x, float *y, float *m , float *s, int *index, const int *maxn, 
	      const float *fitrad, const float *psfrad, const float *vary, const char* magfil,
	      const char* strfil, const int *nreq, const float* maglim);
  
/*! get a psf from an image */
#define GETPSF getpsf_
void GETPSF (const float *array, const int *nx, const int *ny, float *par, float *psf, 
	     const float *opt, const int *nopt, const char* magfil, const char* strfil, const char* psffil, 
	     const char* neifil, const float *global_sky);

/*!  Single-star profile-fitting routine on every stars of an image */
#define DAOPK daopk_
void DAOPK(const float *par, const int *maxpar, const float *psf, const int *maxpsf, 
	   const int *maxexp, const float *f, const float *watch, const float *fitrad, 
	   const float *pererr, const float *proerr, const char* magfil, const char* psffil,
	   const char* pkfil, const float *global_sky);

/*! group stars for multi-star NSTAR */
#define GROUP group_
void  GROUP(const float *par, const int *maxpar, const float *psf, const int *maxpsf, const int *maxexp,
	    const int *maxstr, const float *fitrad, const float *psfrad, const char* magfil, const char* psffil, 
	    const char* grpfil, const float* crit, const float *global_sky);

/*!  Multi-star profile-fitting routine on every stars of an image */
#define NSTAR nstar_
void NSTAR(const float *par, const int *maxpar, const float *psf, const int *maxpsf, const int *maxexp, 
	   const float *data, const int *ncol, const int *nrow, const float *watch, const float *fitrad, 
	   const float *pererr, const float *proerr, const char* grpfil, const char* psffil, const char* fitfil,
	   const float *global_sky);

/*! subtract stars from an image using a scaled PSF */
#define SUBSTR substr_
void SUBSTR(const float *par, const int *maxpar, const float *psf, const int *maxpsf, 
	    const int *maxexp, float *f, const int   *ncol, const int *nrow, const float *watch,
	    const char* psffil, const char* profil, const char* lstfil, const char* subpic);

/*! generates fake stars from PSF and add poisson noise at locations given in a file */
#define ADDSTR addstr_
void ADDSTR(const float *par, const int *maxpar, const float *psf, const int *maxpsf, 
	    const int *maxexp, float *f, const int *ncol, const int *nrow, const float *watch,
	    const char* psffil, const char* addpic, const char* outstm, const float *rmag, const int *nstar,
	    const int *nframe, const int *inseed);

/*! fit all stars group by group with a PSF on an image */
#define ALLSTR allstr_
void ALLSTR(const float *array, const int *nx, const int *ny, float *subt, float *sigma, 
	    const float *fitrad, const float *watch, const float *half, const int *iexp, 
	    const int *center, const int *maxgrp, const float *pererr, const float *proerr, 
	    const float *skyin, const float *skyout, int *istat, const char* psffil, const char* magfil,
	    const char* alsfil, const char* subpic, const float *global_sky);

/*! append two daophot star list file */
#define APPEND append_
void APPEND(const char* ifile1, const char* ifile2, const char* cmbfile);

/*! dump intensities values on screen around a point */
#define DUMP dump_
void DUMP(const float *f, const int *ncol, const int *nrow, const float *coords, const float *size);

/*! fudge pixel intensity values */
#define FUDGE fudge_
void FUDGE(const char* file, float *f, const char* newfil, const int* nbord, 
	   const int *npoly, const int *lx, const int* mx, const int* ly, const int* my, 
	   const float* bright);

/*! select star groups */
#define DAOSLT daoslt_
  void DAOSLT(const int* max, const char* ingrpfile, const char* outgrpfile, 
	    const int *mingrp, const int* maxgrp);

/*! offset a star list file with positions or magnitudes */
#define OFFSET offset_
  void OFFSET(const char* file, const float *delta, const char* offile);

/*! sort star file according to some rule*/
#define SORTER sorter_
void SORTER(const int *maxstr, const float *watch, const int* mode, const char* file, 
	    const char* srtfile, const char *answer);

#ifdef __cplusplus
}
#endif


#endif /* CDAOPHOT__H */
