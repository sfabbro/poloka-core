#ifndef WCSCON__H
#define WCSCON__H


#define WCS_J2000	1	/* J2000(FK5) right ascension and declination */
#define WCS_B1950	2	/* B1950(FK4) right ascension and declination */

/* sys1 and sys2 to be defined using the 2 tags above (much more than those are 
indeed available). eq1 and eq2 may then be 0. theta and phi are ra and dec.
epoch should be zero for most if not all applications . */
/* this code was taken in WCSLIB  http://tdc-www.harvard.edu/software/wcstools/ */

#ifdef __cplusplus
extern "C" {
#endif

void
wcscon (int sys1, int sys2, double eq1, double eq2, double *dtheta, double *dphi, double epoch);

#ifdef __cplusplus
           }
#endif

#endif
