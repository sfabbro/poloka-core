/* quelques fonctions mathematiques et constantes utiles */
/*                                         cmv  23/06/94 */
#ifndef NBCONST_H_SEEN
#define NBCONST_H_SEEN


#ifdef __cplusplus
extern "C" {
#endif

// EA conflit entre les diverses copies de nbmath.h, nbconst.h... A nettoyer...
#undef  Pi   
#undef  Pis2
#undef  DeuxPi 
#undef  SPi   
#undef  S2Pi   
#undef  Rac2   
#undef  Log2   
#undef  LnPi   
#undef  LgPi   
#undef  Ln10   
#undef  DftoDm 
#undef  Hln2pi 
#undef  JourSec 
#undef  AnSec  

#undef  GRAND  
#undef  GRAND2 
#undef  IGRAND 

#define  Pi     (double) (3.1415926535897931)   /* c'est Pi */
#define  Pis2   (double) (1.57079632679489655)  /* c'est Pi/2 */
#define  DeuxPi (double) (6.2831853071795862)   /* c'est 2*Pi */
#define  SPi    (double) (1.7724538509055159)   /* c'est sqrt(Pi) */
#define  S2Pi   (double) (2.5066282746310002)   /* c'est sqrt(2*Pi) */
#define  Rac2   (double) (1.4142135623730950)   /* c'est sqrt(2) */
#define  Log2   (double) (0.30102999566398119)  /* c'est log10(2) */
#define  LnPi   (double) (1.1447298858494002)   /* c'est ln(Pi) */
#define  LgPi   (double) (0.49714987269413385)  /* c'est log10(Pi) */
#define  Ln10   (double) (2.3025850929940456)   /* c'est ln(10) */
#define  DftoDm (double) (1.0857362047581295)   /* c'est 2.5/ln(10) */
#define  Hln2pi (double) (0.91893853320467267)  /* c'est Ln(2*pi)/2 */
#define  JourSec (float) (86400.)    /* nombre de secondes dans 24H */
#define  AnSec  (int) (31557600)     /* nombre de secondes dans 365.25 jours */

#define  GRAND  (float) (1.e+35)
#define  GRAND2 (double) (1.e+35)
#define  IGRAND (int_4) (2147483647)

#ifdef __cplusplus
}
#endif

#endif
