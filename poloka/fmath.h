#ifdef IS_IT_USEFUL
/*    Declaration de fonctions math pour float   */

#ifndef  FMATH_H_SEEN
#define  FMATH_H_SEEN

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif


#ifdef OSF1
#define   fexp(x)       expf(x)
#define   fsqrt(x)      sqrtf(x)
#define   flog10(x)     log10f(x)
#define   flog(x)       logf(x)
#else
#ifdef ULTRIX
float     fexp(float x);            
#define   expf(x)  fexp(x)
float     fsqrt(float x);
float     flog10(float x);
float     flog(float x);
#define   powf(x,y)    ((float) pow((double)(x), (double)(y)))
#define   fabsf(x)    ((float)(fabs((double)(x))))
#define   ceilf(x)  fceil(x)
#define   nintf (int)
#else
#define   fexp(x)     ((float)(exp((double)(x))))
#define   expf(x)     ((float)(exp((double)(x))))
#define   fabsf(x)    ((float)(fabs((double)(x))))
#define   fsqrt(x)    ((float)(sqrt((double)(x))))
#define   sqrtf(x)    ((float)(sqrt((double)(x))))
#define   flog10(x)   ((float)(log10((double)(x))))
#define   log10f(x)   ((float)(log10((double)(x))))
#define   flog(x)     ((float)(log((double)(x))))
#define   logf(x)     ((float)(log((double)(x))))
#define   floorf(x)   ((float) floor((double)(x)))
#define   ceilf(x)    ((float) ceil((double)(x)))
#define   powf(x,y)    ((float) pow((double)(x), (double)(y)))
#define   nintf (int)
#endif
#endif

#ifdef __cplusplus
}
#endif

#endif

#endif /* IS_IT_USEFUL */
