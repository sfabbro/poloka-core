#ifdef IS_IT_USEFUL
#ifndef DEFS_SEEN
#define DEFS_SEEN

/*#define DEBUG*/

#ifdef DEBUG
#define IMGRGCHECK
#define IMGVOIDCHECK
#endif

/* Sur mac, PP utilise fp.h plutot que math.h, et les deux sont incompatibles... */
#ifdef __MWERKS__
//#include <fp.h>
//#define __MATH__
#endif


/********************************************************/
/* Outil d'impressions pour debug */
#define PR(_data_)  {cout<<" "<<#_data_<<" "<<_data_;}
#define CR {cout<<"\n"<<flush;}
/* ********************************************************/

/* Sur quelle machine sommes-nous ?                                       */
/* Il faudrait uniformiser tous les sources avec les memes conventions... */

#if defined(__alpha) || defined(__alpha__) || defined(OSF1) || defined(DECALPHA)
#undef OSF1
#define OSF1
#undef DECALPHA
#define DECALPHA 
#undef __alpha__
#define __alpha__
#endif

#if defined(_IBMR2) && !defined(__GNUC__)
#define __xlC
#endif

#if defined(_AIX) && !defined(AIX)
#define AIX 1
#endif

#if defined(__hpux__) && !defined(HPUX)
#define HPUX  1
#endif


#ifdef ultrix
#define DECULTRIX 1
#endif

#ifdef __ultrix__
#define DECULTRIX 1
#endif

#ifdef THINK_CPLUS
#define __mac__
#endif

#ifdef __MWERKS__
#define __mac__
#endif

/* Quelques variantes du C++ selon le compilateur */

#define USESTRING

#define HAS_VEC_NEW

#ifdef __MWERKS__
#define __ANSI_TEMPLATES__
#undef HAS_VEC_NEW
#define COMPILER_EXCEPTIONS
#define NO_STRSTREAM
#define STREAMPOS_IS_CLASS
#endif

#ifdef __GNUG__
#define __GNU_TEMPLATES__
#endif

#ifdef __DECCXX
#if __DECCXX_VER > 60000000
#define CXX6 1
#define ITER_TAG
/*#define __GNU_TEMPLATES__*/
#define __CXX_PRAGMA_TEMPLATES__
#define COMPILER_EXCEPTIONS
#else
#define CXX5 1
#define __CXX_PRAGMA_TEMPLATES__
#endif

#endif

/* GCC 2.8.0 : exception, string::npos */
#ifdef __GNUC__
#if __GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ >= 8)
#define COMPILER_EXCEPTIONS
#define NPOS string::npos
#endif
#endif

#if defined(__cplusplus) && defined(AIX) && !defined(__GNUG__)
#define __xlC
#undef  COMPILER_EXCEPTIONS
#define COMPILER_EXCEPTIONS
#endif

#ifndef __GNUG__
/*typedef int bool;*/
#define __PRETTY_FUNCTION__ __FILE__ " __LINE__ "
#ifndef __FUNCTION__
#define __FUNCTION__ __FILE__ " __LINE__ "
#endif
#define __attribute__(_x_)
#endif

#ifdef __GNUG__
#define HAS_NAMED_RETURN 0
#else
#define HAS_NAMED_RETURN 0
#endif

#ifdef THINK_CPLUS
#define ARG_FUNC_CPP_C 0
#else
#define ARG_FUNC_CPP_C 1
#endif

#ifdef __GNUG__
#define HAS_EXPLICIT
#endif

#ifdef HAS_EXPLICIT
#define EXPLICIT explicit
#else
#define EXPLICIT
#endif
/* Quelques variantes de la librairie C++.                 */
/* La librairie est AT&T par defaut, sauf si __ANSI_IO__   */

/* $CHECK$ EA : a uniformiser */
/* remplacer HAS_IOS_BIN par un selecteur qui indique si c'est ios::binary ou */
/* ios::bin */

#ifdef __MWERKS__
//#define __ANSI_IO__
#endif

#if defined(THINK_CPLUS)
#define HAS_IOS_BIN 0
#define HAS_IOS_NOCREATE 1
#define FITS_IOS_IN_OPT  ios::in | ios::nocreate
#define FITS_IOS_OUT_OPT ios::out

#elif defined(__MWERKS__)
#define HAS_IOS_BIN 1
#define HAS_IOS_NOCREATE 1
#define NPOS string::npos

#elif defined(__GNUG__)
#define HAS_IOS_BIN 1
#define HAS_IOS_NOCREATE 1
#define FITS_IOS_IN_OPT  ios::in | ios::bin | ios::nocreate
#define FITS_IOS_OUT_OPT ios::out | ios::bin

#else
#define HAS_IOS_BIN 0
#define HAS_IOS_NOCREATE 0
#endif

#ifdef __ANSI_IO__
#define seekg(p_, m_) rdbuf()->pubseekoff(p_,m_)
#define seekp(p_, m_) rdbuf()->pubseekoff(p_,m_)
#endif


#define STR2CH(x) ((string(x)).c_str())

/* Des fonctions qui manquent sur certaines machines  */

#if defined(__DECCXX) 
/* a cause des exceptions dans math.h */
#define exception math_exception
#include <math.h>
#undef exception
#else
#include <math.h>
#endif

#if defined(__DECCXX) 
/*  Pour definir bool   Reza 20/05/97  */
#include <stl_macros>
#endif

#if defined(__xlC)  && defined (__cplusplus)
#include <bool.h>
#endif

#if defined( __DECCXX ) || defined (__aCC__)
/* en fait si c'est RogueWave */
/*  Remis par Reza  20/05/97  */
#define NPOS (size_t)-1
#endif


#ifndef M_PI
#define M_PI 3.1415926535
#endif

#ifdef __mac__
#define hypot(_x_,_y_) sqrt((_x_)*(_x_) + (_y_)*(_y_))
//#define nintf(x) ((int) (x))
#define random() (rand()*65538.0 + rand()*2. + rand())
#define srandom srand
#define initstate(seed, tab, n) srand(seed)
#define nice(_x)
/*#define getpid() 1 */

#ifdef __cplusplus
#include <string.h>
#include <string>
inline bool operator< (string const& s1, string const& s2)
{ return (strcmp(s1.c_str(), s2.c_str()) < 0); }
#endif
#endif



#endif

#endif /* IS_IT_USEFUL */
