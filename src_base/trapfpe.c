
/*

   a piece of code found on 
   http://gcc.gnu.org/onlinedocs/gcc-3.4.1/g77/Floating-point-Exception-Handling.html

   only added the getenv to allow "run time" selection of behavior 
*/


#include <stdlib.h> /* for getenv */

#define _GNU_SOURCE 1
#define __USE_GNU
#include <fenv.h>


static void __attribute__ ((constructor))
trapfpe ()
{
  /* Enable some exceptions.  At startup all exceptions are masked.  */
  
  if (getenv("DUMP_CORE_ON_FPE"))
    feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
     
