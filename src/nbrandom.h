/* fonctions pour generateurs aleatoires */
/*                         cmv  23/06/94 */

#ifndef NBRANDOM_H_SEEN
#define NBRANDOM_H_SEEN

#ifdef IS_IT_USEFUL
#include "defs.h"
#include <stdlib.h>
#ifdef __mac__
#include <unixmac.h>
#endif
#endif /*IS_IT_USEFUL */

//#define drand01() ( drand48() )
#ifdef __cplusplus
extern "C" {
#endif

/*supprime lors du passage a rand
void Ini_Ranf_Quick(long seed_val, int lp);*/
void Ini_Ranf(unsigned short seed_16v[3], int lp);
    /* plus de Get_Ranf, car rand ne peremet pas de recuperer la seed
void Get_Ranf(unsigned short seed_16v[3], int lp);*/
  //void Auto_Ini_Ranf(int lp);

  //void SetGauRange(double range);
  // float NorRand1(void);
  // double GauRnd1(double am, double s);
float NorRand(void);
double GauRnd(double am, double s);
int NormCo(double *a,double *b
          ,double mx,double my,double sx,double sy,double ro);
void NormGau(double *x,double *y
            ,double mx,double my,double sa,double sb,double teta);

#ifdef IS_IT_USEFUL
TIREALEA *init_tirage_alea(int nbin,double xmin,double xmax,double (*fonc) (double));
double tirage_alea(TIREALEA *alea);
int end_tirage_alea(TIREALEA *alea);
#endif

#ifdef __cplusplus
}
#endif

#endif /* IS_IT_USEFUL */

