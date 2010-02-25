
//#include "defs.h"
#include <stdlib.h>
#include <stdio.h>
//#include <string.h>
//#include <sys/time.h>
#include <math.h>
//#include "machine.h"
//#include "nbmath.h"
#include "nbrandom.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846 
#endif

#define DeuxPi (2.*M_PI)

#define drand01() ( (rand()*1.)/(RAND_MAX*1.0))


#define frand01() ( (float) drand01() )
#define ranf01() drand01()

#define frandpm1() ( 2. * frand01() - 1.)
#define drandpm1() ( 2. * drand01() - 1.)
#define ranfpm1() drandpm1()

#ifdef IS_IT_USEFUL

struct tirage_alea {
  int Nbin;
  double Min,Max,Lbin;
  double *Tab;
};
typedef struct tirage_alea TIREALEA;
#endif




//static double GAU_RANGE=6.;

/* 
++ 
  Module 	Tirages aleatoires (C)
  Lib	LibsUtil
  include	nbrandom.h
--
*/

/*
++
  frand01()
	tirage aleatoire entre 0 et 1, retourne float  
  drand01()
	tirage aleatoire entre 0 et 1, retourne double  
  rand01()
	 c'est le defaut: drand01()
  frandpm1()
	tirage aleatoire entre -1 et 1, retourne float
  drandpm1()
	tirage aleatoire entre -1 et 1, retourne double
  ranfpm1()
	c'est le defaut: drandpm1()
--
*/

/*=========================================================================*/
/*
++
  void Ini_Ranf_Quick(long seed_val, int lp)
	Initialisation rapide du generateur (drand48) par un entier
	de 32 bits de type long (cf srand48).
--
*/

/*commentee lors du passage a rand au lieu de rand48*/

/*void Ini_Ranf_Quick(long seed_val, int lp)
{
if(lp) printf("Ini_Ranf_Quick: %d\n",(int) seed_val);
srand48(seed_val);
return;
}*/

/*=========================================================================*/
/*
++
  void Ini_Ranf(unsigned short seed_16v[3], int lp)
	Initialisation complete du generateur (drand48) par
	48 bits (cf seed48).
--
*/
void Ini_Ranf(unsigned short seed_16v[3], int lp)
{
if(lp) printf("Ini_Ranf: %d %d %d\n"
             ,seed_16v[0],seed_16v[1],seed_16v[2]);
/*seed48(seed_16v);*/
srand(seed_16v[0]);
return;
}


float NorRand(void)
{
double x,A,B;

LAB10:
A = drand01();
if ( A == 0. ) goto LAB10;
B = drand01();
x = sqrt(-2.*log(A))*cos(DeuxPi*B);
return( (float) x );
}


/*=========================================================================*/
/*
++
  double GauRnd(double am, double s)
	Generation aleatoire gaussienne de centre "am" et de sigma "s".
--
*/
double GauRnd(double am, double s)
{
double x,A,B;

LAB10:
A = drand01();
if ( A == 0. ) goto LAB10;
B = drand01();
x = am + s * sqrt(-2.*log(A))*cos(DeuxPi*B);
return(x);
}


/*==========================================================================*/
/*
++
  void NormGau(double *x,double *y,double mx,double my,double sa,double sb,double teta);
	Tirage de 2 nombres aleatoires x et y distribues sur une gaussienne 2D
	de centre (x=mx,y=my), de sigmas grand axe et petit axe (sa,sb)
	et dont le grand axe fait un angle teta (radian) avec l'axe des x.
--
*/
/*
++
| - La densite de probabilite (normalisee a 1) sur laquelle on tire est:
| N*exp[-0.5*{ (A/sa)**2+(C/sc)**2 }],  N=1/(2Pi*sa*sc)
| ou A et B sont les coordonnees selon le grand axe et le petit axe
| et teta = angle(x,A), le resultat subit ensuite une rotation d'angle teta.
| - La matrice des covariances C des variables A,B est:
|   | sa^2   0   |
|   |            |  et det(C) = (1-ro^2)*sa^2*sb^2
|   |  0    sb^2 |
| - La distribution x,y resultante est:
| N*exp[-0.5*{[(dx/sx)^2-2*ro/(sx*sy)*dx*dy+(dy/sy)^2]/(1-ro^2)}]
| ou N est donne dans NormCo et sx,sy,ro sont calcules a partir
| de sa,sc,teta (voir fonctions paramga ou gaparam). La matrice des
| covariances des variables x,y est donnee dans la fonction NormCo.
--
*/
void NormGau(double *x,double *y
            ,double mx,double my,double sa,double sb,double teta)
{
double c,s,X,Y;

LAB10:
 s = drand01();
 if ( s == 0. ) goto LAB10;
s = sqrt(-2.*log(s));
c = DeuxPi * drand01();

X = sa*s*cos(c);
Y = sb*s*sin(c);

c = cos(teta); s = sin(teta);
*x = mx + c*X - s*Y;
*y = my + s*X + c*Y;
}

/*==========================================================================*/
/*
++
  int NormCo(double *x,double *y,double mx,double my,double sx,double sy,double ro)
	Tirage de 2 nombres aleatoires x et y distribues sur une gaussienne 2D
	de centre (mx,my), de coefficient de correlation rho (ro) et telle que
	les sigmas finals des variables x et y soient sx,sy (ce sont
	les valeurs des distributions marginales des variables aleatoires x et y
	c'est a dire les sigmas des projections x et y de l'histogramme 2D
	de la gaussienne). Retourne 0 si ok.
--
*/
/*
++
| - La densite de probabilite (normalisee a 1) sur laquelle on tire est:
|   N*exp[-0.5*{[(dx/sx)^2-2*ro/(sx*sy)*dx*dy+(dy/sy)^2]/(1-ro^2)}]
|     avec dx = x-mx, dy = y-my et N = 1/[2Pi*sx*sy*sqrt(1-ro^2)]
| - Dans ce cas la distribution marginale est (ex en X):
|   1/(sqrt(2Pi)*sx) * exp[-0.5*{dx^2/sx^2}]
| - La matrice des covariances C des variables x,y est:
|   |   sx^2      ro*sx*sy |
|   |                      |  et det(C) = (1-ro^2)*sx^2*sy^2
|   | ro*sx*sy      sy^2   |
| - La matrice inverse C^(-1) est:
|   |   1/sx^2      -ro/(sx*sy) |
|   |                           | * 1/(1-ro^2)
|   | -ro/(sx*sy)      1/sy^2   |
--
*/
/*
++
| - Remarque:
| le sigma que l'on obtient quand on fait une coupe de la gaussienne 2D
| en y=0 (ou x=0) est: SX0(y=0) = sx*sqrt(1-ro^2) different de sx
|                      SY0(x=0) = sy*sqrt(1-ro^2) different de sy
| La distribution qui correspond a des sigmas SX0,SY0
| pour les coupes en y=0,x=0 de la gaussienne 2D serait:
|   N*exp[-0.5*{ (dx/SX0)^2-2*ro/(SX0*SY0)*dx*dy+(dy/SY0)^2 }]
| avec N = sqrt(1-ro^2)/(2Pi*SX0*SY0) et les variances
| des variables x,y sont toujours
|  sx=SX0/sqrt(1-ro^2), sy=SY0/sqrt(1-ro^2)
--
*/
int NormCo(double *x,double *y
          ,double mx,double my,double sx,double sy,double ro)
{
double a,b,sa;
if( ro <= -1. || ro >= 1. ) return(1);
LAB10:
 b = drand01();
 if ( b == 0. ) goto LAB10;
b = sqrt(-2.*log(b));
a = DeuxPi * drand01();
sa = sin(a);

*x = mx + sx*b*(sqrt(1.-ro*ro)*cos(a)+ro*sa);
*y = my + sy*b*sa;

return(0);
}

/*==========================================================================*/
/*
++
  Titre	Exemple d'utilisation des aleatoires avec initialisation.
--
*/
/*
++
|   #include "nbrandom.h"
|   
|   void main() {
|   long i,ini=123456789;
|   unsigned short seed[3];
|   
|   printf(" 1./ ==> test nitialisation par un long\n");
|   Ini_Ranf_Quick(ini,1);
|   for(i=0;i<10;i++) printf("%d  -> %f\n",i,ranf01());
--
*/
/*
++
|   
|   printf("\n 2./ ==> test initialisation par tableau de 3 unsigned short\n");
|   Ini_Ranf_Quick(ini,1);
|   for(i=0;i<5;i++) printf("%d  -> %f\n",i,ranf01());
|   Get_Ranf(seed,1);
|   for(i=5;i<10;i++) printf("%d  -> %f\n",i,ranf01());
|   Ini_Ranf(seed,1);
|   for(i=5;i<10;i++) printf("%d  -> %f\n",i,ranf01());
|   Get_Ranf(seed,1);
--
*/
/*
++
|   
|   printf("\n 3./ ==> test initialisation automatique\n");
|   Auto_Ini_Ranf(2);
|   for(i=0;i<5;i++) printf("%d  -> %f\n",i,ranf01());
|   i=0; while(i<10000000) i++;
|   Auto_Ini_Ranf(2);
|   for(i=0;i<5;i++) printf("%d  -> %f\n",i,ranf01());
|   i=0; while(i<10000000) i++;
|   Auto_Ini_Ranf(2);
|   for(i=0;i<5;i++) printf("%d  -> %f\n",i,ranf01());
|   }
--
*/
/*
++
|    1./ ==> test initialisation par un long
|   Ini_Ranf_Quick: 123456789
|   0  -> 0.052468
|   1  -> 0.025444
|   2  -> 0.099272
|   3  -> 0.436130
|   4  -> 0.327740
|   5  -> 0.821202
|   6  -> 0.560493
|   7  -> 0.018157
|   8  -> 0.872758
|   9  -> 0.652496
--
*/
/*
++
|   
|    2./ ==> test initialisation par tableau de 3 unsigned short
|   Ini_Ranf_Quick: 123456789
|   0  -> 0.052468
|   1  -> 0.025444
|   2  -> 0.099272
|   3  -> 0.436130
|   4  -> 0.327740
--
*/
/*
++
|   Get_Ranf: 36117 51106 21478
|   5  -> 0.821202
|   6  -> 0.560493
|   7  -> 0.018157
|   8  -> 0.872758
|   9  -> 0.652496
--
*/
/*
++
|   Ini_Ranf: 36117 51106 21478
|   5  -> 0.821202
|   6  -> 0.560493
|   7  -> 0.018157
|   8  -> 0.872758
|   9  -> 0.652496
|   Get_Ranf: 16576 62373 42761
--
*/
/*
++
|   
|    3./ ==> test initialisation automatique
|   Auto_Ini_Ranf: date 887117206 s 868138 10^-6 sec seed=826006868:
|   ... njours=10267 nj23=9 buf=826006868.13800001
|   Ini_Ranf_Quick: 826006868
|   0  -> 0.798860
|   1  -> 0.342478
|   2  -> 0.401300
|   3  -> 0.442912
|   4  -> 0.170912
--
*/
/*
++
|   Auto_Ini_Ranf: date 887117207 s 188779 10^-6 sec seed=826007188:
|   ... njours=10267 nj23=9 buf=826007188.77900004
|   Ini_Ranf_Quick: 826007188
|   0  -> 0.455599
|   1  -> 0.811427
|   2  -> 0.703880
|   3  -> 0.409569
|   4  -> 0.390399
--
*/
/*
++
|   Auto_Ini_Ranf: date 887117207 s 489750 10^-6 sec seed=826007489:
|   ... njours=10267 nj23=9 buf=826007489.75
|   Ini_Ranf_Quick: 826007489
|   0  -> 0.567094
|   1  -> 0.893156
|   2  -> 0.975995
|   3  -> 0.531331
|   4  -> 0.834354
--
*/
/*==========================================================================*/

/*==========================================================================*/
/* 
++ 
  Module 	Tirages aleatoires selon une fonction (C)
  Lib	LibsUtil
  include	nbrandom.h
--
*/
/*
++
  TIREALEA *init_tirage_alea(int nbin,double xmin,double xmax,double (*fonc) (double))
	Initialise la structure qui va permettre le tirage aleatoire
	d'un nombre compris entre xmin et xmax selon la
	distribution fonc (histo de nbin bins)
--
*/

#ifdef IS_IT_USEFUL

TIREALEA *init_tirage_alea(int nbin,double xmin,double xmax,double (*fonc) (double))
{
int sof,i;
double x;
struct tirage_alea *t;

if ( xmax-xmin<0.) return(NULL);

if(nbin<=3) nbin=50;

sof = sizeof(struct tirage_alea);
if( (t = malloc(sof) ) == NULL ) {
  printf("impossible d'allouer *tirage_alea par malloc \n");
  return(NULL);
}

t->Nbin=nbin; t->Min=xmin; t->Max=xmax; t->Lbin=(xmax-xmin) /nbin;

sof = nbin * sizeof(double);
if( (t->Tab = malloc(sof) ) == NULL ) {
  printf("impossible d'allouer *tirage_alea.Tab par malloc \n");
  return(NULL);
}

x = xmin + .5*t->Lbin;
t->Tab[0] =  fonc(x);
for(i=1;i<nbin;i++) {
  x = xmin + (i+.5)*t->Lbin;
  t->Tab[i] = t->Tab[i-1] + fonc(x);
}

for(i=0;i<nbin-1;i++)  t->Tab[i] /= t->Tab[nbin-1];
t->Tab[nbin-1] = 1.;

return(t);
}

/*==========================================================================*/
/*
++
  double tirage_alea( TIREALEA *alea )
	tirage aleatoire d'un nombre compris entre xmin et xmax
	selon la fonction fonc (cf init_tirage_alea).
--
*/
double tirage_alea( TIREALEA *alea )
{ 
int i,ibin = -1;
double z,t1,t2,x1,x2,t;

z=drand01();
/* protections z<=0 ou z>=1 */
if( z <= 0. ) return ( alea->Min );
if( z >= 1. ) return ( alea->Max );
/* cas z <= tab[0] */
if(z <= alea->Tab[0]) {
  t = alea->Min + (alea->Lbin/2.)/alea->Tab[0] * z;
  return (t);
}

/* recherche du premier bin plus grand que z */
for(i=0;i<alea->Nbin;i++) {
  ibin=i;
  if ( z < alea->Tab[i] ) break;
}

/* extrapolation pour trouver la valeur du tirage aleatoire */
if( ibin == alea->Nbin-1 ) ibin--;
t1=alea->Tab[ibin];
x1 = alea->Min + (ibin+0.5) * alea->Lbin;
t2=alea->Tab[ibin+1];
x2 = x1 + alea->Lbin;
t = x1 + (x2-x1)/(t2-t1) *(z-t1);
if ( t < alea->Min ) t = alea->Min;
if ( t > alea->Max ) t = alea->Max;
return(t);
}

/*==========================================================================*/
/*
++
  int end_tirage_alea( TIREALEA *alea )
	De-allocation de la structure qui a permis le tirage aleatoire.
--
*/
int end_tirage_alea( TIREALEA *alea )
{ 
if ( alea != NULL ) { free(alea); return(0);}
  else return(-1);
}

#endif
