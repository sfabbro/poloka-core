#include <stdlib.h> // for random
#ifdef __sun
extern "C" {
extern long random(void);
}
#endif
 
#include <unistd.h>
#include <sys/time.h> 


#include "image.h"

void RandomSeed()
{
  struct timeval now;
  long nj,nj23;


 unsigned int seed=0;

 double buf;

 gettimeofday (&now,0);

 /* dans 32 bits signes on met environ 23 jours a 1/1000 de seconde pres! */
 /* Nombre de jours depuis l'origine */
 nj = (long) now.tv_sec / 86400;
 /* Nombre de jours depuis le dernier jour multiple de 23 jours */
 nj23 = nj % 23;
 /* nombre de secondes depuis le dernier jour multiple de 23 jours */
 buf = (double) (nj23*86400 + (now.tv_sec-nj*86400));
 /* nombre de milliemes de secondes depuis ... */
 buf = buf*1000. +  now.tv_usec/1000.;

 seed = (unsigned int) buf;
 srandom(seed);
 return;
}

double my_rndm()
{
  return ((1.0*random())/(RAND_MAX+1.0));
}

// tirage de x selon une loi gaussienne normee centree
// re-pompee de EROS
float NormalGaussRand()
{
double x,A,B;
A = my_rndm();
while ( A == 0. ) A = my_rndm() ;
B = my_rndm();
x = sqrt(-2.*log(A))*cos(2. * M_PI *B);
return( (float) x );
}

// ajouter a un pixel un bruit gaussien de sigma^2=valeur du pixel
// en tenant compte du gain et du read-out noise
float Gaussian_Noisy(float pixel, float Gain, float ROnoise)
{
  double ronsq = ROnoise/Gain;  ronsq *= ronsq;
  double sigma = ronsq+fabs(pixel)/Gain;
  sigma = sqrt(sigma);
  return(pixel + sigma * NormalGaussRand());
}

// tirage de x et y selon une loi gaussienne centree
// en (0,0) de parametres sig_a, sig_b, theta
// re-pompee de EROS
void NormGauss(double & x, double & y, double sig_a, double sig_b, double theta)
{
  double s = my_rndm();
  while ( s < 1.e-20 ) s = my_rndm();
  s = sqrt(-2.*log(s));
  double c = 2* M_PI * my_rndm();

  double X = sig_a*s*cos(c);
  double Y = sig_b*s*sin(c);

  c = cos(theta); s = sin(theta);
  x = c*X - s*Y;
  y = s*X + c*Y;
}


void ImgAddNoise(Image & img,float Gain, float ROnoise)
{
  int size = img.Nx() * img.Ny() ;
  float * p = img.begin();
  for (int i=size ; i>0 ; i--)
    {
      *p = Gaussian_Noisy(*p,Gain,ROnoise) ;
      ++p;
    }
}


