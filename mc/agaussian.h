// This may look like C code, but it is really -*- C++ -*-
#ifndef AGAUSSIAN__H
#define AGAUSSIAN__H

//#error 'coucou'

#include <stdlib.h>
#include <math.h>  



class Image;


//====================== quick gaussian simulation for debug
class  AGaussian {
public:
  double xc ;
  double yc ;
  double sigma_x ;
  double sigma_y ;
  double rho ;
  double flux ;
  double fond ;

  AGaussian();
  double Norma()const {return( flux * sqrt(1.-rho*rho)/(2*M_PI*sigma_x*sigma_y));}
  double ExpValue(double xin, double yin)const ;
  double ExpIntegValue(int i, int j) const ;
  double Value(double xin, double yin)const {return(Norma()*ExpValue(xin, yin)+fond);}
  double IntegValue(int i, int j)const {return(Norma()*ExpIntegValue(i, j)+fond);}
  void AddToImage(Image &image, Image * psat=NULL, double saturation=1) const;
  //private :
      double window_size() const ;  
};


#include "basestar.h"

void AddListWGaussianToImage(double sigmax, double sigmay, double rho,
			     BaseStarList *list, 
			     Image & dest,  
			     Image * psat=NULL, 
			     double satlevel=-1) ;










#endif
