#ifndef ELLIPTICAPER__H
#define ELLIPTICAPER__H

#include <iostream>


//#include "persistence.h"
using namespace std;

class Image;
class SEStar;

// definition elliptic aperture
class Elliptic_Aperture
{
public:
    double xc, yc;
    double cxx, cyy,cxy ; // redite de parametres SE qu'on ne recupere pas.
    double kron_factor; // in A and B units
    double fixradius; // in pixels
    double background ;
    double flux; // in ADUS
    double eflux;
    int nbad;
    int ncos; // number of bad pixels flagged as cosmics
    double fcos; // total flux of these pixels
    bool is_circle;

    bool IsCircle(){return(is_circle); }

    Elliptic_Aperture() { xc=0. ; yc=0. ; cxx =0 ; cyy = 0 ; cxy=0; kron_factor=0; fixradius=-1;background = 0. ;flux=0; eflux=0; nbad=0; ncos=0; fcos=0;is_circle=false;  }
    void dump(ostream & pr){pr << " xc : " << xc << " yc : " << yc << " cxx : " << cxx << " cyy : " << cyy << " cxy : " << cxy << " kron_factor : " << kron_factor << " fix radius : " << fixradius << " bck : " << background << " flux : " << flux << " eflux : " << eflux << " nbad : " << nbad << endl;}
    void computeflux(const Image& I, const Image& W, const Image *pC, const Image *pS, const double Gain, int segmentation_number);
    //ovoid computeflux(const Image& I, const Image& W, const Image *pC, const double Gain);
    void SetParameters(const SEStar & star,const double dilatation, 
		       const double RadMin, const double Radius, bool fromabtheta=false);
    double SqEllipticDistance(double x1, double y1 );
   double SqEllipticDistance(double x_or, double y_or, double x1, double y1 );
  };



#endif /* ELLIPTICAPER__H */
