#ifndef ELLIPTICAPER__H
#define ELLIPTICAPER__H

#include <iostream>

#include "fastifstream.h"
#include <string.h>
//#include "persistence.h"
using namespace std;

class Image;
class SEStar;
class AperSEStar;

// definition elliptic aperture
class Elliptic_Aperture
{
public:
    double xc, yc;
    double a,b, angle; // redite parametre SE
    double cxx, cyy,cxy ; // redite de parametres SE qu'on ne recupere pas.
    double kron_factor; // in A and B units
    double fixradius; // in pixels
    double background ;
    double flux; // in ADUS
    double eflux;
    int nbad;
    int ncos; // number of bad pixels flagged as cosmics
    double fcos; // total flux of these pixels
    int  is_circle; //-1=non, 1=oui
    int is_good ;  //-1=non, 1=oui

    bool IsCircle() const {return(is_circle>0); }
    bool IsGood(){return(is_good>0); }

    Elliptic_Aperture() { xc=0. ; yc=0. ; a=0; b=0; angle=0;cxx =0 ; cyy = 0 ; cxy=0; kron_factor=0; fixradius=-1;background = 0. ;flux=0; eflux=0; nbad=0; ncos=0; fcos=0;is_circle=-1; is_good=1;  }

    void dump(ostream & pr){pr << " xc : " << xc << " yc : " << yc << " a: " << a << " b: " << b << " angle: " << angle <<  " cxx : " << cxx << " cyy : " << cyy << " cxy : " << cxy << " kron_factor : " << kron_factor << " fix radius : " << fixradius << " bck : " << background << " flux : " << flux << " eflux : " << eflux << " nbad : " << nbad << "is_circle : " << is_circle << "is_good : " << is_good << endl;}

  

    void WriteHeader_(ostream & pr, const char* i = NULL) const ;

 
    void writen(ostream & pr) const;


    void read(fastifstream& r );


    void computeflux(const Image& I, const Image& W, const Image *pC, const Image *pS, const double Gain, int segmentation_number, double scale_fact=1);

    //ovoid computeflux(const Image& I, const Image& W, const Image *pC, const double Gain);
    // RadMin et Radius = PHOT_AUTOAPER_[1] 
    void SetParameters(const SEStar & star,const double dilatation, 
		       const double RadMin, const double Radius, bool fromabtheta=false);
    void SetParameters_Aper(const AperSEStar & star,
			    const double kron_scale_factor, //2.5
			    const double kron_radius_min, //3.5
					   const double RadMin, const double Radius);
    double SqEllipticDistance(double x1, double y1 ) const ;
    double NormalizedEllipticDistance(double x1, double y1) const ;

   double SqEllipticDistance(double x_or, double y_or, double x1, double y1 ) const ;


  };



#endif /* ELLIPTICAPER__H */
