#ifndef  DODETECTION__H 
#define  DODETECTION__H

#include "image.h"
#include "candidatestar.h"
#include "candstar.h"
double
Noise_Aperture(Image  const & img, double x, double y, double radius, 
	       double Fond, double SigFond);
int 
NewCvDetection( Image & img , Image   & mask,
		CandidateStarList & stl,  DatDetec & datdet, 
		double & seeing, double & sigweigth);

bool
IsNotBadPix( Image const & img , Image const & mask, int x , int y , 
	     float  rad_bad);

void 
Flux_Coord(Image const & img , double x, double y , double rad_flux,
	   float fond, double & xbar, double & ybar, double & flux );

double Flux_Aperture(const Image   & img, const double x, const double y, const double radius, const double Fond, int &Npix);

// NRL: not implemented !
//double
//FondLocal( Image const & image ,float xcentre , float ycentre,
//           float rayon1 , float rayon2 , float fraclow,
//           float frachigh , int lp);

// NRL: not implemented !
//void
//Sig_Position(Image const & img ,double const &rad_flux, 
//		      double & xbar, double & ybar, const double & fond, 
//		      double & Sigx, double & Sigy, double & Sigxy);
     
void Barycentre(Image const & img , double x, double y , double rad_flux,
		double fond, double & xbar, double & ybar, double & flux );
#endif
