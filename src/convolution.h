#ifndef CONVOLUTION__H
#define CONVOLUTION__H
// This may look like C code, but it is really -*- C++ -*-
#include "image.h"

// procedures pour la convolution




// Pour faire le noyau de convolution. C'est interne en fait.
// On ne peut prendre qu'un noyau gaussienne integree pour l'instant.

Image
Convole_MakeNoyau(double r, int l) ;
double *  
Convole_MakeNoyau_1D(double r, int l) ;


Image
Convole_MakeNoyau(double rx, double ry, int l) ;


// procedures de convolution

// convolution de source, resultat dans a. Le mask sert a ne convoler que les pixels
// qui y sont mis a 1.
// l est la demi-taille de la fenetre de convolution
// r est le sigma de la PSF utilisee pour convoluer, 
// rx et ry sont le sigmax et sigmay (rho =0) de la PSF utilisee pour convoluer.
// convole1DX convolue dans la direction X  (PSF 1 D de sigma r),
// convole1DY convolue dans la direction Y  (PSF 1 D de sigma r), 
// convole convolue dans les 2 directions. 


void 
Convole(Image const& src, double r, int l , 
	Image const& mask, Image& a ) ;

void 
Convole(Image const& src, double r, int l , Image & a ) ;

void 
Convole(Image const& src, double rx,double ry, int l , Image & a ) ;
void 
Convole(Image const& src, double rx,double ry, int l , 
	Image const& mask, Image & a ) ; 

void 
Convole1DX(Image const& src, double r, int l , 
	   Image  const& mask, Image & a ) ;

void 
Convole1DX(Image const& src, double r, int l , Image & a ) ;

void 
Convole1DY(Image const& src, double r, int l , 
	   Image const& mask, Image& a ) ;

void 
Convole1DY(Image const& src, double r, int l , Image& a ) ;


int
Demi_Fenetre_Optimale(double sigma, double precision=1.e-5);
int
Demi_Fenetre(int l , double sigma, double precision=1.e-5);

//! Convolve Image "In" with a kernel along axes. Precision drives the Window size. returns the kernel
Image ConvoleGauss1D1D( const Image &In, const double SigKernX, 
		       const double SigKernY, const double Precision,
		       Image &Out);

//! Same as above but convole (correctly) "in place". returns the kernel.
Image ConvoleGauss1D1D( Image &InOut, const double SigKernX, 
		       const double SigKernY, const double Precision);
#endif /* CONVOLUTION__H */
