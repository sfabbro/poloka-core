// This may look like C code, but it is really -*- C++ -*-
#ifndef SIMFITVIGNET__H
#define SIMFITVIGNET__H

#include <countedref.h>
#include <daophotpsf.h>
#include <reducedimage.h>

#include "photstar.h"
#include "vignet.h"

//
//! \file simfitvignet.h
//! \brief a Vignet with its associated tabulated kernels, psfs, residuals.
//!  Primarly to use for SimultaneousFit.
//

class TabulatedPsf : public Kernel {

public:
  
  //! empty constructor allocate nothing
  TabulatedPsf() {}
  
  //! allocate psf and derivatives of half-sizes Hx and Hy
  TabulatedPsf(const int Hx, const int Hy) : Kernel(Hx,Hy), Dx(Hx,Hy), Dy(Hx,Hy)  {}
  
  //! allocate psf and derivatives of a Radius
  TabulatedPsf(const int Radius) : Kernel(Radius), Dx(Radius), Dy(Radius)  {}
  
  //! allocate psf and derivatives, and fill them with DaoPsf value on that Pt.
  TabulatedPsf(const Point& Pt, const DaoPsf& Psf, const int Radius);

  //! allocate psf and derivatives, and fill them with DaoPsf value on that Pt
  TabulatedPsf(const Point& Pt, const DaoPsf& Dao, const Window& Rect);

  // default destructor, copy constructor and assigning operator are OK  

  //! tabulated derivative of the psf with x
  Kernel Dx;

  //! tabulated derivative of the psf with y
  Kernel Dy;

  void Resize(const int Hx, const int Hy);

  void Tabulate(const Point& Pt, const DaoPsf& Dao, const int Radius);

  void Tabulate(const Point& Pt, const DaoPsf& Dao, const Window& Rect);

  //! rescale psf and derivative by a factor
  void Scale(const double& scale);

  //! return the current norm of the psf
  double Norm() const { return sum(); }

  //! enable "cout << TabulatedPsf << endl"
  //friend ostream& operator << (ostream & stream, const TabulatedPsf& Psf);
};


class SimFitRefVignet : public Vignet {

public:
  
  //!
  CountedRef<DaoPsf> psf;

  //! build a rough initial galaxy to start the iterative fit
  void makeInitialGalaxy();

  //! empty constructor allocate nothing
  SimFitRefVignet() { DaoPsf *p = 0; psf = p; }

  //! extract vignet image and weight from a ReducedImage of a max radius of Nfwhm*Fwhm pixels.
  //! also initialize Psf if any, but assumes Kern not present.
  SimFitRefVignet(const ReducedImage *Rim, const int Radius);

  //! extract vignet image and weight from a ReducedImage of a max radius of Nfwhm*Fwhm pixels.
  //! also initialize Psf if any, but assumes Kern not present.
  SimFitRefVignet(const PhotStar *Star, const ReducedImage *Rim, const int Radius);

  // default destructor, copy constructor and assigning operator are OK

  //! tabulated PSF and its derivatives to allow fast computation
  TabulatedPsf Psf;

  //! the galaxy underneath at its resolution
  Kernel Galaxy;

  //!
  void UpdatePsfResid();

  //!
  void Load(const PhotStar *RefStar);

  void Resize(const int Hx, const int Hy);

  //! enable "cout << SimFitRefVignet << endl"
  //friend ostream& operator << (ostream& stream, const SimFitRefVignet& myVignet);

};

//! A Vignet with tabulated PSF, and Kernel
class SimFitVignet : public Vignet {

public:

  //! empty constructor allocate nothing
  SimFitVignet() : FitFlux(false) {}

  SimFitVignet(const ReducedImage *Rim);

  //! extract vignet image and weight from a ReducedImage of a max radius of Nfwhm*Fwhm pixels
  //! loads the Kern with the Reference, and also builds the proper Psf.
  SimFitVignet(const PhotStar *Star, const ReducedImage *Rim, const SimFitRefVignet& Ref);
  

  // default destructor, copy constructor and assigning operator are OK

  bool FitFlux;

  bool DontConvolve;

  //! a kernel to convolve a reference to match the current vignet
  Kernel Kern;

  //! tabulated PSF and its derivatives to allow fast computation
  TabulatedPsf Psf;

  void Resize(const int Hx, const int Hy);

  //! 
  void BuildKernel(const ReducedImage *Ref);

  //! update psf and residuals assuming the model = Star->flux[Kern*RefPsf] + Kern*RefGal + Star->sky
  void UpdatePsfResid(const SimFitRefVignet& Ref);

  //! update psf and residuals assuming the model = Star->flux[Kern*RefPsf] + Star->sky
  void UpdatePsfResid(const TabulatedPsf& RefPsf);

  //! update residuals assuming the model = Star->flux*Psf + Kern*RefGal + Star->sky
  void UpdateResid(const Kernel& RefGal);

  //! update residuals assuming the model = Star->flux*Psf + Star->sky
  void UpdateResid();

  //! enable "cout << SimFitVignet << endl"
  friend ostream& operator << (ostream & stream, const SimFitVignet& myVignet);
};


    
#endif // SIMFITVIGNET__H
