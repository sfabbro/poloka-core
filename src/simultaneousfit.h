// This may look like C code, but it is really -*- C++ -*-
#ifndef SIMULTANEOUSFIT__H
#define SIMULTANEOUSFIT__H

const int FitFlux = 1;
const int FitPos  = 2;
const int FitGal  = 4;
const int FitSky  = 8;

#include "imagelist.h"

#include <lafnames.h> 
#include LA_GEN_MAT_DOUBLE_H
#include LA_VECTOR_DOUBLE_H 

/*! 
  \file simultaneousfit.h
  \brief A class to fit pixels simultaneously from a set of vignets.

  The model on a given vignet to be fitted is 
  Model(i,j) = flux*PSF(i-xc,j-yc) + [Gal X Kernel](i,j) + sky
  where user chooses what to fit (flux, position, galaxy pixels and constant sky). 
  The flux and sky are free to change on each vignets, but the galaxy pixels 
  and the point source position are fixed at the best seeing resolution and photometric ratio vignet.
  PSF and Kernel should be given and will not be fitted.
*/

class VignetFit;
class Vignet;

typedef enum {LevenbergMarquardt,
	      NewtonRaphson} MinimMethod;

typedef ImageList<VignetFit>::iterator VignetFitIterator;
typedef ImageList<VignetFit>::const_iterator VignetFitCIterator;

//! Simultaneous fitting of Vignets
class SimultaneousFit : public  ImageList<VignetFit> {

  bool refill;
  // flag to check whether we fit the point source fluxes
  bool fit_flux;
  // flag to check whether we fit the galaxy
  bool fit_gal;
  // flag to check whether we fit the sky 
  bool fit_sky;
  // flag to check whether we fit the point source position
  bool fit_pos;
  // b,
  LaVectorDouble Vec;
  // A,
  LaGenMatDouble Mat;
  LaGenMatDouble MatGal;
  // matrix sizes: number of fluxes to fit, pixels, and total matrix size
  int nfx,nfy,hfx,hfy;
  // matrix indices where thigs are stored
  int fluxstart,galstart,skystart;
  int xind,yind;
  int fluxend,galend,skyend;
  // scaling factor 
  double scale,minscale;
  int nparams,ndata,iter;
  double chi2;
  // returns the galaxy matrix index given pixel (i,j)
  inline int galind(const int i, const int j) const;
  // matrix A and vector b filling routines
  void fillFluxFlux();
  void fillFluxPos();
  void fillFluxGal();
  void fillFluxSky();
  void fillPosPos();
  void fillPosGal();
  void fillPosSky();
  void fillGalGal();
  void fillGalSky();
  void fillSkySky();

  bool printlevel;

public:
  SimultaneousFit();
  ~SimultaneousFit();
  MinimMethod Minim;
  //! fill the entire matrix and vectors
  void FillMatAndVec();
  //! solve the linear equation system
  bool Solve(const double &lambda=0,const bool invert=false);
  //! iterate on solution and solve the system
  bool IterateAndSolve(const int MaxIterations, const double Epsilon=0.01);
  //! assign the VignetFit star 
  void AssignStar();
  //! compute total chi2 and residuals of the fit
  double ComputeChi2();
  //! apply corrections of the last solution with a scale factor
  bool ApplyCorrections(const double Factor=1);
  //! update all vignets from last model fitted and resid
  void MakeResidsAndWeights();
  //! redo all vignets integrated and convolved psfs and derivatives at the last fitted position
  void MakePsfs();
  //! fill, fit and shit
  void DoTheFit(const double &MaxScale=1);
  //! set what you want to fit
  void SetWhatToFit(const int ToFit = FitFlux | FitGal);
  //! get the minimum scaling factor to resize the vignets
  void FindMinimumScale(const double &WorstSeeing);
  //! build a rough initial galaxy to start the iterative fit
  void MakeInitialModel();
  //! resize all the vignets of a scale factor
  void Resize(const double &ScaleFactor);
  //! pointer to the best seeing vignet
  VignetFit *VignetRef;
  //! pointer to the reconstructed galaxy centered double image
  Vignet *galaxy;
  //! compute galaxy flux consistently with the rest
  double GetGalaxyFlux(double &VarGalFlux) const;
  //! compute zero flux as weighted average of zero nights
  double GetZeroFlux(double &VarZeroFlux) const;
  //! write
  void write(const string &MiddleName) const;
  //! returns total chi2 of the fit
  double Chi2() const { return chi2;}
  void GetChi2(double &newchi2, const double &oldchi2, double &lambda);
  //! returns number of degrees of freedom
  int Dof() const { return nparams-ndata;}
  //! returns the current scale factor for vignets
  double Scale() const {return scale;}
};

#endif // SIMULTANEOUSFIT__H



