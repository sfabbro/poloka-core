// This may look like C code, but it is really -*- C++ -*-
#ifndef VIGNETFIT__H
#define VIGNETFIT__H

#include "vignet.h"
/*!
  \file vignetfit.h
  \brief a sub-image with its associated tabulated kernels, psfs, ...
  Mainly to use for SimultaneousFit
 */
//! A Vignet with tabulated PSFs, Weights and Kernels 
class PhotStar;
class Night;

class VignetFit : public Vignet {
public:
  VignetFit();
  VignetFit(PhotStar *aStar, const Night *aNight,
	    const Image &Source, const Kernel &aKernel,
	    const int HRefX, const int HRefY);
  VignetFit(PhotStar *aStar, const Night *aNight, const string &Name);
  ~VignetFit(){};

  const Night *night;
  PhotStar *star;

  bool IsStarHere;
  bool IsRefResolution;

  double flux, sky, chi2, residsigma;

  Kernel Kern;
  Kernel Psf;
  Kernel DpDx;
  Kernel DpDy;
  Kernel Weight;
  Kernel Chi2;
  Kernel Resid;
  Kernel Model;
  string Name() const;
  void MakeNormalizedPsf();
  void MakeConvolvedPsf(const Kernel &PsfRef, const Kernel &DpDxRef, const Kernel &DpDyRef);
  void MakeInitialWeight();
  void UpdateWeight(const bool robustify=false);
  void MakeResid();
  void MakeChi2();
  void MakeModel(const Vignet *Galaxy=NULL);
  void QuickPhotometry(const double &Nfwhm);
  double WeightedAperture(double &VarAper);

  void dump(ostream &Stream=cout) const;
#ifndef SWIG
  friend ostream& operator << (ostream & stream, const VignetFit& myVignet) 
    { myVignet.dump(stream); return stream;}
#endif

  void writeAllFits(const string &MiddleName) const;
  void readAllFits(const string &MiddleName);
};

#endif // VIGNETFIT__H
