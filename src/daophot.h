// This may look like C code, but it is really -*- C++ -*-
//  daophot.h
//
// Last change: Time-stamp: <07 Mar 2003 at 19:52:07 by Sebastien Fabbro>
//
//
#ifndef DAOPHOT__H
#define DAOPHOT__H
#include <string>
#include "daophotpsf.h"

using namespace std;
class ReducedImage;
const int NMINVARPSF = 10;   // the minimum number of stars to build a variable PSF.
const int NOPT = 30;  

class AperOpt {
public:
  AperOpt(const double &Fwhm);
  AperOpt();
  ~AperOpt();
  float *Radius;
  int Naper;
  float InnerSky;
  float OutterSky;
  void Write(const string ApFile="photo.opt") const;
};

//! dump the set options on screen
void DumpDaophotOptions(const float *Opt);

void WriteDaophotOptions(const float *Opt, const string FileName="daophot.opt", 
			 const int First=0, const int Last=NOPT);

//! returns the order of the spatial variations of the PSF depending on the number of stars
int DaophotPsfVariability(const int nstars);

//! a wrapper to most DAOPHOT routines
class Daophot {

private:

  int open;
  const char **opt_name;
  float *opt_min, *opt_max, *opt_default, *opt, *data;  
  float global_sky, lowbad, threshold;
  string rootname;
  DaoPsf *psf;

  void init_opt();
  void check_opt();

public:
  Daophot();
  Daophot(const string &FitsName);
  Daophot(const ReducedImage &Rim);
  ~Daophot();

  void WriteOptions(const string FileName="daophot.opt", 
		    const int First=0, const int Last=NOPT-1) const;

  void Option(const string OptFile="daophot.opt") const;
  void SetOptions(const ReducedImage &Rim);
  void AddStars();
  void AllStar(const string& FileWithStarsToFit) const;
  void Attach(const string FitsName);

  void Find() const;
  void Group() const;
  void Nstar() const;
  void Peak() const;
  void PeakFit(SEStar &Star) const;
  void Pick() const;
  void Photometry() const;
  void Psf(const int Variability = 0, const bool Manual=false);
  void Sky(float &SkyMean, float &SkyMedian, float &SkyMode, float &SkySigma) const;
  void SubStar() const;

  void IterPsf(const ReducedImage &Rim, const bool Manual=false);
  void PrecisePsf(const ReducedImage &Rim);
};

#endif // DAOPHOT__H
