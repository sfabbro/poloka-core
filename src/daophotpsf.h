// This may look like C code, but it is really -*- C++ -*-
//  daophotpsf.h
//
// Last change: Time-stamp: <05 Mar 2003 at 15:00:15 by Sebastien Fabbro>
//
//
#ifndef DAOPHOTPSF__H
#define DAOPHOTPSF__H

#include "reducedimage.h"
#include "sestar.h"
#include "dimage.h"

class DaoPsf {

  
public:
  DaoPsf(const string &FileName);
  DaoPsf(const DbImage &DbIm);
  DaoPsf() : param(0),table(0) {}

  ~DaoPsf();

  float *param, *table;
  float psfmag, bright, xpsf, ypsf;
  int npsf, npar, nexp, nfrac, type;
  double radius;

  string Type() const;
  void Allocate();
  bool read(const string &FileName);
  void dump(ostream &Stream=cout) const;
  double Xpsf() const {return xpsf;}
  double Ypsf() const {return ypsf;}
  double HwhmX() const {return param[0];} 
  double HwhmY() const {return param[1];}
  double ThetaXY() const;
  double Radius() const {return radius;}
  double PsfMag() const {return psfmag;}
  double Value(int i, int j, const double &Xc, const double &Yc, 
	       double &DpDx, double &DpDy) const;
  void MakeStar(Kernel &ImStar, const BaseStar &Star) const;
};

//! 
class PsfConditions {
private:
  SEStarList catalog;
  double fw_min, fw_max, mx_min, mx_max, my_min, my_max, 
    flux_min, flux_max, fm_max, cstar_min, neigh_rad, edges_rad;
  int flag_min, flag_max, flagbad_max;
  Frame usableFrame;
public:
  PsfConditions();
  PsfConditions(const ReducedImage &Rim);
  bool AreVerified(const SEStar* star) const;
};

//! a class to select PSF stars from different ways
class PsfStars : public SEStarList {
private:
  const ReducedImage *rim;
public:
  //! just reads in the decent objects from a ReducedImage
  PsfStars(const ReducedImage &Rim);
  ~PsfStars(){}
  //! removes the stars with bad flags on the neighbor file created during PSF building
  size_t FilterNei();
  //! removes the stars from the file created by ALLSTAR, given some cuts
  size_t FilterAls(const double ChiMax=3.0, const double SharpMin=-0.5, const double SharpMax=0.5);
  //! removes the stars which do not verify the conditions Cond
  size_t FilterCond(const PsfConditions &Cond);
  //! keeps the stars in the set of images in RimList
  size_t FilterMultiIm(const ReducedImageList &RimList, const ReducedImage *ref, const double MaxDist=2);
};

//! Update a ReducedImage with its associated DAOPHOT PSF file
bool UpdateSeeingFromDaoPsf(ReducedImage &Rim);

#endif // DAOPHOTPSF__H
