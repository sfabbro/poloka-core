// This may look like C code, but it is really -*- C++ -*-
#ifndef DAOPHOTPSF__H
#define DAOPHOTPSF__H

#include <persistence.h>
#include <reducedimage.h>
#include <sestar.h>

//! returns the order of the spatial variations of the PSF depending on the number of stars
int DaoPsfVariability(const int Nstars);


//! a tabulated psf with residuals
class DaoPsf : public RefCount {

  friend class Daophot; 

  float *param, *table;
  float  psfmag, bright, xpsf, ypsf, scale;
  int    npsf, npar, nexp, nfrac, type;
  double radius;
  
  static const int MAXPSF;      // maximum PSF array size allowed by DAOPHOT
  static const int MAXPAR;      // maximum number of PSF parameters
  static const int MAXEXP;      // maximum coeff for interpolating spatial variable table

public:

  //! empty constructor make sure array are not pointing to anything
  DaoPsf() : param(0), table(0) {}

  //! constructor read a daophot psf file
  DaoPsf(const string &FileName);

  //! constructor read a daophot psf file of a DbImage
  DaoPsf(const DbImage &DbIm);
    
  //! free memory
  ~DaoPsf();

  //! allocate the tabulated array of residuals and parameters array
  void Allocate(const int Npsf, const int Npar);

  //! return the analytical type of the PSF
  string Type() const;

  //! useful number to compute residuals
  double Xpsf() const {return xpsf;}

  //! useful number to compute residuals 
  double Ypsf() const {return ypsf;}

  //! return the half-width at half max in x-direction
  double HwhmX() const {return param[0];} 

  //! return the half-width at half max in y-direction
  double HwhmY() const {return param[1];}

  //! return the gyro angle
  double ThetaXY() const;

  //! return the radius of the psf. above it is zero.
  double Radius() const {return radius;}

  //! return the "instrumental" zero point . see daophot manuals
  double PsfMag() const {return psfmag;}

  //! compute the psf value in a pixel (i,j) from a point (Xc,Yc) and its derivatives
  double Value(const int i, const int j, const double &Xc, const double &Yc, 
	       double &DpDx, double &DpDy) const;

  //! read a daophot psf style file
  bool read(const string &FileName);

  //! dump DaoPsf properties
  void dump(ostream &Stream=cout) const;
};

#endif // DAOPHOTPSF__H
