// This may look like C code, but it is really -*- C++ -*-
#ifndef VIGNET__H
#define VIGNET__H

#include "dimage.h"
#include "countedref.h"

class Image;
/*! \file 
    \brief A resizeable sub-image, with centered coordinates
*/
//! Similar to Kernel, but resizeable and contains some kind of Point
class Vignet : public Kernel, public RefCount {
protected:
  int hx,hy;
public:
  Vignet(const double &X, const double &Y, const Image &Source,  const int HMaxX, const int HMaxY);
  Vignet(const double &X, const double &Y, const int HMaxX, const int HMaxY);
  Vignet(const string &FileName);
  Vignet();
  ~Vignet(){}
  //! coordinates of the center pixel where the Vignet was grabbed from
  int ic,jc;
  //! real coordinates of the star centroid relative to the down-left corner of center pixel (ic,jc)
  double dxc,dyc; 
  //! image coordinates of the down-left corner pixel where the Vignet was grabbed
  int istart,jstart;
  //! image coordinates of the upper-right corner pixel where the Vignet was grabbed
  int iend,jend;
  void Resize(const int HSizeX, const int HSizeY);
  void Resize(const double &ScaleFactor);
  void Detect(const double &PosThresh, const double &NegThresh);
  double Aperture(double &VarAper, const double &Radius, 
		  const double &VarPix, const double &Sky) const;
  double WeightedAperture(double &VarAper, const Kernel &Model, 
			  const double &VarPix, const double &Sky) const;
  void WeightedRecentroid(double &Flux, const double &Sky,
			  const double &SigX, const double &Sigy);
  int Hx() const {return hx;}
  int Hy() const {return hy;}
  void readFits(const string &FileName);
  void writeFits(const string &FileName) const;
};


#include "imagelist.h"
typedef ImageList<Vignet> VignetList;
typedef VignetList::iterator       VignetIterator;
typedef VignetList::const_iterator VignetCIterator;

#endif // VIGNET__H
