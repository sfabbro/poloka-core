#ifndef HISTO1D__H
#define HISTO1D__H

#include <math.h>

#include "persistence.h"

class Histo1d {
  CLASS_VERSION(Histo1d,1);
  #define Histo1d__is__persistent
 private:
  float *data;
  int nx;
  float minx;
  float scalex;

 public:
  Histo1d() {}
  Histo1d(int Nx, float Minx, float Maxx);
  void Fill(float X, float weight=1.)
      {int bin = int(floor((X - minx)*scalex)); 
       if (bin<0 || bin>=nx) return; 
       data[bin] += weight;
      }
  const float *array() const { return data;}
  int Nx() const { return nx;}
  int BinAt(double X) const;
  double BinCenter(int bin) const;
  double Scale() const {return scalex;}
  double Minx() const { return minx;}
  /* returns the contents and set X to the abcissa */
  double  MaxBin(double &X) const;
  double BinWidth() const { return (1./scalex);}
  void ZeroBins(double Xmin,double Xmax);
  void dump(ostream &stream) const;
  friend ostream& operator << (ostream &stream, const Histo1d &h)
    { h.dump(stream); return stream;}
  
  ~Histo1d() { delete [] data;}
};

#endif /* HISTO1D__H */
