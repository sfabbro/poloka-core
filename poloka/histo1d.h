#ifndef HISTO1D__H
#define HISTO1D__H

#include <iostream>
#include <cmath>


class Histo1d {

 private:
  float *data;
  int nx;
  float minx;
  float scalex;

 public:
  Histo1d() { data=NULL;}
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
  void dump(std::ostream &stream) const;
  friend std::ostream& operator << (std::ostream &stream, const Histo1d &h)
    { h.dump(stream); return stream;}
  
  ~Histo1d() { if (data) delete [] data;}
};

#endif /* HISTO1D__H */
