// This may look like C code, but it is really -*- C++ -*-
//

class PhotStar;

class DaoPsfPhot {

  float *data;
  float radius, perr, pkerr, x, y,  deltax, deltay, scale, errmag, sky;
  int lx,ly,nx,ny;
  
public:

  DaoPsfPhot() : data(0) {}

  DaoPsfPhot(const ReducedImage& Rim);

  ~DaoPsfPhot() { if (data) delete [] data; }

  void PeakFit(PhotStar *Star);

  void operator () (PhotStar *Star);
};



