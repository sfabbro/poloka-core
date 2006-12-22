#ifndef HISTO2D__H
#define HISTO2D__H


class Histo2d {
 
 private:
  float *data;
  int nx,ny;
  float minx,miny;
  float scalex, scaley;

  bool indices(const double &X, const double &Y, int &ix, int &iy) const;

 public:

  Histo2d() {data=NULL;}
  Histo2d(int nx, float minx, float maxx, int ny, float miny,float maxy);

  void Fill(float x, float y, float weight=1.);

  double MaxBin(double &x, double &y) const ;

  void BinWidth(double &Hdx, double &Hdy) const { Hdx = 1./scalex; Hdy = 1./scaley;}

  double BinContent(const double &X, const double &Y) const;

  void ZeroBin(const double &X, const double &Y);

  ~Histo2d() { if (data) delete [] data;}
};

#endif /* HISTO2D__H */
