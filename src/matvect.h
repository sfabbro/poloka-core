#ifndef MATVECT__H
#define MATVECT__H

#include <iostream>


#define MATVECT_CHECK_BOUNDS

/* obviously the later class should be borrowed somewhere  */
class Mat {
  
 private:
  double *data;
  unsigned int nx,ny;
  
  
 public:
  
  Mat() : data(NULL), nx(0), ny(0) {};  
  Mat(const unsigned int NX, const unsigned int NY);
  Mat(const Mat& other);
  ~Mat() { delete [] data;}
  
  void allocate(const unsigned int NX, const unsigned int NY);
  
  double operator () (const unsigned int i, const unsigned int j) const;
  double& operator () (const unsigned int i, const unsigned int j);
  
  unsigned int SizeX() const { return nx;}
  unsigned int SizeY() const { return ny;}
  
  void Zero() {memset(data, 0, nx*ny*sizeof(double));};
  void Identity();

  const double* Data() const {return data;};
  double* NonConstData() {return data;};

  void dump(std::ostream& Stream) const;
  friend std::ostream& operator << (std::ostream &stream, const Mat &m)
    { m.dump(stream); return stream;}

  // operators
  Mat operator +(const Mat& Right) const;
  Mat operator -(const Mat& Right) const;
  //Mat operator *(const Mat& Right) const;
  
  void operator +=(const Mat& Right) const;
  void operator -=(const Mat& Right) const;
  //void operator *=(const Mat& Right) const;
  
  Mat operator *(const double Right) const;
  friend Mat operator *(const double Left, const Mat &Right);
  void operator *=(const double Right);
};

class Vect {

 private:
  
  double *data;
  unsigned int n;
  
  
 public:
  
  Vect() : data(NULL), n(0) {};
  Vect(const unsigned int N);
  Vect(const Vect&);
  ~Vect() { delete [] data;}

  void allocate(const unsigned int N);
  
  double operator () (const unsigned int i) const;
  double& operator () (const unsigned int i);
  
  unsigned int Size() const { return n;}
  
  void Zero() {memset(data, 0, n*sizeof(double));};

  const double* Data() const {return data;};
  double* NonConstData() {return data;};

  void dump(std::ostream& Stream) const;
  friend std::ostream& operator << (std::ostream &stream, const Vect &v)
    { v.dump(stream); return stream;}
};

#endif /*MATVECT__H */
