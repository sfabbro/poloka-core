#ifndef MATVECT__H
#define MATVECT__H

#include <iostream>
#include <string>

class Vect;
class Mat;

#define MATVECT_CHECK_BOUNDS

// Routines for solving linear systems + inversion 
//==================================================================



// solving linear system A.X = B 
// Uses lapack dposv_
// Matrix A is assumed to be symmetric (you'll get a core dump if it is not
// (actually this is not compulsory with dposv (see 'man dposv') but we do it
// for consistency).
// You just need to fill half (n*(n+1)/2 parameters) of the matrix
// if you have filled matrix M parameters for which y>=x (with M(x,y)=...), use UorL="L"
// else use UorL="U"
// Matrix A is modified in this routine (also B which is at the end the solution X)
int cholesky_solve(Mat &A, Vect &B, char* UorL = "L");

// Inverts matrix A using the factorization done in cholesky_solve
// Uses lapack dpotri_
// Matrix A is assumed to be symmetric (you'll get a core dump if it is not
// (actually this is not compulsory (see 'man dptri') but we do it
// for consistency).
// This routine must be called after cholesky_solve, make sure the value of UorL
// is the same as that used with dposv
int cholesky_invert(Mat &A, char* UorL = "L"); // when cholesky_solve is called first


// Mat and Vect classes
//====================================================================


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

  
  friend std::ostream& operator << (std::ostream &stream, const Mat &m)
    { m.writeASCII(stream); return stream;}
  
  // get a block of this matrix as a new matrix
  // size (x_max-x_min+1)*(y_max-y_min+1) (both min and max included)
  // remember  0 <= x < nx ,  0 <= y < ny
  Mat SubBlock
    (unsigned int x_min,unsigned int x_max,unsigned int y_min,unsigned int y_max) const;
  
  Mat WithoutRows(unsigned int y_min,unsigned int y_max) const;
  Mat WithoutColumns(unsigned int x_min,unsigned int x_max) const;


  // i/o in fits for matrices
  int readFits(const std::string &FitsName);
  int writeFits(const std::string &FitsName) const;
  int readASCII(std::istream& Stream);
  int readASCII(const std::string &FileName);
  int writeASCII(std::ostream& Stream) const;
  int writeASCII(const std::string &FileName) const;
  
  void Symmetrize(const char* UorL = "L");
  
  // inverts a symetric posdef matrix using DSINV  CERNLIB's routine
  int SymMatInvert();

  // operators
  Mat operator +(const Mat& Right) const;
  Mat operator -(const Mat& Right) const;
  Mat operator *(const Mat& Right) const;
  Mat operator *(const Vect& Right) const;
  Mat & operator =(const Mat& Right);
  
  void operator +=(const Mat& Right);
  void operator -=(const Mat& Right);
  void operator *=(const Mat& Right);
  
  Mat operator *(const double Right) const;
  friend Mat operator *(const double Left, const Mat &Right);
  void operator *=(const double Right);

  operator double() const;
  operator Vect() const;
  Mat transposed() const;

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

   int writeASCII(std::ostream& Stream) const;
   int readASCII(std::istream& Stream);


  void dump(std::ostream& Stream) const;
  friend std::ostream& operator << (std::ostream &stream, const Vect &v)
    { v.dump(stream); return stream;}

  // operators
  Vect operator +(const Vect& Right) const;
  Vect operator -(const Vect& Right) const;
  Vect & operator =(const Vect& Right);

  double operator *(const Vect& Right) const; // scalar product
  
  void operator +=(const Vect& Right);
  void operator -=(const Vect& Right);
  
  

  Vect operator *(const double Right) const;
  friend Vect operator *(const double Left, const Vect &Right);
  void operator *=(const double Right);

  Mat transposed() const;
  operator Mat() const;
  operator double() const;
};

#endif /*MATVECT__H */
