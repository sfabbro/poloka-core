// This may look like C code, but it is really -*- C++ -*-
#ifndef LAPACKUTILS__H
#define LAPACKUTILS__H

class LaGenMatDouble;
class LaVectorDouble;
#include <string>
using namespace std;


// so far uses LU factorization to invert. SVD is on the way to be implemented. See below.
int DivideAndConquerSolve( LaGenMatDouble& A, LaVectorDouble& x, const LaVectorDouble& b);
int LUSolveAndInvert(LaGenMatDouble &A, LaVectorDouble &x, const LaVectorDouble &b);
int SpdSolveAndInvert(LaGenMatDouble &A, LaVectorDouble &x, const LaVectorDouble &b);
//! solve and replace without copying
int SpdSolveAndInvert(LaGenMatDouble& A, LaVectorDouble &b, const bool invert=true);

int SpdSolve(LaGenMatDouble& A, LaVectorDouble &b);
int SpdInvert(LaGenMatDouble& A);

//! A += B for matrices
inline void operator += (LaGenMatDouble &A, const LaGenMatDouble &B);

//! compute numerical residuals
double residual(LaGenMatDouble &A, const LaVectorDouble &x, const LaVectorDouble& b);

//! few routines to write matrices and vectors as FITS files
void writeMatFits(const LaGenMatDouble &mat, const string &name);
void writeVecImage(const LaVectorDouble &vec, int nx, int ny, const string &name);
void readMatFits(LaGenMatDouble &mat, const string &name);


#endif //LAPACKUTILS__H
