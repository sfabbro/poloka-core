#include <iomanip>
#include "lapackutils.h"
#include "dimage.h"
#include <lafnames.h> 
#include <lapack.h>
#include LA_GEN_MAT_DOUBLE_H
#include LA_GEN_FACT_DOUBLE_H
#include LA_VECTOR_DOUBLE_H   
#include LA_SPD_MAT_DOUBLE_H
#include <blas++.h>
#include LA_SOLVE_DOUBLE_H
#include LA_UTIL_H

// this following linear solver don't come out of lapack++. 
// here we need another solver to have access to the covariance matrix
extern "C" 
{
  // inversion for general matrix using LU factorization
  void F77NAME(dgetri)(integer *N, doublereal *A, integer *lda, 
		       integer * ipiv, doublereal *b, integer *ldb, integer *info);
  
  // inversion for symm.pos.def matrix using Cholesky factorization
  void F77NAME(dpotri)(char *uplo, integer *N, doublereal *A, integer *lda, integer *info);
  
  // sophisticated routine for Cholesky factorization
  void F77NAME(dposvx)(char *fact, char *uplo, integer *N, integer *nhrs, doublereal *A, integer *lda,
		       doublereal *AF, integer *ldaf, char *equed, doublereal *s, doublereal *b, 
		       integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
		       doublereal *berr,  doublereal *work, integer *iwork, integer *info);

  // least square SVD routine
  void F77NAME(dgelsd) (integer *M, integer *N, integer *NRHS, doublereal *A, integer *LDA, 
			doublereal *B, integer *LDB, doublereal *S, doublereal *RCOND, integer *RANK,
			doublereal *WORK, integer *LWORK, integer *IWORK, integer *INFO);

}


int DivideAndConquerSolve( LaGenMatDouble& A, LaVectorDouble& x, const LaVectorDouble& b)
{  
  x.copy(b);
  integer N = A.size(0);
  long int ldx = x.inc() * x.gdim(0);
  integer nhrs=1;
  LaGenMatDouble toDo = A;
  long int lda = A.inc(0) * A.gdim(0);
  doublereal *s = new doublereal[ldx];  
  doublereal rcond = -1;
  integer lwork = -1;//i don't want to figure this one out and seems to 
  // be rather important . that's probaly why this routine does not work
  doublereal *work  = new doublereal(N);
  integer *iwork = new integer[3*N];
  integer info;
  F77NAME(dgelsd) (&N, &N, &nhrs, &toDo(0,0), &lda, 
		   &x(0), &ldx, s, &rcond, &N,
		   work, &lwork, iwork, &info);
  cerr << "cond=" << rcond << endl;
  cerr << "info=" << info << endl;
  cerr << x(LaIndex(0,9)) << endl;
  delete [] work;
  delete [] s;
  delete [] iwork;
  return info;
}

// compute numerical residuals
double residual(LaGenMatDouble &A, const LaVectorDouble &x, const LaVectorDouble& b)
{
  int M = A.size(0);
  int N = A.size(1);
  (A*x).info();
  b.info();
  cout << setiosflags(ios::floatfield);
  cout << " Numerical residuals: " << endl;
  cout << "  Norm(A*x-b) = " << Norm_Inf(A*x-b) << endl;
  cout << "  Norm(A) = " << Norm_Inf(A) << endl;
  cout << "  Norm(x) = " << Norm_Inf(x) << endl;
  cout << "  Macheps = " << Mach_eps_double() << endl;
  
  if (M>N)
    {
      LaVectorDouble Axb = A*x-b;
      LaVectorDouble R(M,1);

      Blas_Mat_Trans_Vec_Mult(A, Axb, R);
      return Norm_Inf(R) / 
	(Norm_Inf(A)* Norm_Inf(x) * N * Mach_eps_double());
    }
  else
    {
      return Norm_Inf(A*x-b ) /
	( Norm_Inf(A)* Norm_Inf(x) * N * Mach_eps_double());
    }
}

void operator += (LaGenMatDouble &A, const LaGenMatDouble &B)
{
  int M = A.size(0);  int N = A.size(1);
  if (M != B.size(0) || N != B.size(1)) 
    {cerr << " trying to add non-conformant matrices" << endl; return;}
  for (int i=0;  i<M; i++) for(int j=0; j<N; j++) A(i,j) += B(i,j);
}


// few routines to check matrices and vector are properly filled
void writeMatFits(const LaGenMatDouble &mat, const string &name)
{
  DImage im(mat.size(0), mat.size(1));
  for (int i=0; i<mat.size(0); ++i) 
    for (int j=0; j<mat.size(1); ++j) 
      im(i,j) = mat(i,j);
  im.writeFits(name);
}

void readMatFits(LaGenMatDouble &mat, const string &name)
{
  DImage im(name);
  mat.resize(im.Nx(), im.Ny());
  for (int i=0; i<mat.size(0); ++i) 
    for (int j=0; j<mat.size(1); ++j) 
      mat(i,j) = im(i,j);
}

void writeVecImage(const LaVectorDouble &vec, int nx, int ny, const string &name)
{
  DImage im(nx,ny);
  for (int m=0; m<vec.size(); ++m) 
    {
      int i = m / ny;
      int j = m % ny;
      im(i,j) = vec(m);
    }
  im.writeFits(name);
}


// so far uses LU factorization to invert. SVD is on the way to be implemented. See below.
int LUSolveAndInvert( LaGenMatDouble& A, LaVectorDouble& x, const LaVectorDouble& b)
{  

  if ( A.inc(0) != 1 || A.inc(1) != 1)
    cerr << " LU Linear solver error : A non-contiguous" << endl;
  if (x.size() != b.size())
      cerr << " LU Linear solver error : x and b not same size " << endl;
  x.copy(b);            // will throw exception if not conformant
  if (A.size(0) != A.size(1))
    cerr << " LU Linear solver: Square matrix expected." << endl; 
  long int info;
  int M = A.size(0);
  long Ml = M;
  long int K = 1;
  long int lda = A.inc(0) * A.gdim(0);
  long int ldx = x.inc() * x.gdim(0);
  LaVectorLongInt ipiv(M);        
  //  LaGenMatDouble old_A=A;
  F77NAME(dgesv) (&Ml, &K, &A(0,0), &lda, &ipiv(0), &x(0), &ldx, &info);
  if (info!=0) 
    {
      cerr << " WARNING: Factorization failure . info =" << info << "(>0 : U is singular)" << endl;
      return info;
    }
  long int lwork = ldx;
  LaVectorDouble work(b);
  F77NAME(dgetri) (&Ml, &A(0,0), &lda, &ipiv(0), &work(0), &lwork, &info);
  if (info!=0) cerr << " WARNING: Inversion failure . info = " << info << endl;
  //  LaGenMatDouble Ident = A*old_A;
  //  for (int i=0;i<Ml;++i)cout << Ident(i,i) << endl;
  return info;
}

// LaSpdMatDouble is not really convenient to use. Use LaGenMatDouble, but Spd linear solver
// Note also Lapack++ does not have anything for inverting. So recode the whole schmol

int SpdSolveAndInvert(LaGenMatDouble& A, LaVectorDouble &x, const LaVectorDouble &b)
{  
  char uplo = 'L';
  integer N = A.size(0), nhrs = 1, info = 0;
  integer lda = A.inc(0) * A.gdim(0);
  integer ldb = b.inc() * b.gdim(0);
  LaSpdMatDouble toInvert(N,N);
  for (int i=0; i<N; ++i) for (int j=i; j<N; ++j) toInvert(i,j) = A(j,i);

#ifdef EXPERT //not totally understood yet.
  char fact = 'E';
  char equed;
  // need a copy bcause lapack will change it
  LaVectorDouble bcopy(ldb);
  bcopy.copy(b);
  doublereal *s = new doublereal[N];
  doublereal *work = new doublereal[3*N];
  integer *iwork = new integer[N];
  doublereal rcond, ferr, berr;
  LaSpdMatDouble AF(N,N);
  cout << " Factorize " << endl;
  //this routine factorize, equilibriate, and solve but something is yet to be implemented
  F77NAME(dposvx)(&fact, &uplo, &N, &nhrs, &toInvert(0,0), &lda, &AF(0,0), &lda, 
		  &equed, s, &bcopy(0), &ldb, &x(0), &ldb, &rcond, &ferr, &berr, 
		  work, iwork, &info);
  cout << " Equilibrated ? " << equed << endl;
  cout << " Forward error " << ferr << endl;
  cout << " Backward error " << berr << endl;
  cout << " Condition number = " << rcond << endl;
  delete [] s;
  delete [] work;
  delete [] iwork;
#endif

  F77NAME(dposv)(&uplo, &N, &nhrs, &toInvert(0,0), &lda, &x(0), &ldb, &info);
  if (info!=0)
    {
      cerr << " WARNING: Factorization failure . info =" << info <<  " (>0 is not pos.def)" << endl;
      return info;
    }

  cout << " Now invert using previous factorization" << endl;
  F77NAME(dpotri)(&uplo, &N, &toInvert(0,0), &lda, &info);
#ifdef DEBUG
  cout << " Checking " << endl;
  LaGenMatDouble Ident(N,N);
  Blas_Mat_Mat_Mult(toInvert,A,Ident);
  double total=0;
  for (int i=0;i<N;++i) total +=Ident(i,i);
  cout << " check should be 1: " << total/N << endl;
#endif
  if (info!=0)
    {
      cerr << " WARNING: Inversion failure . info = " << info << endl;
      return info;
    }
  A = LaGenMatDouble(toInvert);
  return info;
}

int SpdSolveAndInvert(LaGenMatDouble& A, LaVectorDouble &b, const bool invert)
{  
  char uplo = 'L';
  integer N = A.size(0), nhrs = 1, info = 0;
  integer lda = A.inc(0) * A.gdim(0);
  integer ldb = b.inc()  * b.gdim(0);
  F77NAME(dposv)(&uplo, &N, &nhrs, &A(0,0), &lda, &b(0), &ldb, &info);
  if (info != 0) 
    {
      cerr << " WARNING: Factorization failure . info =" << info <<  " (>0 is not pos.def)" << endl;
      return info;
    }
  if (invert)
    {
      F77NAME(dpotri)(&uplo, &N, &A(0,0), &lda, &info);
      if (info!=0) cerr << " WARNING: Inversion failure . info = " << info << endl;
    }
  return info;
}


int SpdSolve(LaGenMatDouble& A, LaVectorDouble &b)
{  
  char uplo = 'L';
  integer N = A.size(0), nhrs = 1, info = 0;
  integer lda = A.inc(0) * A.gdim(0);
  integer ldb = b.inc() * b.gdim(0);
  F77NAME(dposv)(&uplo, &N, &nhrs, &A(0,0), &lda, &b(0), &ldb, &info);
  if (info != 0) 
    cerr << " WARNING: Factorization failure . info =" << info <<  " (>0 is not pos.def)" << endl;
  return info;
}

int SpdInvert(LaGenMatDouble& A)
{  
  char uplo = 'L';
  integer  N = A.size(0), info = 0;
  integer lda = A.inc(0) * A.gdim(0);
  F77NAME(dpotri)(&uplo, &N, &A(0,0), &lda, &info);
  if (info!=0) cerr << " WARNING: Inversion failure . info = " << info << endl;
  return info;
}


#ifdef SVD // singular value decomposition to be implemented. 
extern "C"
{
  // SVD via divide and conquer factorization fortran lapack routine
  void F77NAME(dgesdd)( char *jobz, integer *m, integer *n,
			doublereal *a, integer *lda, double *s, 
			doublereal *u, integer *ldu, 
			doublereal *vt, integer *ldvt, 
			doublereal *work, integer *lwork, 
			integer *info);
}

int SVDLinearSolver(LaGenMatDouble& A, LaVectorDouble& x, const LaVectorDouble& b)
{  

  if ( A.inc(0) != 1 || A.inc(1) != 1)
    cerr << " SVD Linear solver error : A non-contiguous" << endl;
  if (x.size() != b.size())
      cerr << " SVD Linear solver error : x and b not same size " << endl;
  x.copy(b); 
  long int info;
  char jobvz = 'A';
  int N = A.size(0);
  LaVectorDouble sigma[N];
  LaGenMatDouble U(N,N);
  LaGenMatDouble Vt(N,N);
  long int lda = A.inc(0) * A.gdim(0);
  long int ldx = x.inc() * x.gdim(0);
  F77NAME(dgesdd) (&jobz, &N, &N, &A(0,0), &lda, &sigma(0), &U, &ldu, 
		   &Vt, &ldvt, &work, &lwork, &info);  
  long int lwork = ldx;
  LaVectorDouble work(b);
  
  return info;
}




#endif
