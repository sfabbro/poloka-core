#ifndef VUTILS__H
#define VUTILS__H

/*! \file 
   \brief utilities around vector and matrix algebra
   */

/*! Diagonalize Real Symmetric Matrix matrix(n*n) 
 * memory for eigenvectors of size n*n must be allocated by the user
 * memory for eigenvalues of size  n must be allocated by the user
 * \return 0 if ok
 */
int DiagonalizeRealSymmetricMatrix(int n,double * matrix, double * eigenvectors , double * eigenvalues);
 

//! return the median of the array. The array is mixed up. 
/*!  Use next function if you want to keep your array untouched  */
double DArrayMedian(double *array, const int size);

//! same as above but does not mix up the array. (calls previous one on a copy of its input ... ). 
double DConstArrayMedian(const double *array, const int size);

//!
float FArrayMedian(float *array, const int size);

//! inverts and solve using various CERNLIB routines
int MatSolve(double *A, const int N, double *B);

//! inverts a symetric posdef matrix using DSINV  CERNLIB's routine
int SymMatInv(double *A, const int N);

//!
double ScalProd(const double A[], const double B[], const int N);

//!
void MatVec(const double *M, const int N1, const int N2, const double *V, double *R);

//!
double VecMatVec(const double*V1, const double *M, const int N1, const int N2, const double *V2);

//! returns mean, median and rms of an array.
void Dmean_median_sigma(double *values, const int nval, double &mean,  double &median, double &sigma);

//! returns mean, median and rms of an array
void Fmean_median_sigma(float *values, const int nval, float &mean, float &median, float &sigma);
//! returns median and rms of an array
float Fmedian_sigma(float *values, const int nval, float &sigma);
void DConst_mean_median_sigma(const double *array, const int size, double &mean,  
			      double &median, double &sigma);

//! computes the clipped-mean and sigma with cutting at k-sigma, return mean
double clipmean(double *values, int &nval, double &sigma, const double &k=3.5, const int niter=4);

/*! fit a gaussian on a region about mean of half width k-sigma, return mean
 * if first_evalutation, mean and sigma are guessed with  DConst_mean_median_sigma (using median)
 */
double gaussianfit(const double *values , int nval, double &mean, double &sigma, const double &k=3.5, bool first_evalutation=true);


//template<class T>  T ScalProd(const T A[], const T B[], int N);
#endif /* VUTILS__H */
