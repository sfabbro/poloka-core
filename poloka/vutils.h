#ifndef VUTILS__H
#define VUTILS__H

#include <vector>

/*! \file 
   \brief utilities around vector and matrix algebra
   */

//! return the median of the array. The array is mixed up. 
/*!  Use next function if you want to keep your array untouched  */
double DArrayMedian(double *array, const int size);

//! same as above but does not mix up the array. (calls previous one on a copy of its input ... ). 
double DConstArrayMedian(const double *array, const int size);

//!
float FArrayMedian(float *array, const int size);

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

//! compute median and M.A.D. = median(|x - median(x)|)
//! robust estimator of standard deviation
double median_mad(std::vector<double>& x, double& disp);

//! compute median of a vector
double median_mad(std::vector<double>& x);

//template<class T>  T ScalProd(const T A[], const T B[], int N);
#endif /* VUTILS__H */
