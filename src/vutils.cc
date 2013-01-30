#include <math.h>
#include "vutils.h"
#include <algorithm>
#include <string.h> // memcpy 

using namespace std;

double DArrayMedian(double *array, const int size)
{
  sort(array, array+size);
  return size&1? array[size/2] : (array[size/2-1] + array[size/2])*0.5;
}

double DConstArrayMedian(const double *array, const int size)
{
  double *to_delete = new double[size*sizeof(double)];
  memcpy(to_delete, array, size*sizeof(double));
  double median = DArrayMedian(to_delete, size);
  delete [] to_delete;
  return median;
}
 

float FArrayMedian(float *array, const int size)
{
  sort(array, array+size);
  return size&1? array[size/2] : (array[size/2-1] + array[size/2])*0.5;
}


void Dmean_median_sigma(double *values, const int nval, double &mean,  double &median, double &sigma)
{
  mean =0;
  sigma =0;
  median = DArrayMedian(values, nval);
  for (int i=nval-1; i >= 0 ; --i)
    {
      mean += values[i];
      sigma += values[i]*values[i];
    }
  mean /= double(nval);
  sigma = sigma/double(nval) - mean*mean;
  if (sigma>0)  sigma = sqrt(sigma); else sigma = 0;
}

void DConst_mean_median_sigma(const double *array, const int size, double &mean,  
			      double &median, double &sigma) {
  double *tarray=new double[size];
  memcpy(tarray,array,size*sizeof(double));
  Dmean_median_sigma(tarray,size,mean,median,sigma);
  delete [] tarray;
}

void Fmean_median_sigma(float *values, const int nval, float &mean, float &median, float &sigma)
{
  double dmean =0;
  double dsigma =0;
  median = FArrayMedian(values, nval);
  for (int i=nval-1; i>=0 ; --i)
    {
      dmean += values[i];
      dsigma += values[i]*values[i];
    }
  dmean /= double(nval);
  dsigma = dsigma/double(nval) - dmean * dmean;
  mean = dmean;
  if (dsigma>0)  sigma = sqrt(dsigma); else sigma = 0;
}

float Fmedian_sigma(float *values, const int nval, float &sigma)
{
  double dmean =0;
  double dsigma =0;
  double median = FArrayMedian(values, nval);
  for (int i=nval-1; i>=0 ; --i)
    {
      dmean += values[i];
      dsigma += values[i]*values[i];
    }
  dmean /= double(nval);
  dsigma = dsigma/double(nval) - dmean * dmean;
  if (dsigma>0)  sigma = sqrt(dsigma); else sigma = 0;
  return median;
}



double clipmean(double *values, int &nval, double &sigma, const double &k, const int niter)
{
  double mean,median,sigmat;
  Dmean_median_sigma(values,nval,mean,median,sigmat);
  double clip = k * sigmat;
  double low = median - clip;
  double high = median + clip;
  int j,n=0;
  int nold = nval;
  if (nval==1)
    {
      sigma = 0.0;
      return values[0];
    }
  for (int i=0; i<niter;i++)
    {
      mean = 0;
      sigma = 0;
      n = 0;
      for (j=0;j<nval;++j)
	{
	  if (values[j] > low && values[j] < high ) 
	    {
	      n++;
	      mean += values[j];
	      sigma += values[j]*values[j];
	    }
	}
      mean /= double(n);
      sigma = sigma/double(n) - mean*mean;
      if (sigma >0 ) sigma = sqrt(sigma);
      else {sigma = 0; break;}
       clip = k * sigma;
      if (nold == n) break;
      low = mean - clip;
      high = mean + clip;
      nold = n;
    }
  nval = n;
  return mean;
}




