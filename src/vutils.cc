#include <math.h>
#include "vutils.h"
#include "histo1d.h"
#include "matvect.h"

#define dsinv dsinv_
#define dfact dfact_
#define dfinv dfinv_
#define dfeqn dfeqn_
#define eisrs1 eisrs1_


extern "C" 
{
  void dsinv(int *N, double *A, int *IDIM, int *IFAIL);
  void dfact(int *N, double *A, int *idim, double *r, int *ifail, double *Det, int *jfail);
  void dfinv(int *n, double *A, int *idim, double *r);
  void dfeqn(int *n, double *a, int *idim, double *r, int *k, double *b);
  void eisrs1(int *NM,int *N,double *AR,double *WR,double *ZR,int *IERR,double *WORK);
}


int DiagonalizeRealSymmetricMatrix(int n,double * matrix, double * eigenvectors , double * eigenvalues) {
  double * work = new double[n];
  int ierr = 0;
  
  eisrs1(&n,&n,matrix,eigenvalues,eigenvectors,&ierr,work);
  
  delete [] work;
  return ierr;
}


int MatSolve(double *A, const int N, double *B)
{
  double *r = new double [N];
  int n = N;
  int ifail, jfail;
  double det;
  dfact(&n, A, &n, r, &ifail, &det, &jfail);
  if (ifail == 0)
    {
      int k=1;
      //  dfinv(&n, A, &n,r);
      dfeqn(&n, A, &n, r, &k, B);
    }
  delete [] r;
  return (!ifail);
}


#include <algorithm>
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

int SymMatInv(double *A, const int N)
{
  int n = N;
  int ierr;
  dsinv(&n,A,&n,&ierr);
  return (!ierr);
}

double ScalProd(const double A[], const double B[], const int N)
{
  double res = 0;
  for (int i = 0; i< N; i++) res += A[i]*B[i];
  return res;
}

void MatVec(const double *M, const int N1, const int N2, const double *V, double *R)
{ /* M is N1 x N2  V is N2, R is N1. 
     M(i,j) = M[i*N2+j]*/
  for (int i=0; i<N1; ++i)
    {
      R[i] = ScalProd( M+(i*N2), V, N2);
    }
}

double VecMatVec(const double*V1, const double *M, const int N1, const int N2, const double *V2)
{
  /* M is N1xN2, V1 is N1, V2 is N2 */
  double sum =0;
  for (int i=0;i<N1; ++i)
    sum += V1[i]*ScalProd(M+(i*N2),V2,N2);
  return sum;
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

static double sqr(double x) {
  return x*x;
}

#define DEBUG
double gaussianfit(const double *values, int nval, double &mean, double &sigma, const double &k, bool first_evalutation)
{
  if(first_evalutation) {
    double median;
    DConst_mean_median_sigma(values,nval,mean,median,sigma);
    mean=median;
#ifdef DEBUG
    cout << "gaussianfit::first_eval  " << nval << " " << mean << " " << sigma << endl; 
#endif
  }

  // first count the number of entries in one sigma about mean histo
  double xmin=mean-2*sigma;
  double xmax=mean+2*sigma;
  int nok=0;
  for(int i=0;i<nval;i++) {
    if (values[i]>xmin && values[i]<xmax) nok++;
  }
  
  // choose the number of bins, min= 3 per sigma
  int nbins = int((2*k)*4*max(nok/200,1)+1);
  
  double center=mean;
  mean=0;
  
  Histo1d histo(nbins,-k*sigma,k*sigma);
  for(int i=0;i<nval;i++) {
    histo.Fill(values[i]-center);
  }
#ifdef DEBUG
  cout << histo << endl;
#endif
  
  // fill a log array of histo content and take into account Poisson noise
  double *x = new double[nbins];
  double *y = new double[nbins];
  double *w = new double[nbins];
  const float *content = histo.array();
  int validentries=0;
  double sumweight=0;
  for(int i=0;i<nbins;i++) {
    if(content[i]>0) {
      x[validentries]=histo.BinCenter(i);
      y[validentries]=log(content[i]);
      w[validentries]=content[i];
      sumweight+=w[validentries];
      validentries++;
    }
  }
  // if no weights, quit
  if(sumweight==0){
    sigma=0;
    return 0;
  }
  
  
  // now fit this histo with a parabola  y=a*x*x+b*x+c
  double a = -1./(2*sqr(sigma));
  double b = -2.*a*mean;
  
  // first guess of c
  double sum=0;
  for(int i=0;i<validentries;i++) {
    sum += w[i]*(y[i]-(a*sqr(x[i])+b*x[i]));
  }
  double c = sum/sumweight;
#ifdef DEBUG
  cout << "guess -1 " << mean << " " << sigma << " " << c << endl;
#endif



  // iteratively fit the parabola
  Mat A(3,3); // 3 parameters: a,b,c
  Vect B(3); 
  double dRi[3];
  double res,x2,model,nmean,nsigma,chi2,nchi2;
  chi2=-12;
  for(int iteration=0;iteration<20;iteration++) {
    
    // init
    for(int kk=0;kk<3;kk++) {
      B(kk)=0;
      for(int l=0;l<3;l++) {
	A(kk,l)=0;
      }
    }
    
    // fill matrix
    for(int i=0;i<validentries;i++) {
      x2=x[i]*x[i];
      dRi[0]=-x2; // derivative of redisual with respect a
      dRi[1]=-x[i]; // derivative of redisual with respect to b
      dRi[2]=-1; // derivative of redisual with respect to c
      model=a*x2+b*x[i]+c; // the superbe model
      res = y[i]-model; // residual
      w[i] = exp(model); // re-evaluate weight according to model, w=1/sigma2=1/(dN/N)**2=N
      for(int kk=0;kk<3;kk++) {
	B(kk)+=w[i]*res*dRi[kk];
	for(int l=0;l<3;l++) {
	  A(kk,l)+=w[i]*dRi[l]*dRi[kk];
	}
      }
    }
    
    // solve
    int status=cholesky_solve(A,B,"L");
    if(status!=0) {
      return mean;
    }
    
    // apply solution
    a -= B(0);
    b -= B(1);
    c -= B(2);
    
    // check chi2 for possible over shoot
    nchi2=0;
    for(int i=0;i<validentries;i++) {
      nchi2+=w[i]*sqr(y[i]-(x[i]*(a*x[i]+b)+c));
    }
    if(chi2>0 && nchi2>chi2) { // go back
      // cout << "over shoot" << endl;
      a -= B(0)*0.95;
      b -= B(1)*0.95;
      c -= B(2)*0.95;
    }
    
    if (a>=0) {
      cout << "gaussian fit failure (a=" << a << ")" << endl;
      cout << "return a clipped mean" << endl;
      double *newvalues = new double[nval];
      for(int i=0;i<nval;i++)
	newvalues[i]=values[i];
      mean = clipmean(newvalues,nval,sigma,k);
      sigma = -sigma ; // to check error
      return mean;
    }


    chi2=nchi2;
    nsigma=1./sqrt(-2*a);
    nmean=-b/2./a;
#ifdef DEBUG
    cout << "%inc " << (nsigma/sigma-1.) << " " << (nmean/mean-1.) << " ";   
    cout << "guess " << iteration << " " << nmean << " " << nsigma << " " << c << endl;
#endif
    if( fabs(nmean/mean-1.0)<0.0001 && fabs(nsigma/sigma-1.0)<0.0001) 
      break;
    mean=nmean; sigma=nsigma;
  }
  delete [] x;
  delete [] y;
  delete [] w;
   
  mean += center;
  return mean;
}

#ifdef LA_POUBELLE_DE_LHISTOIRE

float new_clipped_median(float *values, int &nval, float &sigma)
{
  
  float mean = 0.0;
  float median = 0.0;
  float sigmat = 0.0;
  
  sort(values,values+nval);

  int middle = int(nval/2); middle--;
  
  float centerValue = values[middle];

  int nSup = nval - middle;
  int nInf = middle;
  
  float *valuesSup = new float[nSup];
  float *valuesInf = new float[nInf];
  for (int i = middle; i < nval; i++) valuesSup[i-middle] = values[i];
  for (int i = 0; i < middle; i++) valuesInf[i] = centerValue-values[i];

  Fmean_median_sigma(valuesSup, nSup, mean, median, sigmat);
  while ( (fabs(mean-median)/sigmat) > 1/sqrt(nSup))
    {
      nSup--;
      Fmean_median_sigma(valuesSup, nSup, mean, median, sigmat);
    }

  Fmean_median_sigma(valuesInf, nInf, mean, median, sigmat);
  while ( (fabs(mean-median)/sigmat) > 1/sqrt(nInf))
    {
      nInf--;
      Fmean_median_sigma(valuesInf, nInf, mean, median, sigmat);
    }

  float *valuesFinal = new float[nInf+nSup];
  for(int i = 0; i < nInf; i++) valuesFinal[i] = centerValue - valuesInf[i];
  for(int i = 0; i < nSup; i++) valuesFinal[i+nInf] = valuesSup[i];

  Fmean_median_sigma(valuesFinal, nInf+nSup, mean, median, sigma);
  nval = nSup+nInf;

  delete [] valuesInf;
  delete [] valuesSup;
  delete [] valuesFinal;

  return median;
}


float clipped_median(float *values, int &nval, const float k, float &sigma)
{
  if (nval < 2 )
    {
      sigma = 0.;
      return values[0];
    }

  sort(values,values+nval);
  
  double mean=0;
  int j,n,n_old;
  int j_middle = int(nval/2);
  j_middle--;
  int j_beg = j_middle-1, j_end = j_middle+1;

  n = nval;
  do
    {
      n_old = n;
      n = 0;
      double meanT=0.0,VarT =0.0;
      for (j=j_end;j>=j_beg;j--)
	{
	  double value = values[j];
	  {
	    n++;
	    meanT += value;
	    VarT += value*value;
	  }
	}
      
      mean = meanT/double(n);

      sigma = (VarT - double(n)*mean*mean)/double(n-1);
      if (sigma>0 && n>1)  sigma = sqrt(sigma);
      else // all values are (almost) identical. we can stop
	sigma = 0.0;
      
      if (sigma>0 && j_beg>0 && fabs(values[j_beg-1]-mean) < k*sigma)
	{
	  n++;
	  j_beg--;
	}
      else if (sigma>0 && j_end+1<nval && fabs(values[j_end+1]-mean) < k*sigma)
	{
	  n++;
	  j_end++;
	}
    }
  while (n!=n_old);
  
  nval = n;
  return mean;
  
}

#endif





