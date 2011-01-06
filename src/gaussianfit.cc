#include <iostream>
#include "matvect.h"
#include "histo1d.h"
#include "vutils.h"

using namespace std;

static double sqr(double x) {return x*x;};

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
  if (sigma < 1.e-30 ){
    mean =-1;
    sigma =-1 ;
    return(-1) ;
  }
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
    if ( (-2*a) < 1.e-30 ){
      mean =-1;
      sigma =-1 ;
      return(-1) ; 
    }
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
