#include <iostream>
#include <fstream>
#include <matvect.h>

#include "lightcurvepoint.h"
#include "lightcurvepoint_dict.h"
#include "objio.h"
#include "typemgr.h"

using namespace std;




int main(int argc, char **argv)
{
  
  // vecteur de flux
  Vect FluxVec;
  {
    Mat m;
    if(m.readFits("vec_sn.fits")!=0)
      return -1;
    FluxVec = m;
  }
  int nflux = FluxVec.Size();
  
  
  // matrice de covariance
  Mat CovarianceMat;  
  if(CovarianceMat.readFits("pmat_sn.fits")!=0)
    return -1;
  Mat FluxCovarianceMat = CovarianceMat.SubBlock(0,nflux-1,0,nflux-1);
  FluxCovarianceMat.Symmetrize("L");
  Mat A;
  A.readFits("nightmat_sn.fits");
  cout << A << endl;
  
  vector<int> suppressedfluxes;
  while(true) {
  
    Mat FluxWeightMat = FluxCovarianceMat;
    FluxWeightMat.SymMatInvert();
    Mat AtWA = A.transposed()*FluxWeightMat*A;
    Mat AtWA_invert = AtWA; AtWA_invert.SymMatInvert();
    Vect flux_per_night = (AtWA_invert*(A.transposed()*FluxWeightMat))*FluxVec;
    cout << "Mean flux per night" << endl;
    cout << flux_per_night << endl;
    cout << "Covariance matrix" << endl;
    cout << AtWA_invert << endl;
  
    // now compute chi2
    Vect B = FluxVec - A*flux_per_night;
    double chi2 = B.transposed()*FluxWeightMat*B;
    cout << "chi2 = " << chi2 << endl;
    
    int ndf = A.SizeY()-A.SizeX();
    cout << "ndf = " << ndf << endl;
    cout << "chi2/ndf = " << chi2/ndf << endl;
    
    if(chi2/ndf<1.5 || suppressedfluxes.size()>=6)
      break;

    int outlier = -1;
    double flux_chi2;
    double chi2_max = 0;
    for(unsigned int iflux = 0 ;iflux<B.Size();++iflux) {
      
      flux_chi2 = pow(B(iflux),2)*FluxWeightMat(iflux,iflux);
      if(flux_chi2>chi2_max) {
	chi2_max = flux_chi2;
	outlier = iflux;
      }
    }
    
    cout << outlier << " " << sqrt(chi2_max) << endl;
    
    if(sqrt(chi2_max)<3.)
      break;
    suppressedfluxes.push_back(outlier);
    // on vire cet outlier de toutes les matrices
    A = A.WithoutRows(outlier,outlier);
    {
      Mat mFluxVec = FluxVec;
      FluxVec = mFluxVec.WithoutRows(outlier,outlier);
      //cout << FluxVec << endl;
    }
    FluxCovarianceMat = FluxCovarianceMat.WithoutRows(outlier,outlier);
    FluxCovarianceMat = FluxCovarianceMat.WithoutColumns(outlier,outlier);
  }
  return 0;

  //list< CountedRef<LightCurvePoint> > lcpoints;
  //obj_input<xmlstream> oi("lc.xml");
  //oi >> lcpoints;
  //oi.close();
  //cout << lcpoints.size() << endl;
  

  return 0;
}

