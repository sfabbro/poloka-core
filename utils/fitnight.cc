#include <math.h>
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
  bool fitsingleflux=false;
  if(argc>1) {
    if(strcmp(argv[1],"-single")==0) {
      fitsingleflux=true;
      cout << "fit of a single constant flux for all points" << endl;
    }
  }

 
  // vecteur de flux
  Vect FluxVec;
  {
    Mat m;
    if(m.readFits("vec_sn.fits")!=0)
      return -1;
    FluxVec = m;
  }
  int nflux = FluxVec.Size();
  //cout << FluxVec << endl;
  
  // matrice de covariance
  Mat CovarianceMat;  
  if(CovarianceMat.readFits("pmat_sn.fits")!=0)
    return -1;
  //cout << "CovarianceMat" << endl;
  //cout << CovarianceMat << endl;

  Mat FluxCovarianceMat = CovarianceMat.SubBlock(0,nflux-1,0,nflux-1);
  FluxCovarianceMat.Symmetrize("L");
      
  //cout << FluxCovarianceMat << endl;
  Mat A;
  if(!fitsingleflux) {
    A.readFits("nightmat_sn.fits");
  }else{  
    A.allocate(1,nflux);
    for(int i=0;i<nflux;++i)
      A(0,i)=1;
  }
  Mat Abis = A; // we save a copy for output 
  cout << A << endl;
  
  vector<int> suppressedfluxes;
  Vect flux_per_night;
  Mat AtWA_invert;
  while(true) {
  
    Mat FluxWeightMat = FluxCovarianceMat;
    FluxWeightMat.SymMatInvert();
    //cout << "FluxWeightMat" << endl;
    //cout << FluxWeightMat << endl;
    Mat AtWA = A.transposed()*FluxWeightMat*A;
    AtWA_invert = AtWA; AtWA_invert.SymMatInvert();
    flux_per_night = (AtWA_invert*(A.transposed()*FluxWeightMat))*FluxVec;
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
  
  // OUTPUT
  // ============================================================================================


  // save these results in ASCII files
  // ... TODO ...
  {
    ofstream st("flux_per_expo_covmat.dat");
    st << FluxCovarianceMat;
    st.close();
    FluxCovarianceMat.writeFits("flux_per_expo_covmat.fits");
    ofstream st2("flux_per_night_covmat.dat");
    st2 <<  AtWA_invert;
    st2.close();
    AtWA_invert.writeFits("flux_per_night_covmat.fits");
  }
  

  // now we read the lightcurve point list to get julian dates and zero point
  vector< CountedRef<LightCurvePoint> > lcpoints;
  {
    vector< CountedRef<LightCurvePoint> > lcpoints_all;
    obj_input<xmlstream> oi("lc.xml");
    oi >> lcpoints_all;
    oi.close();
    
    // no keep only points with flux!=0 (i.e. fitted points)
    vector< CountedRef<LightCurvePoint> >::iterator it = lcpoints_all.begin();
    for(;it!=lcpoints_all.end();++it) {
      cout << (*it)->julianday << " " << (*it)->flux;
      if((*it)->flux!=0) {
	lcpoints.push_back(*it);
      }
      cout << endl;
    }
  }
  

  cout << lcpoints.size() << endl;
  for(unsigned int expo=0;expo<lcpoints.size();++expo) {
    	cout 
	     << expo << " " 
	     << lcpoints[expo]->flux << " " 
	     << lcpoints[expo]->julianday-2452854.0 << endl;
  }

  // get zero point 
  double zp = lcpoints[0]->zeropoint;
  cout << "zp=" << zp << endl;


  ofstream outputlc("lc_per_night.dat");
  outputlc << "# jd : \n"
	   << "# flux : \n"  
	   << "# eflux : \n"
	   << "# mag : using elixir zero point = " << zp << "\n"
	   << "# emag_minus : useful for drawing\n"
	   << "# emag_plus : useful for drawing\n"
	   << "# zeropoint : elixir zp\n";
  outputlc << "# end \n";
  
  vector< CountedRef<LightCurvePoint> > newlcpoints;
  double jd;
  int nexpo;
  for(unsigned int night = 0; night < flux_per_night.Size(); ++ night) {
    CountedRef<LightCurvePoint> newpoint = new LightCurvePoint();
    newpoint->flux = flux_per_night(night);
    newpoint->eflux = sqrt(FluxCovarianceMat(night,night));
    newpoint->computemag(zp);
    // now get julian day (mean of all exposures)
    jd=0;
    nexpo=0;
    for(unsigned int expo=0;expo<Abis.SizeY();++expo) {
      if(Abis(night,expo)>0.5) {
	jd += lcpoints[expo]->julianday;
	cout << night << " " 
	     << expo << " " 
	     << lcpoints[expo]->julianday-2452854.0 << endl;
	nexpo++;
      }
    }
    newpoint->julianday = jd/nexpo;
    
    outputlc << (*newpoint) << endl;
    newlcpoints.push_back(newpoint);
  }
  outputlc.close();
  
  // save it also in xml
  obj_output<xmlstream> oo("lc_per_night.xml");
  oo << newlcpoints;
  oo.close();

  
  // done !!

  return 0;
}

