#include <iostream>
#include <fstream>
#include <matvect.h>

using namespace std;


int main(int argc, char **argv)
{
  
  Vect FluxVec;
  {
    Mat m;
    if(m.readFits("vec_sn.fits")!=0)
      return -1;
    FluxVec = m;
  }
  int nflux = FluxVec.Size();

  
  Mat CovarianceMat;  
  if(CovarianceMat.readFits("pmat_sn.fits")!=0)
    return -1;
  Mat FluxCovarianceMat = CovarianceMat.SubBlock(0,nflux-1,0,nflux-1);
  FluxCovarianceMat.Symmetrize("L");
  Mat FluxWeightMat = FluxCovarianceMat;
  FluxWeightMat.SymMatInvert();
  // test
  //Mat res = FluxWeightMat*FluxCovarianceMat;
  //res.writeFits("toto.fits");
  //cout << res.SizeX() << " " << res.SizeY() << endl;
  
 
  Mat A;
  A.readFits("nightmat_sn.fits");
  cout << A << endl;
  
  cout << FluxWeightMat.SizeX() << " " << FluxWeightMat.SizeY() << endl;
  cout << A.SizeX() << " " << A.SizeY() << endl;
  Mat AtWA = A.transposed()*FluxWeightMat*A;
  cout << AtWA << endl;
  Mat AtWA_invert = AtWA; AtWA_invert.SymMatInvert();
  
  Vect flux_per_night = (AtWA_invert*(A.transposed()*FluxWeightMat))*FluxVec;
  cout << "Mean flux per night" << endl;
  cout << flux_per_night << endl;
  cout << "Covariance matrix" << endl;
  cout << AtWA_invert << endl;
  
  return 0;
}

