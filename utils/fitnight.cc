#include <iostream>
#include <fstream>
#include <matvect.h>

using namespace std;


int main(int argc, char **argv)
{
  int nflux = 32;
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
  
  Vect FluxVec;
  {
    Mat m;
    if(m.readFits("vec_sn.fits")!=0)
      return -1;
    FluxVec = m;
  }
  Mat A;
  A.readFits("nightmat_sn.fits");
  cout << A << endl;
  
  //cout << FluxVec << endl;
  //Mat AtWA = A.transposed()*FluxWeightMat*A;
  //cout << AtWA << endl;
  
  return 0;
}

