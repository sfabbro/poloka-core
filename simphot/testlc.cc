#include "pixelblock.h"
#include "gtransfo.h"
#include "resampler.h"
#include "array4d.h"
#include "vignette.h"

#include <iostream>
#include "matvect.h"

#ifdef STORAGE
void ConvolveModel(const PixelBlock &M, const Array4D &Coeffs, 
		   PixelBlock &Result)
{
  if ((const IntFrame &) Result != (const IntFrame &) Coeffs)
    {
      cout << " size inconsistency in  ConvolveModel " << endl;
      abort();
    }
  for (int b=Result.ymin; b < Result.ymax; ++b)
    for (int a = Result.xmin; a< Result.xmax; ++a)
      {
	const CoeffBlock &coeffs = Coeffs(a,b);
	double val = 0;
	for (int j=coeffs.ymin; j < coeffs.ymax; ++j)
	  for (int i = coeffs.xmin; i < coeffs.xmax; ++i)
	    {
	      val += M(i,j)*coeffs(i,j);
	    }
	Result(a,b) = val;
      }
}





int main()
{


  int size = 7;
  PixelBlock model(-size, -size, size+1, size+1);
  model.SetVal(1.);
  //  model(-1,-1) = 1;
  // model(3,3) = 1;
  
  GtransfoLinShift shift(0.3,0.3);
  GtransfoLin invShift(shift);
  invShift.invert();
 
  IntFrame shiftFrame(model);
  shiftFrame.CutMargin(ResamplerBoundarySize());
  PixelBlock resamp(shiftFrame);
  ResampleImage(model, &invShift, resamp);


  int hks =0;
  PixelBlock kernel(-hks,-hks, hks+1, hks+1);
  kernel(0,0) = 1;
  //  GaussianFill(kernel, 0.5, 2.0);
  cout << " kernel integral " << kernel.Sum() << endl;
  PixelBlock out;

  ConvolveImage(resamp, kernel, out);


  kernel.WriteFits("kernel.fits");
  out.WriteFits("out.fits");
  cout << " model bounds " << model << " resamp bounds " << resamp 
       << " resamp+conv bounds " << out << endl;


  Array4D resampDer(resamp);
  ResamplerComputeDerivatives(&invShift, resampDer);
  Array4D resamp_transposed;
  resampDer.Transpose(resamp_transposed);
  cout << " resamp_transposed" << resamp_transposed << endl;
  resamp_transposed.Convolve(kernel,resampDer);
  Array4D der;
  resampDer.Transpose(der);
  PixelBlock resultDer((const IntFrame &)der);
  ConvolveModel(model, der, resultDer);
  resultDer.WriteFits("der.fits");
  return EXIT_SUCCESS;
}
#endif


int main(int nargs, char **args)
{
  if (nargs<3)
    {
      cout << args[0] << ' ' << " <mat.fits> <rank> " << endl
	   << " diagonalise la sous-matrice 0:rank de mat.fits " << endl;
      exit(-1);
    }
  string name(args[1]);
  int si = atoi(args[2]);
      

  Mat Ai;
  Ai.readFits(name);
  Mat A(si,si);
  for (int j=0; j < si; ++j)
  for (int i=0; i < si; ++i)
    {
      A(i,j) = Ai(i,j);
      if (j<i) A(i,j) = Ai(j,i);
    }

  A.writeFits("acut.fits");
  Mat vectors(A);
  Vect values(si);
  DiagonalizeRealSymmetricMatrix(A, vectors,values);

  //  Vect test(si);
  // for (int i=0; i < si; ++i) test(i) = vectors(i,20);
  // cout << A*test-test*values(20);

  cout << " valeurs propres " << endl << values << endl;

  cout << " les vecteurs propres apparaissent HORIZONTALEMENT dans eigenvectors.fits" << endl;
  vectors.writeFits("eigenvectors.fits");


  return EXIT_SUCCESS;

}
