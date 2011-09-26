
#include "vignette.h"
#include "lightcurvefile.h"
#include "objecttofit.h"

#include "imagematch.h"


#include "gtransfo.h"
#include "imagepsf.h"
#include "imagepsfserver.h"


#include <cmath> // for floor
static IntPoint nearest_integer_pos(const Point &P)
{
  return IntPoint(int(floor(P.x+0.5)), int(floor(P.y+0.5)));
} 

static double sq(const double &x) { return x*x;}


static int nearest_integer(const double &X)
{
  return int(floor(X+0.5));
}

// should go in Frame class:
static IntFrame NearestIntegerFrame(const Frame &In)
{
  return IntFrame(nearest_integer(In.xMin), 
		  nearest_integer(In.yMin),
		  nearest_integer(In.xMax),
		  nearest_integer(In.yMax));
}


#include "simphotfit.h"
#include "model.h"
#include "resampler.h"
#include "array4d.h"


Vignette::Vignette(SimPhotFit &SPF, const RImageRef &Current):
  ri(Current), simPhotFit(SPF), photomRatio(1.), 
  flux(0), sky(0)
{

  //  mjd = ri->ModifiedJulianDate();
  //  exptime = ri->Exposure();
  //  seeing = ri->GFSeeing();
  //  sesky = ri->BackLevelNoSub();
  //  sigsky = ri->SigmaBack();
  size_n_seeing = SPF.vignette_size_n_seeing;
  const ObjectToFit &obj =  SPF.ObjToFit();
  double mjd = MJD();
  mightFitFlux = (mjd >=obj.JdMin() && mjd <= obj.JdMax());
  imagePSF = FindImagePSF(Current);
}


/* When we want to fit several objects using the same Vignette,
   we have to store the 2 geometrical transformations and the photometric ratio
   to the reference. One way would be to store them into some data structure
   in SimPhotFit, so that SimPhotFit::FindTransfos is fast, and we can just
   forget the earlier Vignette's. 
   Caveat : the ImagePSF loaded in each Vignette will be lost.
*/


bool Vignette::SetGeomTransfos()
{
  //TODO : compute photom ratio.
  // TODO : calculer le decalage entre position posInImage et trFromRef(simPhotFit.ObjectPos());

  GtransfoRef trFromRef;
  GtransfoRef tr2Ref;
  Point posInImage;

  if (!simPhotFit.FindTransfos(ri, tr2Ref, trFromRef, posInImage, photomRatio))
    return false;

  cout << " photom ratio for image " << Name() << " : " << photomRatio << endl;

  intPos = nearest_integer_pos(posInImage);
  // reduced coordinates
  posInStamp = posInImage-intPos;

  // store transformations from and to reduced coordinates
  GtransfoLinShift shiftHere(intPos.x, intPos.y);
  vignette2Model = GtransfoCompose(tr2Ref,&shiftHere);
  GtransfoLin reverseShift(shiftHere.invert());
  model2Vignette = GtransfoCompose(&reverseShift, trFromRef);


  //TODO : put "4" into the datacards : put in size_n_seeing now
  int radius = nearest_integer(size_n_seeing * Seeing()+1);
  IntFrame frame(0,0,1,1); // single pixel frame
  //  frame.CutMargin(-radius); // right size
  frame.CutMargin(-(imagePSF->HSizeX()), -(imagePSF->HSizeY())); // borrow stamp size to PSF
  stampLimits = frame;
  convolvedStampLimits = stampLimits; // until SetKernel eventually updates
  //DEBUG
  cout << ri->Name() << " stampLimits " << stampLimits << endl;



  // actually load pixels
  return ReadPixels();
}

// I did not put this constructor in Frame to avoid "cross talk"
static Frame FloatFrame(const IntFrame &IF)
{
  return Frame(IF.xmin, IF.ymin, IF.xmax, IF.ymax);
}


#include "frame.h"
#include "imageutils.h" // for ApplyTransfo(Frame ...)

void Vignette::ComputeModelLimits(Frame &ModelFrame) const
{
  Frame currentFrame= FloatFrame(convolvedStampLimits);
  // TODO : if oversampling, we have to do something here: enlarge a little bit
  // current frame.
  ModelFrame = ApplyTransfo(currentFrame, *vignette2Model, LargeFrame);
  cerr << " #ComputeModelLimits : convolvedStampLimits=" << convolvedStampLimits << " apres transfo : " << ModelFrame.Nx() << " ResamplerBoundarySize = " << ResamplerBoundarySize() << endl ;
  // now add the resampling overhead:
  // Frame::CutMargin increases the frame size if argument is <0
  ModelFrame.CutMargin(-ResamplerBoundarySize());
  // DEBUG
  cout << ri->Name() << " Model limits " << ModelFrame << endl;

}

#include "matvect.h"

void Vignette::ComputeGalaxyDerivatives(Array4D &GalaxyDerivatives)
{
  // without kernel : 
  // ResamplerComputeDerivatives(vignette2Model, GalaxyDerivatives);
  Array4D resamp(convolvedStampLimits);
  ResamplerComputeDerivatives(vignette2Model, resamp);

  Array4D resamp_transposed;
  resamp.Transpose(resamp_transposed);

  resamp_transposed.Convolve(kernel,resamp);
  resamp.Transpose(GalaxyDerivatives);
}

void Vignette::UpdateResiduals()
{
  if (imagePix.Ntot() == 0) ReadPixels();
  residuals = imagePix;
  if (simPhotFit.HasGalaxy()) // subtract it
    {
      PixelBlock resampModel(convolvedStampLimits);
      ResampleImage(simPhotFit.GalaxyPixels(), vignette2Model, resampModel);
      PixelBlock convolvedGal;
      ConvolveImage(resampModel,kernel, convolvedGal);
      residuals -= convolvedGal;
    }
  if (flux != 0)    // subtract point source from residuals
    {
      // DEBUG
      //  Point oldPos = posInStamp;
      posInStamp = model2Vignette->apply(simPhotFit.ObjectPos(MJD()));
      //      cout << " obj pos in image " << Name() << " old " << oldPos << " new " << posInStamp 
      //	   << 	" delta " << posInStamp-oldPos << endl;
      // DEBUG
      //      cout << " verif " << vignette2Model->apply(posInStamp) - simPhotFit.ObjectPos(MJD()) << endl;

      PixelBlock psf(residuals);// borrow the frame
      GetPSF(psf);
      psf *= flux;
      residuals -= psf;
    }
  residuals -= sky;   // subtract sky

  // compute chi2 and number of chi2 terms this Vignette provides
  nterms = 0;
  chi2 = 0;
  PIXEL_LOOP(weightPix, a, b)
    {
      double w = weightPix(a,b);
      chi2 += sq(residuals(a,b))*w;
      if (w) nterms += 1;
    }
}

int Vignette::KillOutliers(const double NSigCut)
{    
  double sw = 0;
  double sf = 0;
  double sf2 = 0;
  PixelType *pw = weightPix.begin();
  const PixelType *pr = residuals.begin();
  const PixelType *pend = residuals.end();
  for ( ; pr<pend; ++pr, ++pw)
    {
      sw += *pw; sf += *pw * (*pr); sf2 += *pw* sq(*pr);   
    }
  if (sw == 0) return 0;
  double mean = sf/sw;
  double sigma = sf2/sw-sq(mean);
  if (sigma<=0) return 0; // might happen with a single pixel.
  sigma = sqrt(sigma);
  double cuth = mean + NSigCut * sigma;
  double cutl = mean - NSigCut * sigma;
  int outCount = 0;
  for (pr = residuals.begin(), pw = weightPix.begin() ; pr<pend; ++pr, ++pw)
    {
      if (*pr > cuth  || *pr < cutl ) {*pw = 0; outCount++;}
    }
  return outCount;
}
						       

bool Vignette::GetPSF(PixelBlock &PSFPixels, PixelBlock *PSFDerX,
		      PixelBlock *PSFDerY) const
{
  Point posInImage = posInStamp+ intPos;
  ComputePSFPixels(*imagePSF, posInImage, intPos,
		   PSFPixels,PSFDerX,PSFDerY);
  if (PSFDerX && PSFDerY)
    {
      // operate the transformation to express the derivative w.r.t 
      // model coordinates rather than this image coordinates.
      GtransfoLin der;
      model2Vignette->Derivative(posInStamp, der, 1.);
      double dx,dy;
      for (int j = PSFPixels.ymin; j < PSFPixels.ymax; ++j)
      for (int i = PSFPixels.xmin; i < PSFPixels.xmax; ++i)
	{
	  dx = (*PSFDerX)(i,j)*der.A11() + (*PSFDerY)(i,j)*der.A12(); 
	  dy = (*PSFDerX)(i,j)*der.A21() + (*PSFDerY)(i,j)*der.A22(); 
	  (*PSFDerX)(i,j) = dx;
	  (*PSFDerY)(i,j) = dy;
	}
    }
  return true;
}

	  
int Vignette::CanDo() const
{
  if (weightPix.Sum() == 0) return 0;
  // sum weights at the center of the stamp 
  IntFrame sumFrame;
  sumFrame.CutMargin(-3);
  double sumw = 0;
  sumFrame *= (const IntFrame &)weightPix;
  for (int j=sumFrame.ymin; j < sumFrame.ymax; ++j)
    for (int i = sumFrame.xmin; i < sumFrame.xmax; ++i)
      sumw += weightPix(i,j);
  int toDo = FIT_GALAXY + FIT_SKY;
  if (sumw > 0 && mightFitFlux) toDo += (FIT_FLUX + FIT_POS);
  return toDo;
}


void Vignette::FillAAndB(Mat &A, Vect &B, const int ToDo)
{
  //DEBUG
  //  cout << " Vignette::FillAAndB " << ri->Name() << endl;


  
  /* structure of the routine :
     there are 4 type of parameters (galaxy/flux/sky/position)
     and hence 10 types of terms to plug into the matrix 
     (gal-gal, gal-flux, gal-sky, .... ,pos-pos).

     The routine makes use of the ordering of the parameters (galaxy <
     flux < sky < position) in order to fill only the upper half of the
     matrix.  

     The filling of the B vector (the right hand side) is done
     together with the diagonal blocks.


     To validate the code, we just checked that when the problem
     is linear, (i.e. fixed position), a second iteration does not change
     anything. We also checked that when at the minimum, iterating the fit
     with any parameter class frozen does not change the parameters.

  */

  // collect indices to address A and B. 
  /*  If ToDo does not indicate that we should use a certainindex, 
      it is set to -1, and a (graceful) crash should follow if we use 
      the senseless index.
  */
  int fluxIndex = simPhotFit.FluxIndex(this);
  int skyIndex = simPhotFit.SkyIndex(this);
  int posIndex = simPhotFit.PosIndex();

  const PixelBlock &galaxy = simPhotFit.GalaxyPixels(); 
  Array4D modelDer(stampLimits);
  PixelBlock psf;
  PixelBlock psfDx;
  PixelBlock psfDy;
  if ((ToDo & (FIT_FLUX | FIT_POS)))
    {
      /* these "allocate" are not scrictly necessary: the PSF reader
	 "Allocate"'s at the PSF modelled size. If we turn to that,
	 we have to check the pixel loops involving the PSF. */
      psf.Allocate(stampLimits);
      psfDx.Allocate(stampLimits);
      psfDy.Allocate(stampLimits);
      GetPSF(psf, &psfDx, &psfDy);
    }
  


  if (ToDo & FIT_GALAXY) 
    {
      ComputeGalaxyDerivatives(modelDer);

      // Gal-Gal terms
      PIXEL_LOOP(weightPix, a, b)
	{
	  double w = weightPix(a,b);
	  if (w == 0) continue;
	  double res = residuals(a,b);
	  const CoeffBlock &block = modelDer(a,b);
	  int niblock = block.xmax-block.xmin;
	  for (int j1=block.ymin; j1 < block.ymax; ++j1)
	    for (int i1=block.xmin; i1 < block.xmax; ++i1)
	      {
		int mati1 = galaxy.PixelIndex(i1,j1); 
		double hi1 = block(i1,j1)*w;
#define FASTWAY
#ifndef FASTWAY
		for (int j2=block.ymin; j2 <= j1 ; ++j2)
		  {
		    for (int i2=block.xmin; i2 < block.xmax; ++i2)
		      {
			int mati2 = galaxy.PixelIndex(i2,j2);
			A(mati2,mati1) += hi1*double(block(i2,j2));
		      }
		  }
#else
		for (int j2=block.ymin; j2 <= j1 ; ++j2)
		  {
		    int imin2 = block.xmin;
		    const CoeffType *pb = &block(imin2, j2);
		    double *pa = &A(galaxy.PixelIndex(imin2,j2), mati1);
		    for (int k = niblock; k ; k--)
		      {
			*pa += hi1 * double(*pb);
			pa++; pb++;
		      }
		  }
#endif
		// B filling - Gal 
		B(mati1) += res*hi1;
	      }
	  }
      // End of Gal-Gal terms
      // Gal-Flux terms
      if (ToDo & FIT_FLUX)
	{
	  PIXEL_LOOP(weightPix, a, b)
	  {
	    double w = weightPix(a,b);
	    if (w == 0) continue;
	    double psfVal = w*psf(a,b);
	    const CoeffBlock &block = modelDer(a,b);
	    PIXEL_LOOP(block, i, j)
	    {
	      int mati = galaxy.PixelIndex(i,j); 
	      A(mati,fluxIndex) += block(i,j)*psfVal;
	    }
	  }
	}
      //End of Gal-Flux
      // Gal-sky
      if (ToDo & FIT_SKY)
	{
	  PIXEL_LOOP(weightPix, a, b)
	  {
	    double w = weightPix(a,b);
	    if (w == 0) continue;
	    const CoeffBlock &block = modelDer(a,b);
	    PIXEL_LOOP(block, i, j)
	    {
	      int mati = galaxy.PixelIndex(i,j); 
	      A(mati,skyIndex) += w*block.at(i,j);
	    }
	  }
	}
      // end of Gal -sky
      // Gal-pos
      if (ToDo & FIT_POS)
	{
	  PIXEL_LOOP(weightPix, a, b)
	  {
	    double w = weightPix(a,b);
	    if (w == 0) continue;
	    double dx = w*flux*psfDx(a,b);
	    double dy = w*flux*psfDy(a,b);
	    const CoeffBlock &block = modelDer(a,b);
	    PIXEL_LOOP(block, i, j)
	    {
	      int mati = galaxy.PixelIndex(i,j);
	      A(mati,posIndex)   += dx*block(i,j);
	      A(mati,posIndex+1) += dy*block(i,j);
	    }
	  }
	}// End of Gal-Pos      
    } // End of if FIT_GALAXY
  if (ToDo & FIT_SKY)
    {
      // Sky-Sky
      A(skyIndex,skyIndex) += weightPix.Sum();
      B(skyIndex) += ScalProd(weightPix,residuals);
      //cout << "ScalProd(weightPix,residuals) weightPix residuals " << ScalProd(weightPix,residuals) << " " << weightPix << " " << residuals << endl ; 
 // B term - Sky
      // DEBUG
      if (isnan(B(skyIndex))) 
	{
	  cout << "B(skyIndex)=" << B(skyIndex)	<< " " << skyIndex << endl;   
	  abort();
	}
     // Sky - Flux
      if (ToDo & FIT_FLUX)
	A(skyIndex, fluxIndex) += ScalProd(psf,weightPix);
      // Sky - Pos
      if (ToDo & FIT_POS)
	{
	  A(skyIndex, posIndex)   += flux*ScalProd(weightPix,psfDx);
	  A(skyIndex, posIndex+1) += flux*ScalProd(weightPix,psfDy);
	}
    } // end of FIT_SKY

	
  if (ToDo & FIT_FLUX)
    {
      // Flux- Flux terms.
      double suma  = 0;
      double sumb  = 0;
      PIXEL_LOOP(weightPix, a, b)
      {
	double psfVal = psf(a,b);
	double w = weightPix(a,b);
	suma += w*psfVal*psfVal;
	sumb += w*psfVal*residuals(a,b);
      }
      A(fluxIndex, fluxIndex) += suma;
      B(fluxIndex) += sumb; // B term - Flux
      // End of Flux- Flux
      // Flux - Pos
      if (ToDo & FIT_POS)
	{
	  double sumx = 0;
	  double sumy = 0;

	  PIXEL_LOOP(weightPix, a, b)
	  {
	    double psfVal = psf(a,b);
	    double w = weightPix(a,b);
	    sumx += w*psfVal*psfDx(a,b);
	    sumy += w*psfVal*psfDy(a,b);
	  }
	  A(fluxIndex,posIndex)   += sumx;
	  A(fluxIndex,posIndex+1) += sumy;
	} // End of Flux-Pos
    } // End of FIT_FLUX
  if (ToDo & FIT_POS)
    {
      double sumxx  = 0;
      double sumxy  = 0;
      double sumyy  = 0;
      double sumbx  = 0;
      double sumby  = 0;
      PIXEL_LOOP(weightPix, a, b)
      {
	double w = weightPix(a,b);
	double dx = flux*psfDx(a,b);
	double dy = flux*psfDy(a,b);
	sumxx += w*dx*dx;
	sumxy += w*dx*dy;
	sumyy += w*dy*dy;
	double res = residuals(a,b);
	sumbx += w*dx*res;
	sumby += w*dy*res;
      }
      A(posIndex  , posIndex  ) += sumxx;
      A(posIndex  , posIndex+1) += sumxy;
      A(posIndex+1, posIndex+1) += sumyy;
      B(posIndex  ) += sumbx;
      B(posIndex+1) += sumby;
    } // end of FIT_POS
}



#include "imagepsf.h"

bool Vignette::SetKernel()
{
  simPhotFit.FindKernel(*imagePSF, intPos, 
			vignette2Model, model2Vignette, kernel);
  //DEBUG

  convolvedStampLimits = stampLimits;
  convolvedStampLimits.CutMargin(-HalfKernelSize());  
  cerr << " #SetKernel : convolvedStampLimits = " << convolvedStampLimits << endl ;

  return true;
}

bool Vignette::ReadPixels()
{
  //TODO : check that stampLimits is entirely inside the fitsimage, if not,
  // update stampLimits
  imagePix.Allocate(stampLimits);
  if (!imagePix.ReadFromFits(ri->FitsName(), int (intPos.x), int (intPos.y)))
    return false;
  imagePix *= 1./photomRatio;

  weightPix.Allocate(stampLimits);
  if (!weightPix.ReadFromFits(ri->FitsWeightName(), int (intPos.x), int (intPos.y))) return false;
  weightPix *= photomRatio*photomRatio;

  PixelBlock saturPix ;
  saturPix.Allocate(stampLimits);
  if (saturPix.ReadFromFits(ri->FitsSaturName(), int (intPos.x), int (intPos.y))) // il y a une carte de satur
    {
      PixelType *pw = weightPix.begin();
      PixelType *ps = saturPix.begin();
      double sum = 0 ;
      for( ; pw < weightPix.end() && ps < saturPix.end() ; pw++, ps++)
	{
	  *pw *= (1-*ps) ;
	  sum += *ps ;
	}

      has_saturated_pixels=(sum>0);   
      n_saturated_pixels=sum;    
    }
  else
    {
      // on dit rien ???
      //return false ;
    }



  // print de DEBUG des vignettes
  //weightPix.WriteFits(Name()+"_weightpix.fits");
  //imagePix.WriteFits(Name()+"_imagepix.fits");
  
  return true;

}


void Vignette::Write(const string &Directory) const
{
  string genericName = Directory+ri->Name();
  imagePix.WriteFits(genericName+".fits");  
  //imagePSF.WriteFits(genericName+".psf.fits");
  weightPix.WriteFits(genericName+".weight.fits");
  residuals.WriteFits(genericName+".res.fits");
  kernel.WriteFits(genericName+".kernel.fits");
  // write residuals of refPSf*kernel to thisPSF.
  PixelBlock refPSF;
  simPhotFit.RefPSFPixels(refPSF);
  PixelBlock resampRefPSF(refPSF);
  ResampleImage(refPSF, vignette2Model, resampRefPSF);
  PixelBlock conv;
  ConvolveImage(resampRefPSF, kernel, conv);
  PixelBlock thisPSF(conv); // borrow the Frame
  GetPSF(thisPSF);
  thisPSF.WriteFits(genericName+".psfconv.fits");
  PixelBlock diff = thisPSF-conv;
  diff.WriteFits(genericName+".psfconvres.fits");
}




#ifdef STORAGE
/* this routine is probably temporary, because derivatives are in principle
   only necessary in FillAAndB. In UpdateResiduals, we should use "direct"
   resampling and convolution rather than use these derivatives, for sake
   of efficiency. Before such tools actually exist, we just use these 
   derivatives in UpdateResiduals */
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
#ifndef FAST_WAY
	for (int j=coeffs.ymin; j < coeffs.ymax; ++j)
	  for (int i = coeffs.xmin; i < coeffs.xmax; ++i)
	    {
	      val += M(i,j)*coeffs(i,j);
	    }
#else
	// unroll the "i" loop
#endif
	Result(a,b) = val;
      }
}
#endif /* STORAGE */
  



