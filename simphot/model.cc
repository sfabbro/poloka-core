#include "model.h"
#include "lightcurvefile.h"
#include "pmstar.h"

#include <cmath> // for floor
static IntPoint nearest_integer_pos(const Point &P)
{
  return IntPoint(int(floor(P.x+0.5)), int(floor(P.y+0.5)));
} 


#include "matvect.h"
#include "imagematch.h"
#include "gtransfo.h"
#include "psfstar.h"
#include "listmatch.h" // for ListMatchCollect
#include "fitsimage.h"

#ifdef STORAGE 

/*********************** MStar is a utility class used in PSFStarListMatch ****************/



//! utility class. May go into dicstar.h
struct MStar : public BaseStar
{
  CountedRef<DicStar> original;
  MStar(const DicStar *S, const string &XName, const string &YName, 
	const string &FluxName) : 
    BaseStar(S->getval(XName), S->getval(YName), 
	     S->getval(FluxName)), original(S) {};
  double getval(const string &S) const { return original->getval(S);} 

  void writen(ostream &s) const;
  string WriteHeader_(ostream & stream, const char*i) const;

};


string MStar::WriteHeader_(ostream & stream, const char*i) const
{
  string format = BaseStar::WriteHeader_(stream,i);
  format += " ";
  format +=original->WriteHeader_(stream,i);
  return format;
}

void MStar::writen(ostream &s) const
{
  BaseStar::writen(s);
  original->writen(s);
}




class MStarList : public StarList<MStar>
 {
 public :
   MStarList() {};
   MStarList(const DicStarList &L, const string &XName, const string &YNname, 
	     const string &FluxName);
};

MStarList::MStarList(const DicStarList &L, const string &XName, 
		     const string &YName, const string &FluxName)
{
  for (DicStarCIterator i = L.begin(); i != L.end(); ++i)
    {
      const DicStar &s = **i;
      push_back(new MStar(&s, XName, YName, FluxName));
    }
}
#endif /* STORAGE */

static StarMatchList *PSFStarsListMatch(const DbImage &DbImage1, 
					const DbImage &DbImage2,
					const Gtransfo *One2Two, 
					const double &Cut)
{
  PSFStarList l1(DbImage1.Dir()+"/psfstars.list");
  PSFStarList l2(DbImage2.Dir()+"/psfstars.list");

  StarMatchList *matches = ListMatchCollect((const BaseStarList &)l1,
					    (const BaseStarList &)l2, 
					    One2Two, Cut);
  cout << " matched " << matches->size() << endl;
  //DEBUG
  //  matches->write("matches.list");

  return matches;
}

/******************************* GTransfoServer class ***********************/
/*! The purpose of the GtransfoServer is to store the transformations 
between images, in order to avoid rematching, when we fit many objects in a single job.
Rather than modifying all the constructors (SimPhotFit, Vignette, Model), in order
to change only the object coordinate (and recompute everything that depends on it), 
it was found easier to store the transformations
*/

struct GtransfoPair
{
  GtransfoRef One2Two;
  GtransfoRef Two2One;
  GtransfoPair(const Gtransfo* T1, const Gtransfo* T2) : One2Two(T1), Two2One(T2) {};
};


class GtransfoServer {
 private :

  map<string,GtransfoPair> trMap;
  typedef map<string,GtransfoPair>::iterator iterator;
  typedef map<string,GtransfoPair>::const_iterator const_iterator;

 public:

  const GtransfoPair* FindTransfos(const DbImage* I1, const DbImage *I2) const
  {
    string key=I1->Dir()+"@"+I2->Dir();
    const_iterator i = trMap.find(key);
    if (i != trMap.end()) return &i->second;
    return NULL;
  }

  void StoreTransfos(const DbImage* I1, const DbImage *I2, const Gtransfo* One2Two, const Gtransfo* Two2One)
  {
    trMap.insert(pair<string,GtransfoPair>(I1->Dir()+"@"+I2->Dir(),GtransfoPair(One2Two,Two2One)));
  }

};

#include "photoratio.h" // for TLSPhotomRatio et al


struct PhotomRatioServer : public map<string,double>
{
  double FindPhotomRatio(const ReducedImage *Ref, const ReducedImage *Current,
			 const GtransfoRef Transfo2Ref, double Cut)
  {
    // Check if we already have it : 
    string key=Current->Name()+"@"+Ref->Name();
    const_iterator it = find(key);
    if (it != end()) return it->second;

    // we do not have it already to compute it ... 
    PSFStarList l1(Current->Dir()+"/psfstars.list");
    // TODO (perhaps) : store the ref list (but file systems have caches)
    PSFStarList l2(Ref->Dir()+"/psfstars.list");

    StarMatchList *matches = ListMatchCollect((const BaseStarList &)l1,
					      (const BaseStarList &)l2, 
					      Transfo2Ref, Cut);
    double sig;
    //TODO : put the sig cut into datacards
    double photomRatio = TLSPhotoRatio(*matches, sig, 5.);
    delete matches;
    (*this)[key] = photomRatio;
    return photomRatio;
  };
};

static PhotomRatioServer ThePhotomRatioServer;


static GtransfoServer TheGtransfoServer;

/********************************* Model Class *****************************/

Model::Model(const LightCurveFile &LCF, const Point &ObjectPos,
	     const double OverSampling) :
  refImage(LCF.GeomRef()), objectPosInImage(ObjectPos), 
  overSampling(OverSampling), refPSF(*refImage,false), 
  useStoredTransfos(LCF.UseStoredTransfos())
{
  //hSizeX = 0;
  //hSizeY = 0;
  pmStar = LCF.FindNearestPmStar(objectPosInImage);
  // if farther than 2 pixels, it is not it.
  if (pmStar && pmStar->Distance(objectPosInImage)>2) pmStar = NULL;
  refMJD = LCF.RefMJD();
  if (pmStar)
    cout << " INFO, using proper motions x,y, pmx, pmy , refdate " 
	 << pmStar->x << ' ' << pmStar->y << ' ' 
	 <<pmStar->pmx << ' ' << pmStar->pmy << ' ' 
	 <<refMJD << endl;
  refPix = nearest_integer_pos(ObjectPos);
  objectPos = objectPosInImage-refPix;
  hasGalaxy = false;
}


#include "sestar.h"


Point Model::ProperMotionOffset(const double &MJDate) const
{
  if (pmStar)
    {
      return Point(pmStar->pmx*(MJDate - refMJD),
		   pmStar->pmy*(MJDate - refMJD));
    }
  return Point(0,0);
}

static std::string transfo_file_name(const ReducedImage *Ref, const ReducedImage *Cur)
{
  return Cur->Dir()+"/transfoTo"+Ref->Name()+".dat";
}

#include "reducedutils.h"

bool Model::FindTransfos(const RImageRef Current, 
			 GtransfoRef &Transfo2Ref,
			 GtransfoRef &TransfoFromRef, 
			 Point &ObjectPosInCurrent,
			 double &PhotomRatio)
{  // should test if we already have stored transfos ... 
  // ... when they are actually storable!

  const GtransfoPair *tfPair = TheGtransfoServer.FindTransfos(Current, refImage);
  if (tfPair)
    {
      Transfo2Ref = tfPair->One2Two;
      TransfoFromRef = tfPair->Two2One;
    }
  else
    {
      if (useStoredTransfos)
	{
	  Transfo2Ref = GtransfoRead(transfo_file_name(refImage,Current));
	  FitsHeader refHead(refImage->FitsName());
	  TransfoFromRef = InversePolyTransfo(*Transfo2Ref, Frame(refHead), 0.0002);
	}
      else
	{
	  // minor tricks to avoid reloading the same ref catalog again and again
	  if (seRef.empty()) LoadForMatch(*refImage, seRef);
	  BaseStarList seCur; LoadForMatch(*Current, seCur);
	  Transfo2Ref = FindTransfo(seCur, seRef, *Current, *refImage);
	  TransfoFromRef = FindTransfo(seRef, seCur, *refImage, *Current);
	  if (!Transfo2Ref) return false;
	}
      TheGtransfoServer.StoreTransfos(Current, refImage, Transfo2Ref, TransfoFromRef);
    }

  // before altering transformations (to reduced coordinates), match
  // PSF stars for photom ratio.
  //TODO put the distance cut (1) into the datacards
  PhotomRatio = ThePhotomRatioServer.FindPhotomRatio(refImage, Current,Transfo2Ref, 1);
  
  ObjectPosInCurrent = TransfoFromRef->apply(objectPosInImage+ProperMotionOffset(Current.ModifiedJulianDate()));
  GtransfoLinShift shift(-refPix.x, -refPix.y);
  Transfo2Ref = GtransfoCompose(&shift, Transfo2Ref);

  GtransfoLin reverseShift(shift.invert());
  TransfoFromRef = GtransfoCompose(TransfoFromRef, &reverseShift);

  //DEBUG
  //cout << " Model::FindTransfos(), expect 0,0 " 
  //     << Transfo2Ref->apply(ObjectPosInCurrent)-objectPos << endl;

  if (overSampling != 1)
    {
      GtransfoLinScale scale(1./overSampling);
      Transfo2Ref = GtransfoCompose(&scale, Transfo2Ref);
      scale.invert();
      TransfoFromRef = GtransfoCompose(TransfoFromRef, &scale);
    }

  return true;
  
}



bool Model::Solve(Mat &A, Vect &B, const string &U_or_L, 
		  const bool FittingGalaxy, const int HalfKernelSize) const
{
  if (B.size() == 0) 
    {
      cout << "  Model::Solve cannot solve without parameters !" << endl;
      return false;
    }
  for (unsigned int j = 0; j < A.SizeY(); ++j)
  for (unsigned int i = 0; i < A.SizeX(); ++i)
    {
      if (isnan(A(i,j)))
	{
	  cout << " Model::Solve : A(" << i << ',' << j << ") = nan" << endl;
	  cout << " giving up " << endl;
	  return false;
	}
    }
  for (unsigned int i = 0; i < B.Size(); ++i)
    {
      if (isnan(B(i)))
	{
	  cout << " Model::Solve : B(" << i << ") = nan" << endl;
	  cout << " giving up " << endl;
	  return false;
	}
    }


#define SMOOTH_EDGES
#ifdef SMOOTH_EDGES
  if (FittingGalaxy)
    {
      //unsmoothed area:
      IntFrame unsmoothed((const IntFrame &) galaxyPixels);

      // because we have convolution AND resampling, we need to add 1      
      int smoothingLength = HalfKernelSize+1;
      cout << " Model::Solve : smoothing galaxy edges over " 
	   << smoothingLength << endl; 
      unsmoothed.CutMargin(smoothingLength);
      int central_index = galaxyPixels.PixelIndex(1,1);
      double weight = A(central_index, central_index)*0.01;
      
      double chi2 = 0;
      for (int j=galaxyPixels.ymin; j < galaxyPixels.ymax; ++j)
	for (int i=galaxyPixels.xmin; i < galaxyPixels.xmax; ++i)
	  {
	    if (unsmoothed.IsInside(i,j)) continue;
	    int in = i;
	    int jn = j;
	    if (i-galaxyPixels.xmin <= smoothingLength) in  = i+1;
	    if (galaxyPixels.xmax -i <= smoothingLength) in  = i-1;
	    if (j-galaxyPixels.ymin <= smoothingLength) jn  = j+1;
	    if (galaxyPixels.ymax -j <= smoothingLength) jn  = j-1;
	    if (in == i && jn == j) abort(); // should never happen!

	    int nindex = galaxyPixels.PixelIndex(in,jn);
	    int index = galaxyPixels.PixelIndex(i,j);
	    
	    A(nindex,nindex) += weight;
	    A(index,index) += weight;
	    A(index, nindex) -= weight;
	    A(nindex, index) -= weight;
	    double res = galaxyPixels(i,j)-galaxyPixels(in,jn); 
	    B(index)-=weight*res;
	    B(nindex)+=weight*res;
	    chi2 += res*res*weight;	    
	  }
      cout << " chi2 smoothing " << chi2 << endl;;
    }// end if (Fitting galaxy)
  return (cholesky_solve(A,B,U_or_L.c_str()) == 0);
#endif
  //#define CUT_EDGES
#ifdef CUT_EDGES
  this is bugged: it also removes fluxes, sky, position!

  IntFrame unsmoothed((const IntFrame &) galaxyPixels);
  // because we have convolution AND resampling, we may need to add 1
  // Let's try to add 1!
  int smoothingLength = HalfKernelSize+1;
  unsmoothed.CutMargin(smoothingLength); 
  PixelBlock center(unsmoothed);
  int nprime = center.Ntot();
  Mat Ap(nprime,nprime);
  Vect Bp(nprime);
  for (int j1=center.ymin; j1 < center.ymax; ++j1)
    for (int i1=center.xmin; i1 < center.xmax; ++i1)
      {
	int k1 = center.PixelIndex(i1,j1);
	int l1 = galaxyPixels.PixelIndex(i1,j1);
	for (int j2=center.ymin; j2 < center.ymax; ++j2)
	  for (int i2=center.xmin; i2 < center.xmax; ++i2)
	    {
	      int k2 = center.PixelIndex(i2,j2);
	      int l2 = galaxyPixels.PixelIndex(i2,j2);
	      Ap(k1,k2) = A(l1,l2);
	    }
	Bp(k1) = B(l1);
      }
  Ap.writeFits("ap.fits");
  bool rc =(cholesky_solve(Ap,Bp,U_or_L.c_str()) == 0);
  if (rc)
  for (int j1=center.ymin; j1 < center.ymax; ++j1)
    for (int i1=center.xmin; i1 < center.xmax; ++i1)
      {
	int k1 = center.PixelIndex(i1,j1);
	int l1 = galaxyPixels.PixelIndex(i1,j1);
	B(l1) = Bp(k1);
      }
  return rc;

#endif
#if (!defined(CUT_EDGES) && !defined(SMOOTH_EDGES))
  return (cholesky_solve(A,B,U_or_L.c_str()) == 0);
#endif
}



#include "imagepsf.h"
#include "resampler.h"


static double ShiftedScalarProduct(const PixelBlock &Psf1,
				   const int Dx1, const int Dy1, 
				   const PixelBlock &Psf2,
				   const int Dx2, const int Dy2)
{
  IntFrame fr1(Psf1);
  fr1 = fr1.Shift(Dx1,Dy1);
  IntFrame fr2(Psf2);
  fr2 = fr2.Shift(Dx2,Dy2);
  IntFrame overlap(fr1*fr2);
  double val = 0;
  PIXEL_LOOP(overlap, i,j)
  {
    val += Psf1(i - Dx1, j - Dy1)*Psf2(i - Dx2, j - Dy2);
  }
  return val;
}



void Model::RefPSFPixels(PixelBlock &PSFPixels) const
{
  ComputePSFPixels(refPSF, objectPos+refPix, refPix, PSFPixels);
}
  



bool Model::FindKernel(const ImagePSF &CurrentPSF,
		       const IntPoint &CurrentIntOffset,
		       const Gtransfo *Vignette2Model,
		       const Gtransfo *Model2Vignette,
		       PixelBlock &Kernel)
{
  // get PSF of "current" image 
  PixelBlock currentPSFPixels;
  Point currentPos(Model2Vignette->apply(objectPos));
  ComputePSFPixels(CurrentPSF, currentPos+CurrentIntOffset,
		   CurrentIntOffset, currentPSFPixels);

  //get PSF of the REF
  PixelBlock refPSFPixels_tmp;
  ComputePSFPixels(refPSF, objectPos+refPix, refPix, refPSFPixels_tmp);
  PixelBlock refPSFPixels((const IntFrame&)refPSFPixels_tmp);
  ResampleImage(refPSFPixels_tmp, Vignette2Model, refPSFPixels);
  //DEBUG
  cout << " psf : ************ " << endl;
  cout << " psf : ref obj pos " << objectPos << endl;
  double mx,my,mx2,my2,mxy;
  refPSFPixels_tmp.Moments(mx,my,mx2,my2,mxy);
  cout << " ref psf avant : xc yc x2c y2c xyc " << mx << ' ' << my << ' ' << mx2 
       << ' ' << my2 << ' ' << mxy << endl;
  cout << " psf : tf centroide " << Model2Vignette->apply(Point(mx,my)) << endl;
  cout << " psf : tf obj pos " << currentPos << endl;
  refPSFPixels.Moments(mx,my,mx2,my2,mxy);
  cout << " ref psf apres : xc yc x2c y2c xyc " << mx << ' ' << my << ' ' << mx2 
       << ' ' << my2 << ' ' << mxy << endl;
  // end DEBUG							  

  //find out the kernel size....
  double xc,yc,x2c,y2c,xyc;
  currentPSFPixels.Moments(xc,yc,x2c,y2c,xyc);
  double xr,yr,x2r,y2r,xyr;
  refPSFPixels.Moments(xr,yr,x2r,y2r,xyr);

  // DEBUG
  cout << " psf : current obj pos " << currentPos << endl;
  cout << " psf : xc yc x2c y2c xyc " << xc << ' ' << yc << ' ' << x2c 
       << ' ' << y2c << ' ' << xyc << endl;
  cout << " psf : xr yr x2r y2r xyr " << xr << ' ' << yr << ' ' << x2r 
       << ' ' << y2r << ' ' << xyr << endl;

  double sigmaDiff = 0.5*(sqrt(fabs(x2r-x2c))+sqrt(fabs(y2r-y2c)));
  //TODO put 3 in the datacards
  int ks = int(floor(2*sigmaDiff+1));
  cout << " kernel half size " << ks << endl;
  Kernel.Allocate(-ks,-ks, ks+1, ks+1);
  Kernel.SetVal(0.);

  // end DEBUG
  int npix = Kernel.Ntot();
  Mat A(npix,npix);
  Vect B(npix);
  PIXEL_LOOP(Kernel, ik1, jk1)
  {
    int index1 = Kernel.PixelIndex(ik1,jk1);
    for (int jk2 = Kernel.ymin; jk2 <= jk1; ++jk2)
      for (int ik2 = Kernel.xmin; ik2 < Kernel.xmax; ++ik2)
	{
	  int index2 = Kernel.PixelIndex(ik2, jk2);
	  A(index2,index1) = ShiftedScalarProduct(refPSFPixels,
						  ik1,jk1,
						  refPSFPixels,
						  ik2,jk2);
	}
    B(index1) =  ShiftedScalarProduct(refPSFPixels,
				      ik1,jk1,
				      currentPSFPixels,0,0);
  }
  if (cholesky_solve(A,B,"U")== 0)
    {
      PixelType *pk = Kernel.Data();
      for (int k = 0; k < npix; ++k) pk[k] = B(k);
      // DEBUG
      cout << " Kernel sum " <<  Kernel.Sum() << endl;
      return true;
    }
  else
    {
      exit(-1);
      return false;
    }
}

void ComputePSFPixels(const ImagePSF &PSF, 
		      const Point &PosInImage,
		      const IntPoint &IntOffset, 
		      PixelBlock &PSFPixels,
		      PixelBlock *PSFPixelsXDer,
		      PixelBlock *PSFPixelsYDer)
{
  /* if caller did not provide any prefered size, use the one over
     which the PSF was computed */
  if (PSFPixels.Ntot() == 0)
    {
      // ordering of parameters different from IntFrame constructor (sorry):
      int hx = PSF.HSizeX();
      int hy = PSF.HSizeY();
      IntFrame psfFrame(-hx,-hy, hx+1, hy+1);
      PSFPixels.Allocate(psfFrame);
      if (PSFPixelsXDer)  PSFPixelsXDer->Allocate(psfFrame);
      if (PSFPixelsYDer)  PSFPixelsYDer->Allocate(psfFrame);
    }
  else // allocate derivatives if needed
    {
      if (PSFPixelsXDer && 
	  (const IntFrame &) *PSFPixelsXDer != (const IntFrame &) PSFPixels)
	PSFPixelsXDer->Allocate(PSFPixels);
      if (PSFPixelsYDer && 
	  (const IntFrame &) *PSFPixelsYDer != (const IntFrame &) PSFPixels)
	PSFPixelsYDer->Allocate(PSFPixels);
    }      
  if (PSFPixelsXDer)
    {
      Vect der(2);
      PIXEL_LOOP(PSFPixels, i, j)
      {
	PSFPixels(i,j) = PSF.PSFValue(PosInImage.x, PosInImage.y, 
				      i + IntOffset.x, j + IntOffset.y,
				      &der);
	(*PSFPixelsXDer)(i,j) = der(0); 
	(*PSFPixelsYDer)(i,j) = der(1); 
      }
    }
  else
    { 
      PIXEL_LOOP(PSFPixels,i,j)
      {
	PSFPixels(i,j) = PSF.PSFValue(PosInImage.x, PosInImage.y, 
				      i + IntOffset.x, j + IntOffset.y, 
				      NULL);
      }
    }
}

