#include <iostream>
#include <string>


#include "transformedimage.h"
#include "gtransfo.h"
#include "fitsimage.h"
#include "wcsutils.h"
#include "fileutils.h"
#include "imageutils.h"
#include "reducedutils.h"

string TransformedName(const string &ToTransform, const string &Ref)
{
  return "T_"+Ref+ToTransform;
}



/****************** ImageGtransfo ***********************/

#include "imagematch.h"



ImageGtransfo::ImageGtransfo(const Gtransfo *TransfoFromRef, 
			     const Gtransfo *TransfoToRef, 
			     const Frame &OutputImageSize, 
			     const string &GeomRefName)
{
  geomRefName = GeomRefName;
  transfoFromRef = (TransfoFromRef) ? TransfoFromRef->Clone() : NULL ;
  transfoToRef = (TransfoToRef) ? TransfoToRef->Clone() : NULL;
  outputImageSize = OutputImageSize; 
  scaleFactor = 1.;
  if (transfoFromRef) 
    {scaleFactor = 1./sqrt(fabs(transfoFromRef->Jacobian(outputImageSize.Nx()/2, outputImageSize.Ny()/2)));}
}

const GtransfoRef ImageGtransfo::FromRef() const
{
  return transfoFromRef;
}


ImageTransfo* ImageGtransfo::Clone() const
{
  return new ImageGtransfo(transfoFromRef, transfoToRef, outputImageSize, geomRefName);
  }

ImageGtransfo::ImageGtransfo(const ReducedImage &Ref, const ReducedImage& ToAlign,float min_match_ratio)
{
  //transfoFromRef = NULL;
  //transfoToRef = NULL;
  geomRefName = Ref.Name();
  if (!ImageListMatch(Ref, ToAlign, transfoFromRef, transfoToRef, min_match_ratio))
    {
      cerr << " Could not match lists from " << Ref.Name() << " and " <<  ToAlign.Name() << endl;
    }
  outputImageSize = Ref.PhysicalSize();
  // use the same transfo as the one used to actually transform the image.
  if (transfoFromRef)
    {scaleFactor = 1./sqrt(fabs(transfoFromRef->Jacobian(outputImageSize.Nx()/2, outputImageSize.Ny()/2)));}
  else scaleFactor = 1.;
  cout << " Scale factor for rebinning = " << scaleFactor << endl;
  double pixRatio = ToAlign.PixelSize()/Ref.PixelSize();
  if ( fabs( pixRatio / scaleFactor - 1.) > 0.05 ) 
    { 
      cerr << " WARNING : header pixel sizes ratio=" << pixRatio<< " doesn't match scale factor="
	   << scaleFactor << " \n for images " << ToAlign.Name() << " and " << Ref.Name() << endl;
    }
}



ImageGtransfo::ImageGtransfo()
{
  //transfoFromRef = NULL;
  //transfoToRef = NULL;
  scaleFactor = 1; //? 
}


bool ImageGtransfo::IsValid() const
{
  return (transfoFromRef != NULL);
}

void ImageGtransfo::dump(ostream &s) const
{
  s << "*** ImageGtransfo ***" << endl; 
  s << " Geom Ref Name : " << GeomRefName() << endl;
  s << " transfoFromRef : " << endl << *TransfoFromRef() 
    << " transfoRoRef : " << endl << *TransfoToRef();
}

/* DefaultVal could be found in Source. but this routine is also used for boolean images where 
   DefaultVal has to be 0. */

#include <time.h> // for CPU printouts
void ImageGtransfo::TransformImage(const FitsImage &Original, FitsImage& Transformed, 
				   const ReducedImage *Source, ReducedImage *Result, 
				   double DefaultVal) const
{
#ifdef DEBUG
  cout  << " > ImageGtransfo::TransformImage(...) " << endl;
#endif

  /* register the image */
  clock_t tstart = clock();
  Image &transformedImage = Transformed;
  int nx_aligned = int(outputImageSize.Nx());
  int ny_aligned = int(outputImageSize.Ny());
  
  transformedImage = GtransfoImage(Original,*transfoFromRef, int(nx_aligned),
				   int(ny_aligned), DefaultVal, 3);

  transformedImage *= (1/(scaleFactor*scaleFactor)); // account for the average Jacobian.
  
  /* write in the header the frame coordinates */
  Frame frame(dynamic_cast<const FitsHeader&> (Original));
  cout << " Before scaling " << frame << endl;
  //frame = frame.ApplyTransfo(*transfoToRef);
  //  GtransfoLin *lintransfoToRef = GtransfoToLin(transfoToRef);
  //  frame = frame.ApplyTransfo(*lintransfoToRef);

  /* it's better to do it this way due to eventual wrong third order transfo terms*/
  GtransfoLin lintransfoToRef = transfoFromRef->
    LinearApproximation(Frame(transformedImage).Center(), 
			min(nx_aligned,ny_aligned)).invert();
  frame = ApplyTransfo(frame,lintransfoToRef);
  
  // watch it does not go outside image  
  frame *= Frame(dynamic_cast<const Image&> (Transformed));
  cout << " After scaling  " << frame << endl;
  
  Transformed.AddOrModKey("SCALFACT",scaleFactor,"Scaling factor when transforming");

  
  /* The update of the WCS consists in copying the WCS of the geomRef
     in the Transformed image, provided the first one is accurate */
  {
    bool writtenOutWCS = false;
    ReducedImage geomRef(geomRefName);
    if (geomRef.IsValid() && geomRef.HasImage())
      {
	FitsHeader geomHead(geomRef.FitsName());
	writtenOutWCS = CopyWCS(geomHead, Transformed);
      }
    else
      {
	Gtransfo* wcsSource;
	if (WCSFromHeader(Original, wcsSource))
	  {
	    Gtransfo* outWCS = GtransfoCompose(wcsSource, transfoFromRef);
	    TanPix2RaDec *tanWCS = dynamic_cast<TanPix2RaDec *>(outWCS);
	    if (tanWCS)
	      {
		TanWCS2Header(Transformed, *tanWCS);
		writtenOutWCS = true;
	      }
	    if (outWCS) delete outWCS;
	  }
	if (wcsSource) delete wcsSource;
      }
    if (!writtenOutWCS)
      {
	cout << " NO WCS in output image for " 
	     << Transformed.FileName() << std::endl;
      }
  }
  
  // update usable part
  frame.WriteInHeader(Transformed);

  // write a correct DATASEC in header if any
  if (Transformed.HasKey("DATASEC"))
    {
      char datasec[80];
      sprintf(datasec,"[1:%d,1:%d]", nx_aligned,ny_aligned);
      Transformed.AddOrModKey("DATASEC",string(datasec)," updated by ImageGtransfo");
    }

  // update RA and DEC
  //  if we copy the WCS we should not update RA and Dec
  //  string ra = Original.KeyVal("TOADRASC");
  //  string dec = Original.KeyVal("TOADDECL");
  //  double epoch = Original.KeyVal("TOADEQUI");
  //  Transformed.AddOrModKey("TOADRASC",ra.c_str(),"Modified RA after transformation");
  // Transformed.AddOrModKey("TOADDECL",dec.c_str(),"Modified DEC after transformation");
  // Transformed.AddOrModKey("TOADEQUI",epoch,"Modified epoch after transformation");
  
  if (Source && Result) 
    {
      // update pixel size *scale
      double pixsize = Source->PixelSize();
      Result->SetPixelSize(pixsize/scaleFactor,"Pixel size scaled after transformation");
      double scale2 = scaleFactor*scaleFactor;
      // update background level /scale^2
      double backlev = Source->BackLevel();
      Result->SetBackLevel(backlev/scale2,"Background level after transformation");
      // could add at least interpolation error at this stage  
      Result->SetFlatFieldNoise(Source->FlatFieldNoise()/scaleFactor,"Flatfield noise sigma after transformation");
      Result->SetReadoutNoise(Source->ReadoutNoise()/scaleFactor, "Readout noise after transformation");
      Result->SetSigmaBack(Source->SigmaBack()/scaleFactor, "Background fluctuations sigma after transformation");
      // update sky value (!= background if sky subtracted)
      double skylev = Source->OriginalSkyLevel(); 
      Result->SetOriginalSkyLevel(skylev/scale2,"Sky level after transformation");
	
      // update saturation level /scale^2
      double origsatur = Source->OriginalSaturation();  
      Result->SetOriginalSaturation(origsatur/scale2," original saturation after transformation");
      double satur = Source->Saturation();  
      Result->SetSaturation(satur/scale2,"Saturation after transformation");
      
      // update seeing *scale 
      double seeing = Source->Seeing();
      // assume rebinning is equivalent to convolving with a gaussian of 
      // sigma = (pixel size)/sqrt(12) = stdev of step function
      seeing = sqrt(seeing*seeing*scale2 + scale2/12);
      cout << " Expected seeing " << seeing << endl;
      Result->SetSeeing(seeing,"seeing in pixel sigma after transformation");
    }
  clock_t tend = clock();
  cout << " CPU for resampling " 
     << float(tend-tstart)/float(CLOCKS_PER_SEC) << endl;
}

void ImageGtransfo::TransformWeightImage(const FitsImage &Original, 
					 FitsImage& Transformed) const
{
  /* register the image */
  clock_t tstart = clock();
  Image inputWeights = Original;
  {
    Pixel min,max;
    inputWeights.MinMaxValue(&min,&max);
    cout << Original.FileName() << " min, max " << min << ' ' << max << endl;
  }
  double eps = 1e-15;
  inputWeights += eps; // so that it can be inverted
  inputWeights = 1/inputWeights; // inputWeights is now a variance map
  {
    Pixel min,max;
    inputWeights.MinMaxValue(&min,&max);
    cout << "after invertion" << " min, max " << min << ' ' << max << endl;
  }

  int nx_aligned = int(outputImageSize.Nx());
  int ny_aligned = int(outputImageSize.Ny());
  Image &transformedVariance = Transformed;
  transformedVariance = 
    GtransfoImage(inputWeights,*transfoFromRef, int(nx_aligned),
			       int(ny_aligned), 0, 3, true);
  transformedVariance *= 1./(scaleFactor*scaleFactor);
  
  {
    Pixel min,max;
    transformedVariance.MinMaxValue(&min,&max);
    cout << "after transfo" << " min, max " << min << ' ' << max << endl;
  }

  // transform back to weights and put low pixels (value comparable to eps) to 0;
  Pixel bad_pix = 10.*eps;
  Pixel *pend = transformedVariance.end();
  for (Pixel *p = transformedVariance.begin(); p < pend; ++p)
    {
      if (*p <= 0) { *p = 0; continue;} 
      *p = 1./(*p); // come back to weights
      if ((*p) < bad_pix) *p = 0;
    }

  // very important : preserves that 0 remain 0 after FITS write/read
  Transformed.PreserveZeros(); 

  {
    Pixel min,max;
    Transformed.MinMaxValue(&min,&max);
    cout << Transformed.FileName() << " min, max " << min << ' ' << max << endl;
  }
  

  /* write in the header the frame coordinates */
  Frame frame(dynamic_cast<const FitsHeader&> (Original));
  cout << " Before scaling " << frame << endl;
  //frame = frame.ApplyTransfo(*transfoToRef);
  //  GtransfoLin *lintransfoToRef = GtransfoToLin(transfoToRef);
  //  frame = frame.ApplyTransfo(*lintransfoToRef);

  /* it's better to do it this way due to eventual wrong third order transfo terms*/
  GtransfoLin lintransfoToRef = transfoFromRef->
    LinearApproximation(Point(nx_aligned/2, ny_aligned/2), 
			min(nx_aligned,ny_aligned)).invert();
  frame = ApplyTransfo(frame,lintransfoToRef);

  // watch it does not go outside image  
  frame *= Frame(dynamic_cast<const Image&> (Transformed));
  cout << " After scaling  " << frame;

  Transformed.AddOrModKey("SCALFACT",scaleFactor,"Scaling factor when transforming");

  // Whould in fact update the WCS with a linear approximation 
  // if (HasWCS(Transformed)) RemoveWCS(Transformed)

  // update usable part
  frame.WriteInHeader(Transformed);

  /* The update of the WCS consists in copying the WCS of the geomRef
     in the Transformed image, provided the first one is accurate.
     Use the WCS of ref Image, because weights may not have a WCS ? */

  ReducedImage geomRef(geomRefName);
  if (geomRef.IsValid() && geomRef.HasImage())
    {
      FitsHeader geomHead(geomRef.FitsName());
      CopyWCS(geomHead, Transformed);
    }

  clock_t tend = clock();
  cout << " CPU for resampling weights " 
     << float(tend-tstart)/float(CLOCKS_PER_SEC) << endl;
}



void ImageGtransfo::TransformBoolImage(const FitsImage &Original, FitsImage& Transformed) const
{
  /* register the image */
  TransformImage(Original, Transformed, NULL, NULL, 0.);
  Transformed.Simplify(0.1);
  Transformed.ModKey("BITPIX",8);
}

/* this routine could be integrated into SEStar code.
   It does however 3 different things :
   - transform coordinates (and there are untransformed coordinates) .....
   - crop the list to the transformed frame
   - scale a few scores
*/
void ImageGtransfo::TransformCatalog(const SEStarList &Catalog, SEStarList &Transformed) const
{
  Transformed.clear();
  Catalog.CopyTo(Transformed);
  Transformed.ApplyTransfo(*transfoToRef);
  double scale2 = scaleFactor*scaleFactor;
  for (SEStarIterator it=Transformed.begin(); it!=Transformed.end(); )
    {
      SEStar *star = *it;
      if (!outputImageSize.InFrame(*star)) 
	{
	  it = Transformed.erase(it);
	  continue;
	}
      star->Fluxmax() /= scale2;
      star->Fond() /= scale2;
      star->Fwhm() *= scaleFactor;
      star->Kronradius() *= scaleFactor;
      star->Isoarea() *= scale2;
      star->Mxx() *= scale2;
      star->Myy() *= scale2;
      star->Mxy() *= scale2;
      star->A() *= scaleFactor;
      star->B() *= scaleFactor;
      ++it;
    }
}


void ImageGtransfo::TransformAperCatalog(const AperSEStarList &Catalog, AperSEStarList &Transformed) const
{
  // I think the active part of the routine should be in apersestar.cc
  /* should transform also the global values ...... */
  // can't do this, buggy
  // TransformCatalog((const SEStarList &) Catalog, (SEStarList &) Transformed);
  Transformed.clear();
  Catalog.CopyTo(Transformed);
  Transformed.ApplyTransfo(*transfoToRef);
  double scale2 = scaleFactor*scaleFactor;
  for (AperSEStarIterator it=Transformed.begin(); it != Transformed.end(); )
    {
      AperSEStar *star = *it;
      if (!outputImageSize.InFrame(*star)) 
	{
	  it = Transformed.erase(it);
	  continue;
	}
      star->Fluxmax() /= scale2;
      star->Fond() /= scale2;
      star->Fwhm() *= scaleFactor;
      star->Kronradius() *= scaleFactor;
      star->Isoarea() *= scale2;
      star->Mxx() *= scale2;
      star->Myy() *= scale2;
      star->Mxy() *= scale2;
      star->A() *= scaleFactor;
      star->B() *= scaleFactor;
      star->neighborDist *= scaleFactor;
      star->gmxx *= scale2;
      star->gmyy *= scale2;
      star->gmxy *= scale2;
      for (vector<Aperture>::iterator sit=star->apers.begin(); sit != star->apers.end(); ++sit) {
	sit->radius *= scaleFactor;
	sit->ncos = int(ceil(sit->ncos*scale2));
      }
      ++it;
    }
  GlobalVal glob = Transformed.GlobVal();
  double seeing = glob.getDoubleValue("SEEING");
  double ndist  = glob.getDoubleValue("MAXNEIGHBOR_DIST");
  vector<double> rads = glob.getDoubleValues("RADIUS");
  vector<double> shapes = glob.getDoubleValues("STARSHAPE");
  glob.setDoubleValue("SEEING", seeing*scaleFactor);
  glob.setDoubleValue("MAXNEIGHBOR_DIST", ndist*scaleFactor);
  for (size_t i=0; i<rads.size(); ++i) rads[i] *= scaleFactor;
  for (size_t i=0; i<shapes.size(); ++i) shapes[i] *= scaleFactor;
  glob.setDoubleValues("RADIUS", rads);
  glob.setDoubleValues("STARSHAPE", shapes);
}




ImageGtransfo::~ImageGtransfo()
{

}



/********************* TransformedImage ************************************************/
void TransformedImage::init(const ReducedImage &Source,
			    const ImageTransfo *Transfo) 
{
  source = Source.Clone();
  sourceName = Source.Name();
  transfo = Transfo->Clone();
}

TransformedImage::TransformedImage(const string &TransformedName, 
				   const ReducedImage &Source,
				   const ImageTransfo *Transfo) 
       : ReducedImage(TransformedName)
{
#ifdef DEBUG
  cout << " > TransformedImage::TransformedImage(...) TransformedName = " 
       << TransformedName
       << endl;
#endif

  init(Source,Transfo);
  // should we call Create by default here? answer : yes as a trial
  Create("here");
}




TransformedImage::TransformedImage(const string &Name) : ReducedImage(Name)
{
}



void TransformedImage::dump(ostream& s) const
{
  s << "**** TransformedImage : " << Name() << endl;
  s << " source name : " <<  SourceName() << endl;
  transfo->dump(s);
}


const  ImageGtransfoRef TransformedImage::IMAGEGTransfo() const
{ 
  // is the stored transfo an  ImageGtransfo ?
  const ImageGtransfo *test = transfo; 
  if (!test) { cerr << "Dynamic cast failed in TransformedImage::GTransfo " << endl ; return ImageGtransfoRef();} // ie NULL
  ImageGtransfoRef gref(test);
  return(gref);
}


const Gtransfo* TransformedImage::FromRef() const
{
  const ImageGtransfoRef gref = IMAGEGTransfo();
  if (gref == NULL ) return NULL;
  return gref->FromRef();
}

#ifdef STORAGE
//ReducedImageRef TransformedImage::GeometricReference()
//{
//  const ImageGtransfoRef gref = IMAGEGTransfo();
//  if (gref == NULL ) return CountedRef<ReducedImage>(); // soit = NULL
//  geomRef = ReducedImageRead(gref->GeomRefName());
/  return geomRef;
//}

//string TransformedImage::GeomRefName() const
//{
//   const ImageGtransfoRef gref = IMAGEGTransfo();
//  if (gref == NULL) return NULL;
//  return (gref->GeomRefName());
//}
#endif

bool  TransformedImage::MakeCatalog() 
{
  string fileName = ReducedImage::CatalogName();
  if (FileExists(fileName)) return true;

  SEStarList inList(source->CatalogName());
  SEStarList outList;
  cout << " Transforming list from "<< source->Name() << endl;
  transfo->TransformCatalog(inList,outList);
  outList.write(fileName);
  return true;
}


#include <apersestar.h>
bool  TransformedImage::MakeAperCat() 
{
  string fileName = ReducedImage::AperCatalogName();
  if (FileExists(fileName)) return true;

  AperSEStarList inList(source->AperCatalogName());
  AperSEStarList outList;
  cout << " Transforming list from "<< source->Name() << endl;
  transfo->TransformAperCatalog(inList,outList);
  outList.write(fileName);
  return true;
}




bool TransformedImage::MakeFits() 
{
#ifdef DEBUG
  cout << " > TransformedImage::MakeFits() " << endl;
#endif
  string fileName = FitsName();
  if (FileExists(fileName)) return true;
  cout << " Transforming image "<< source->Name() << endl;
  FitsImage inFits(source->FitsName());
  FitsImage outFits(fileName, (FitsHeader &) inFits);
  cout << "============================= " << endl;
  double defaultVal = source->BackLevel();
  cout << " Default value for outside frame " << defaultVal << endl;
  transfo->TransformImage(inFits, outFits, source, this, defaultVal);
  cout << " TransformedImage is produced. Updating parameters. " << endl;  
  return true;
}

bool TransformedImage::MakeWeight() 
{
  string fileName = FitsWeightName();
  if (FileExists(fileName)) return true;
  cout << " Transforming weight image "<< source->Name() << endl;
  if (!source->HasWeight() && ! source->MakeWeight())
    {
      cerr <<" cannot make weight for " << Name() << endl;
      return false;
    }
  FitsImage inFits(source->FitsWeightName());
  FitsImage outFits(fileName, (FitsHeader&) inFits);
  transfo->TransformWeightImage(inFits, outFits);
  {
    FitsImage im(FitsName());
    cout << Name() << " Image/Weight stat :" << ImageAndWeightError(im, outFits) << endl;
  }
  return true;
}



// not easy to make a routine that does the right thing  for the 3 Make{Fits,Satur,Dead}
#include "allreducedimage.h"

bool TransformedImage::MakeDead()
{
  string outFileName = FitsDeadName();
  if (FileExists(outFileName)) return true;
  cout << " Transforming dead image "<< source->Name() << endl;
  string inName = source->FitsDeadName();
  if (!FileExists(inName) && !source->MakeDead()) return false;
  FitsImage inFits(inName);
  FitsImage outFits(outFileName,(FitsHeader&) inFits);
// TODO : pass 'this' as argument to modify e.g. seeing
  transfo->TransformBoolImage(inFits, outFits); 
  return true;
}

bool TransformedImage::MakeCosmic()
{
  string outFileName = FitsCosmicName();
  if (FileExists(outFileName)) return true;
  cout << " Transforming cosmic image "<< source->Name() << endl;
  string inName = source->FitsCosmicName();
  if (!FileExists(inName) && !source->MakeCosmic()) return false;
  FitsImage inFits(inName);
  FitsImage outFits(outFileName,  (FitsHeader &) inFits);
// TODO : pass 'this' as argument to modify e.g. seeing
  transfo->TransformBoolImage(inFits, outFits); 
  return true;
}
bool TransformedImage::MakeSatellite()
{
  string outFileName = FitsSatelliteName();
  if (FileExists(outFileName)) return true;
  cout << " Transforming satellite image "<< source->Name() << endl;
  string inName = source->FitsSatelliteName();
  if (!FileExists(inName) && !source->MakeSatellite()) return false;
  FitsImage inFits(inName);
  FitsImage outFits(outFileName,  (FitsHeader &) inFits);
// TODO : pass 'this' as argument to modify e.g. seeing
  transfo->TransformBoolImage(inFits, outFits); 
  return true;
}

bool TransformedImage::MakeSatur()
{
  string outFileName = FitsSaturName();
  if (FileExists(outFileName)) return true;
  cout << " Transforming satur image "<< source->Name() << endl;
  string inName = source->FitsSaturName();
  if (!FileExists(inName) && !source->MakeSatur()) return false;
  FitsImage inFits(inName);
  FitsImage outFits(outFileName, (FitsHeader&) inFits);
// TODO : pass 'this' as argument to modify e.g. seeing
  transfo->TransformBoolImage(inFits, outFits); 
  return true;
}


TransformedImage::TransformedImage(const TransformedImage &Original) :
  ReducedImage(Original.Name())
{
  this->init(*Original.source, Original.transfo);
}


ReducedImage *TransformedImage::Clone() const
{
  return new TransformedImage(Name(), *source, transfo);
}


TransformedImage::~TransformedImage()
{
}


/******************* utilities ***********************/

void MakeUnionRef(const ReducedImageList& ToAlign, const ReducedImage& Reference, const string& unionName) {
  Frame unionFrame(Reference.UsablePart());
  // enlarge frame to secure resampling extra width
  double hwidth = 5;
  unionFrame.CutMargin(-hwidth);
  cout << " MakeUnionRef: initial frame is: " << unionFrame << endl;
  for (ReducedImageCIterator it = ToAlign.begin(); it != ToAlign.end(); ++it) {
    CountedRef<Gtransfo> reftoim, imtoref;
    ImageListMatch(Reference, **it, reftoim, imtoref);
    Frame frameInRef = ApplyTransfo((*it)->UsablePart(), *imtoref, LargeFrame);
    frameInRef.CutMargin(-hwidth);
    unionFrame += frameInRef;
  } 

  // make an integer frame to avoid resampling
  unionFrame.xMin = floor(unionFrame.xMin);
  unionFrame.yMin = floor(unionFrame.yMin);
  unionFrame.xMax = ceil(unionFrame.xMax);
  unionFrame.yMax = ceil(unionFrame.yMax);

  cout << " MakeUnionRef: final frame is " << unionFrame << endl;

  double dx = unionFrame.xMin;
  double dy = unionFrame.yMin;
  CountedRef<Gtransfo> TransfoFromRef = new GtransfoLinShift(dx, dy);
  CountedRef<Gtransfo> TransfoToRef   = new GtransfoLinShift(-dx, -dy);

  ImageGtransfo transfo(TransfoFromRef, TransfoToRef, unionFrame, Reference.Name());
  
  // now transform image
  TransformedImage transformed(unionName, Reference, &transfo);
  transformed.Execute(DoFits | DoCatalog | DoSatur | DoWeight | DoCosmic);

  // apercat are not handled by the DoCrap 
  string filename = Reference.AperCatalogName();
  if (FileExists(filename)) {
    AperSEStarList inList(filename);
    AperSEStarList outList;
    cout << " Transforming list from " << filename << endl;
    transfo.TransformAperCatalog(inList, outList);
    outList.write(transformed.AperCatalogName());
  }
  filename = Reference.StarCatalogName();
  if (FileExists(filename)) {
    AperSEStarList inList(filename);
    AperSEStarList outList;
    cout << " Transforming list from " << filename << endl;
    transfo.TransformAperCatalog(inList, outList);
    outList.write(transformed.StarCatalogName());
  }
}


int ImagesAlign(const ReducedImageList &ToAlign, const ReducedImage &Reference, 
		ReducedImageList &Aligned, const int ToDo, bool use_wcs, float min_match_ratio)
{
#ifdef DEBUG
  cout << " > ImagesAlign(...) " << endl;
#endif

  Aligned.clear();
  for (ReducedImageCIterator ri = ToAlign.begin(); ri != ToAlign.end(); ++ri)
    {
      const ReducedImage *current = *ri;
      string currentName = current->Name();
      if (!current->IsValid())
	{
	  cerr << "ImagesAlign: was " << currentName << " actually produced ?? " << endl;
	  continue;
        }
      if (*current == Reference) {
	Aligned.push_back(current);
	continue;
      }

      string transformedName = TransformedName(currentName, Reference.Name());
      ImageGtransfo *imTransfo = 0;

      if (use_wcs) {
	Gtransfo* wcs_reference;
	WCSFromHeader(Reference.FitsName(),wcs_reference);
	Gtransfo* wcs_current;
	WCSFromHeader(current->FitsName(),wcs_current);
	Frame frame_reference = Reference.UsablePart();
	Frame frame_current = current->UsablePart();
	
	Gtransfo *direct = GtransfoCompose(wcs_reference->InverseTransfo(0.01,frame_reference),wcs_current);
	Gtransfo *invert = GtransfoCompose(wcs_current->InverseTransfo(0.01,frame_current),wcs_reference);
	imTransfo = new ImageGtransfo(invert,direct,frame_reference, Reference.Name());
	delete wcs_reference;
	delete wcs_current;
	delete direct;
	delete invert;
      }else{
	imTransfo = new ImageGtransfo(Reference,*current,min_match_ratio);
      }
      
      TransformedImage *alignedCurrent =
	new TransformedImage(transformedName, *current, imTransfo);
      // create the Fits images + user requests, as doc says.
      alignedCurrent->Execute(DoFits | ToDo); 
      Aligned.push_back(alignedCurrent);
      delete imTransfo;
    }
  return Aligned.size();
}



/********************** TransformedImageList *************************/
//template class ImageList<TransformedImage>;


#ifdef USE_ROOT

ClassImp(TransformedImage);

/*
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ class ImageGtransfo+;
LINKDEF_CONTENT : #pragma link C++ class ImageTransfo+;
LINKDEF_CONTENT : #pragma link C++ class TransformedImage+;
*/

#include "root_dict/transformedimagedict.cc"

#endif /* USE_ROOT */
