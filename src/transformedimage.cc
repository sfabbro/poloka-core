#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "transformedimage.h"
#include "gtransfo.h"
#include "fitsimage.h"
#include "wcsutils.h"
#include "fileutils.h"
#include "imageutils.h"
#include "reducedutils.h"
#include "allreducedimage.h"
#include "polokaexception.h"
#include "apersestar.h"

string TransformedName(const string &ToTransform, const string &Ref)
{
  return "T_" + Ref + "_" + ToTransform;
}


ImageGtransfo::ImageGtransfo(const GtransfoRef TransfoFromRef, 
			     const GtransfoRef TransfoToRef, 
			     const Frame &OutputImageSize, 
			     const string &GeomRefName)
{
  geomRefName = GeomRefName;
  transfoFromRef = TransfoFromRef;
  transfoToRef = TransfoToRef;
  outputImageSize = OutputImageSize; 
}

double ImageGtransfo::ScaleFactor() const 
{
  if (transfoFromRef)
    return 1./sqrt(fabs(transfoFromRef->Jacobian(outputImageSize.Nx()/2, outputImageSize.Ny()/2)));
  return 1.;
}

const GtransfoRef ImageGtransfo::FromRef() const
{
  return transfoFromRef;
}


ImageTransfo* ImageGtransfo::Clone() const
{
  return new ImageGtransfo(transfoFromRef, transfoToRef, outputImageSize, geomRefName);
}


ImageGtransfo::ImageGtransfo(const ReducedImage &Ref, const ReducedImage& ToAlign)
{
  geomRefName = Ref.Name();
  transfoFromRef = FindTransfo(Ref, ToAlign);
  if (!transfoFromRef)
    throw(PolokaException("ImageGtransfo: could not match lists " + Ref.Name() + " and " + ToAlign.Name()));

  transfoToRef = FindTransfo(ToAlign, Ref);
  if (!transfoToRef)
    throw(PolokaException("ImageGtransfo: could not invert list " + Ref.Name() + " and " + ToAlign.Name()));

  outputImageSize = Ref.PhysicalSize();
  // use the same transfo as the one used to actually transform the image.
  double scale = ScaleFactor();
  cout << " ImageGtransfo: scale factor for rebinning = " << scale << endl;
  double pixRatio = ToAlign.PixelSize()/Ref.PixelSize();
  if ( fabs( pixRatio / scale - 1.) > 0.05 ) 
    { 
      cerr << " ImageGtransfo: WARNING header pixel sizes ratio=" << pixRatio<< " doesn't match scale factor="
	   << scale << " \n for images " << ToAlign.Name() << " and " << Ref.Name() << endl;
    }
}

bool ImageGtransfo::IsValid() const
{
  return (transfoFromRef != NULL);
}

void ImageGtransfo::dump(ostream &s) const
{
  s << "*** ImageGtransfo ***" << endl
    << " Geom Ref Name : " << GeomRefName() << endl
    << " transfoFromRef : " << endl << *TransfoFromRef() 
    << " transfoToRef : " << endl << *TransfoToRef();
}

/* DefaultVal could be found in Source. but this routine is also used for boolean images where 
   DefaultVal has to be 0. */
void ImageGtransfo::TransformImage(const FitsImage &Original, FitsImage& Transformed, 
				   const ReducedImage *Source, ReducedImage *Result, 
				   double DefaultVal) const
{

  // resample the image
  Image &transformedImage = Transformed;
  int nx_aligned = int(outputImageSize.Nx());
  int ny_aligned = int(outputImageSize.Ny());
  
  transformedImage = GtransfoImage(Original,*transfoFromRef, int(nx_aligned),
				   int(ny_aligned), DefaultVal, 3);

   // account for the average Jacobian.
  double scale = ScaleFactor();
  double scale2 = scale*scale;
  transformedImage *= 1./scale2;
  
  // write in the header the frame coordinates
  Frame frame(dynamic_cast<const FitsHeader&> (Original));
  cout << " ImageGtransfo: frame before scaling " << frame;

  // it's better to approximate a linear transfo due to eventual wrong third order transfo terms
  GtransfoLin lintransfoToRef = transfoFromRef->
    LinearApproximation(Frame(transformedImage).Center(),
			min(nx_aligned,ny_aligned)).invert();
  frame = ApplyTransfo(frame,lintransfoToRef);
  
  // watch it does not go outside image  
  frame *= Frame(dynamic_cast<const Image&> (Transformed));
  cout << " ImageGtransfo: frame after scaling  " << frame;
  
  // write a few things in the resampled header
  Transformed.AddOrModKey("SCALFACT",scale,"Scaling factor when transforming");
  Transformed.AddOrModKey("GEOREF", geomRefName, "DbImage geometric reference for alignment");
  
  // update of the WCS: consists of copying the WCS of the geomRef
  // in the Transformed image, provided the first one is accurate
  {
    bool writtenOutWCS = false;
    ReducedImage geomRef(geomRefName);
    if (geomRef.IsValid() && geomRef.HasImage() && geomRef.XSize() == nx_aligned && geomRef.YSize() == ny_aligned)
      {
	FitsHeader geomHead(geomRef.FitsName());
	writtenOutWCS = CopyWCS(geomHead, Transformed);
      }
    else
      {
	GtransfoRef wcsSource = WCSFromHeader(Original);
	if (wcsSource)
	  {
	    GtransfoRef outWCS = GtransfoCompose(wcsSource, transfoFromRef);
	    TanPix2RaDec *tanWCS = dynamic_cast<TanPix2RaDec *>((Gtransfo*)outWCS);
	    if (tanWCS)
	      {
		TanWCS2Header(Transformed, *tanWCS);
		writtenOutWCS = true;
	      }
	  }
      }
    if (!writtenOutWCS)
      {
	cout << " ImageGtransfo: no WCS in output image for "
	     << Transformed.FileName() << endl;
      }
  }
  
  // update usable frame
  frame.WriteInHeader(Transformed);
  
  // update other header values
  if (Source && Result) 
    {
      Transformed.AddOrModKey("DBSOURCE", Source->Name(), "Source DbImage before resampling");
      // pixel size * scale
      double pixsize = Source->PixelSize();
      Result->SetPixelSize(pixsize/scale,"Pixel size scaled after transformation");

      // background level /scale^2
      double backlev = Source->BackLevel();
      Result->SetBackLevel(backlev/scale2,"Background level after transformation");
      // could add at least interpolation error at this stage  
      Result->SetFlatFieldNoise(Source->FlatFieldNoise()/scale,"Flatfield noise sigma after transformation");
      Result->SetReadoutNoise(Source->ReadoutNoise()/scale, "Readout noise after transformation");
      Result->SetSigmaBack(Source->SigmaBack()/scale, "Background fluctuations sigma after transformation");
      // sky value (!= background if sky subtracted)
      double skylev = Source->OriginalSkyLevel(); 
      Result->SetOriginalSkyLevel(skylev/scale2,"Sky level after transformation");	
      // saturation level /scale^2
      double origsatur = Source->OriginalSaturation();  
      Result->SetOriginalSaturation(origsatur/scale2," original saturation after transformation");
      double satur = Source->Saturation();  
      Result->SetSaturation(satur/scale2,"Saturation after transformation");
      // seeing *scale 
      double seeing = Source->Seeing();
      // assume rebinning is equivalent to convolving with a gaussian of 
      // sigma = (pixel size)/sqrt(12) = stdev of step function
      seeing = sqrt(seeing*seeing*scale2 + scale2/12);
      cout << " TransformedImage: expected seeing " << seeing << endl;
      Result->SetSeeing(seeing,"seeing in pixel sigma after transformation");
      Result->SetUsablePart(frame);
    }
}

static const Pixel BAD_WEIGHT = 1e-15;

void ImageGtransfo::TransformWeightImage(const FitsImage &Original, 
					 FitsImage& Transformed) const
{
  // transform weight into variance
  Image inputWeights = Original;
  Pixel minVal, maxVal;
  inputWeights.MinMaxValue(&minVal, &maxVal);
  cout << " ImageGtransfo: original weight   " << Original.FileName() << " min, max " << minVal << ' ' << maxVal << endl;
  Pixel *pend = inputWeights.end();
  for (Pixel *p = inputWeights.begin(); p < pend; ++p)
    *p = 1./(*p + BAD_WEIGHT);

  // resample the variance
  int nx_aligned = int(outputImageSize.Nx());
  int ny_aligned = int(outputImageSize.Ny());
  Image &transformedVariance = Transformed;
  transformedVariance = 
    GtransfoImage(inputWeights,*transfoFromRef, int(nx_aligned),
			       int(ny_aligned), 0, 3, true);
  
  double scale = ScaleFactor();

  // account for Jacobian, transform back to weights and put low pixels (value comparable to BAD_WEIGHT) to 0
  double scale2 = scale*scale;
  Pixel maxvar = 1. / (BAD_WEIGHT*10.);
  pend = transformedVariance.end();
  for (Pixel *p = transformedVariance.begin(); p < pend; ++p)
    (*p > maxvar || *p <= 0) ? *p = 0. : *p = scale2 / *p;

  // very important : preserves that 0 remain 0 after FITS write/read
  Transformed.PreserveZeros(); 
  Transformed.MinMaxValue(&minVal, &maxVal);
  cout << " ImageGtransfo: transformed weight " << Original.FileName() << " min, max " << minVal << ' ' << maxVal << endl;

  // write in the header the frame coordinates
  Frame frame(dynamic_cast<const FitsHeader&> (Original));
  cout << " ImageGtransfo: frame before scaling " << frame;

  // it's better to approximate a linear transfo due to eventual wrong third order transfo terms
  GtransfoLin lintransfoToRef = transfoFromRef->
    LinearApproximation(Point(nx_aligned/2, ny_aligned/2), 
			min(nx_aligned,ny_aligned)).invert();
  frame = ApplyTransfo(frame,lintransfoToRef);
  
  // watch it does not go outside image  
  frame *= Frame(dynamic_cast<const Image&> (Transformed));
  cout << " ImageGtransfo: frame after scaling  " << frame;
  
  // write a few things in the resampled header
  Transformed.AddOrModKey("SCALFACT", scale, "Scaling factor when transforming");

  
  // update of the WCS: consists of copying the WCS of the geomRef
  // in the Transformed image, provided the first one is accurate
  {
    bool writtenOutWCS = false;
    ReducedImage geomRef(geomRefName);
    if (geomRef.IsValid() && geomRef.HasImage() && geomRef.XSize() == nx_aligned && geomRef.YSize() == ny_aligned)
      {
	FitsHeader geomHead(geomRef.FitsName());
	writtenOutWCS = CopyWCS(geomHead, Transformed);
      }
    else
      {
	GtransfoRef wcsSource = WCSFromHeader(Original);
	if (wcsSource)
	  {
	    GtransfoRef outWCS = GtransfoCompose(wcsSource, transfoFromRef);
	    TanPix2RaDec *tanWCS = dynamic_cast<TanPix2RaDec *>((Gtransfo*)outWCS);
	    if (tanWCS)
	      {
		TanWCS2Header(Transformed, *tanWCS);
		writtenOutWCS = true;
	      }
	  }
      }
    if (!writtenOutWCS)
      {
	cout << " ImageGtransfo: no WCS in output image for "
	     << Transformed.FileName() << endl;
      }
  }
  
  // update usable frame
  frame.WriteInHeader(Transformed);  
}

void ImageGtransfo::TransformBoolImage(const FitsImage &Original, FitsImage& Transformed) const
{
  // use the routine above to resample and update header
  TransformImage(Original, Transformed, NULL, NULL, 0.);
  // keep it binary
  Transformed.Simplify(0.1);
  Transformed.ModKey("BITPIX",8);
}

/* this routine could be integrated into SEStar code.
   It does however 3 different things :
   - transform coordinates (and there are untransformed coordinates) .....
   - crop the list to the transformed frame
   - scale a few scores */
void ImageGtransfo::TransformCatalog(const SEStarList &Catalog, SEStarList &Transformed) const
{
  Transformed.clear();
  Catalog.CopyTo(Transformed);
  Transformed.ApplyTransfo(*transfoToRef);
  double scale = ScaleFactor();
  double scale2 = scale*scale;
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
      star->Fwhm() *= scale;
      star->Kronradius() *= scale;
      star->Isoarea() *= scale2;
      star->Mxx() *= scale2;
      star->Myy() *= scale2;
      star->Mxy() *= scale2;
      star->A() *= scale;
      star->B() *= scale;
      ++it;
    }
}

void ImageGtransfo::TransformAperCatalog(const AperSEStarList &Catalog, AperSEStarList &Transformed) const
{
  // I think the active part of the routine should be in apersestar.cc
  // TransformCatalog((const SEStarList &) Catalog, (SEStarList &) Transformed); // buggy on runtime
  Transformed.clear();
  Catalog.CopyTo(Transformed);
  Transformed.ApplyTransfo(*transfoToRef);
  double scale = ScaleFactor();
  double scale2 = scale*scale;
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
      star->Fwhm() *= scale;
      star->Kronradius() *= scale;
      star->Isoarea() *= scale2;
      star->Mxx() *= scale2;
      star->Myy() *= scale2;
      star->Mxy() *= scale2;
      star->A() *= scale;
      star->B() *= scale;
      star->neighborDist *= scale;
      star->gmxx *= scale2;
      star->gmyy *= scale2;
      star->gmxy *= scale2;
      for (vector<Aperture>::iterator sit=star->apers.begin(); sit != star->apers.end(); ++sit) {
	sit->radius *= scale;
	sit->ncos = int(ceil(sit->ncos*scale2));
      }
      ++it;
    }
  GlobalVal glob = Transformed.GlobVal();
  double seeing = glob.getDoubleValue("SEEING");
  double ndist  = glob.getDoubleValue("MAXNEIGHBOR_DIST");
  vector<double> rads = glob.getDoubleValues("RADIUS");
  vector<double> shapes = glob.getDoubleValues("STARSHAPE");
  glob.setDoubleValue("SEEING", seeing*scale);
  glob.setDoubleValue("MAXNEIGHBOR_DIST", ndist*scale);
  for (size_t i=0; i<rads.size(); ++i) rads[i] *= scale;
  for (size_t i=0; i<shapes.size(); ++i) shapes[i] *= scale;
  glob.setDoubleValues("RADIUS", rads);
  glob.setDoubleValues("STARSHAPE", shapes);
}

void TransformedImage::init(const ReducedImage &Source,
			    const ImageTransfo *Transfo) 
{
  source = Source.Clone();
  transfo = Transfo->Clone();
}

TransformedImage::TransformedImage(const string &TransformedName, 
				   const ReducedImage &Source,
				   const ImageTransfo *Transfo) 
  : ReducedImage(TransformedName)
{
  init(Source,Transfo);
  // should we call Create by default here? answer : yes as a trial
  Create("here");
}

TransformedImage::TransformedImage(const TransformedImage &Original)
  : ReducedImage(Original.Name())
{
  this->init(*Original.source, Original.transfo);
}


TransformedImage::TransformedImage(const string &Name) : ReducedImage(Name)
{
  ReducedImage *ref = 0;
  ReducedImage *source = 0;
  if (HasImage())
    {
      FitsHeader head(FitsName());
      string refName = head.KeyVal("GEOREF");
      ref = new ReducedImage(refName);
      string sourceName = head.KeyVal("DBSOURCE");
      source = new ReducedImage(sourceName);
     }
  if (!ref || !source)
    {
      vector<string> dbNames;
      DecomposeString(dbNames, Name,"_");
       if (dbNames.size() == 3)
	{
	  ref = new ReducedImage(dbNames[1]);
	  source = new ReducedImage(dbNames[2]);
	}
      else if (dbNames.size() == 2 && ref)
	{
	  RemovePattern(dbNames[1],ref->Name());
	  source = new ReducedImage(dbNames[1]);
	}
      else
	{
	  cerr << " TransformedImage: " << Name << " is an old transformedimage, no persistence\n";
	}
    }
  // wont work with other ImageTransfo than ImageGtransfo (there's none yet anyway)
  if (ref && source) 
    transfo = new ImageGtransfo(*ref,*source);
  if (ref) delete ref;
  if (source) delete source;
}


ReducedImageRef TransformedImage::Clone() const
{
  return new TransformedImage(Name(), *source, transfo);
}

void TransformedImage::dump(ostream& s) const
{
  s << "**** TransformedImage : " << Name() << endl;
  if (source)
    s << " source name : " <<  source->Name() << endl;
  if (transfo)
    transfo->dump(s);
}


const ImageGtransfoRef TransformedImage::IMAGEGTransfo() const
{ 
  // is the stored transfo an  ImageGtransfo ?
  const ImageGtransfo *test = transfo; 
  if (!test)
    throw(PolokaException("TransformedImage: dynamic cast failed in IMAGEGtransfo"));
  ImageGtransfoRef gref(test);
  return(gref);
}


const GtransfoRef TransformedImage::FromRef() const
{
  const ImageGtransfoRef gref = IMAGEGTransfo();
  if (gref == NULL ) return GtransfoRef();
  return gref->FromRef();
}


bool TransformedImage::MakeCatalog() 
{
  string outFileName = CatalogName();
  if (FileExists(outFileName)) return true;

  if (!source || !transfo)
    throw(PolokaException("TransformedImage: need both source and transfo to produce " + outFileName));
  if (!source->HasCatalog() && !source->MakeCatalog())
    throw(PolokaException("TransformedImage: cannot make catalog for " + Name()));
 
  SEStarList inList(source->CatalogName());
  SEStarList outList;
  cout << " TransformedImage: transforming SExtractor catalog of "<< source->Name() << endl;
  transfo->TransformCatalog(inList,outList);
  outList.write(outFileName);
  return true;
}

bool TransformedImage::MakeAperCat() 
{
  string outFileName = AperCatalogName();
  if (FileExists(outFileName)) return true;

  if (!source || !transfo)
    throw(PolokaException("TransformedImage: need both source and transfo to produce " + outFileName));
  if (!source->HasAperCatalog() && !source->MakeAperCat())
    throw(PolokaException("TransformedImage: cannot make aper catalog for " + Name()));

  AperSEStarList inList(source->AperCatalogName());
  AperSEStarList outList;
  cout << " TransformedImage: transforming aperture catalog of "<< source->Name() << endl;
  transfo->TransformAperCatalog(inList, outList);
  outList.write(outFileName);
  return true;
}

bool TransformedImage::MakeStarCat() 
{
  string outFileName = StarCatalogName();
  if (FileExists(outFileName)) return true;

  if (!source || !transfo)
    throw(PolokaException("TransformedImage: need both source and transfo to produce " + outFileName));
  if (!source->HasStarCatalog() && !source->MakeStarCat())
    throw(PolokaException("TransformedImage: cannot make star catalog for " + Name()));

  AperSEStarList inList(source->StarCatalogName());
  AperSEStarList outList;
  cout << " TransformedImage: transforming star aperture catalog from "<< source->Name() << endl;
  transfo->TransformAperCatalog(inList, outList);
  outList.write(outFileName);
  return true;
}

bool TransformedImage::MakeFits() 
{
  string outFileName = FitsName();
  if (FileExists(outFileName)) return true;

  if (!source || !transfo)
    throw(PolokaException("TransformedImage: need both source and transfo to produce " + outFileName));
  if (!source->HasImage() && !source->MakeFits())
    throw(PolokaException("TransformedImage: cannot make fits for " + Name()));

  cout << " TransformedImage: transforming calibrated image of "<< source->Name() << endl;
  FitsImage inFits(source->FitsName());
  FitsImage outFits(outFileName, (FitsHeader &) inFits);
  double defaultVal = source->BackLevel();
  cout << " TransformedImage: default value for outside frame " << defaultVal << endl;
  transfo->TransformImage(inFits, outFits, source, this, defaultVal);
  return true;
}

bool TransformedImage::MakeWeight() 
{
  string outFileName = FitsWeightName();
  if (FileExists(outFileName)) return true;

  if (!source || !transfo)
    throw(PolokaException("TransformedImage: need both source and transfo to produce " + outFileName));
  if (!source->HasWeight() && !source->MakeWeight())
    throw(PolokaException("TransformedImage: cannot make weight for " + Name()));

  cout << " TransformedImage: transforming weight image of "<< source->Name() << endl;
  FitsImage inFits(source->FitsWeightName());
  FitsImage outFits(outFileName, (FitsHeader&) inFits);
  transfo->TransformWeightImage(inFits, outFits);

  return true;
}

// not easy to make a routine that does the right thing  for the 3 Make{Fits,Satur,Dead}
bool TransformedImage::MakeDead()
{
  string outFileName = FitsDeadName();
  if (FileExists(outFileName)) return true;

  if (!source || !transfo)
    throw(PolokaException("TransformedImage: need both source and transfo to produce " + outFileName));
  if (!source->HasDead() && !source->MakeDead())
    throw(PolokaException("TransformedImage: cannot make dead for " + Name()));

  cout << " TransformedImage: transforming dead image of "<< source->Name() << endl;
  string inName = source->FitsDeadName();
  FitsImage inFits(inName);
  FitsImage outFits(outFileName,(FitsHeader&) inFits);
  transfo->TransformBoolImage(inFits, outFits); 
  return true;
}

bool TransformedImage::MakeCosmic()
{
  string outFileName = FitsCosmicName();
  if (FileExists(outFileName)) return true;

  if (!source || !transfo)
    throw(PolokaException("TransformedImage: need both source and transfo to produce " + outFileName));
  if (!source->HasCosmic() && !source->MakeCosmic())
    throw(PolokaException("TransformedImage: cannot make cosmic for " + Name()));

  cout << " TransformedImage: transforming cosmic image of "<< source->Name() << endl;
  string inName = source->FitsCosmicName();
  FitsImage inFits(inName);
  FitsImage outFits(outFileName,  (FitsHeader &) inFits);
  transfo->TransformBoolImage(inFits, outFits); 
  return true;
}
bool TransformedImage::MakeSatellite()
{
  string outFileName = FitsSatelliteName();
  if (FileExists(outFileName)) return true;

  if (!source || !transfo)
    throw(PolokaException("TransformedImage: need both source and transfo to produce " + outFileName));
  if (!source->HasSatellite() && !source->MakeSatellite())
    throw(PolokaException("TransformedImage: cannot make satellite for " + Name()));


  cout << " TransformedImage: transforming satellite image of "<< source->Name() << endl;
  string inName = source->FitsSatelliteName();
  FitsImage inFits(inName);
  FitsImage outFits(outFileName,  (FitsHeader &) inFits);
  transfo->TransformBoolImage(inFits, outFits); 
  return true;
}

bool TransformedImage::MakeSatur()
{
  string outFileName = FitsSaturName();
  if (FileExists(outFileName)) return true;

  if (!source || !transfo)
    throw(PolokaException("TransformedImage: need both source and transfo to produce " + outFileName));
  if (!source->HasSatur() && !source->MakeSatur())
    throw(PolokaException("TransformedImage: cannot make satur for " + Name()));

  cout << " TransformImage: transforming satur image of "<< source->Name() << endl;
  string inName = source->FitsSaturName();
  FitsImage inFits(inName);
  FitsImage outFits(outFileName, (FitsHeader&) inFits);
  transfo->TransformBoolImage(inFits, outFits); 
  return true;
}


int ImagesAlign(const ReducedImageList &ToAlign, const ReducedImage &Reference, 
		ReducedImageList &Aligned, const int ToDo, bool WcsOnly)
{
  Aligned.clear();
  for (ReducedImageCIterator ri = ToAlign.begin(); ri != ToAlign.end(); ++ri)
    {
      const ReducedImage *current = *ri;
      string currentName = current->Name();
      if (!current->IsValid())
	{
	  cerr << "ImagesAlign: was " << currentName << " actually produced ?" << endl;
	  continue;
        }
      if (*current == Reference) {
	Aligned.push_back(current);
	continue;
      }

      string transformedName = TransformedName(currentName, Reference.Name());

      // skip transfo is the TransformedImage exists
      if (FileExists(DbImage(transformedName).FitsImageName(Calibrated))) {
	Aligned.push_back(new TransformedImage(transformedName));
	continue;
      }

      ImageGtransfo *imTransfo = 0;

      if (WcsOnly)
	{
	  GtransfoRef wcs_reference = WCSFromHeader(Reference.FitsName());
	  GtransfoRef wcs_current = WCSFromHeader(current->FitsName());
	  Frame frame_reference = Reference.UsablePart();
	  Frame frame_current = current->UsablePart();	  
	  GtransfoRef direct = GtransfoCompose(wcs_reference->InverseTransfo(0.01,frame_reference),wcs_current);
	  GtransfoRef invert = GtransfoCompose(wcs_current->InverseTransfo(0.01,frame_current),wcs_reference);
	  imTransfo = new ImageGtransfo(invert,direct,frame_reference, Reference.Name());
	}
      else
	{
	  imTransfo = new ImageGtransfo(Reference, *current);
	}
      
      TransformedImage *alignedCurrent = new TransformedImage(transformedName, *current, imTransfo);
      // create the Fits images + user requests, as doc says.
      alignedCurrent->Execute(DoFits | ToDo);
      Aligned.push_back(alignedCurrent);
      delete imTransfo;
    }
  return Aligned.size();
}
