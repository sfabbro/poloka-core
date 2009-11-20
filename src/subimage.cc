#include <stdio.h>

#include "fitsimage.h"
#include "subimage.h"
#include "frame.h"
#include "gtransfo.h"
#include "imageutils.h"


SubImage::SubImage(const string &Name, const string &LargeImageName, 
		   const Frame &SubFrame)
  : ReducedImage(Name)
{
  subFrame = SubFrame;
  largeImageName = LargeImageName;
  // get the actual size of the large image
  ReducedImage large(largeImageName);
  FitsHeader head(large.FitsName());
  largeFrame = Frame(head, WholeSizeFrame);
  // create the sub if needed
  Create("here");
}


static string frame_to_string( const Frame &AFrame, int &Imin, int &Imax)
{
  Imin = int(AFrame.xMin+0.5)+1;
  Imax = int(AFrame.xMax+0.5)+1;
  int jmin = int(AFrame.yMin+0.5)+1;
  int jmax = int(AFrame.yMax+0.5)+1;
  char cstring[256];
  sprintf(cstring,"[%d:%d,%d:%d]",Imin,Imax,jmin,jmax);
  return string(cstring);
}


#include "wcsutils.h"

static void ShiftWCSOrigin(const FitsHeader &In, FitsHeader &Out,
			   const int Dx, const int Dy)
{
  /* should not modify CRPIX{1,2} : cfitsio does it on the fly when
     adressing the image via fileName[xmin:xmax, ymin:ymax]. this is
     done in fits_select_image_section. So we do not put this routine
     (specific from the way we use cfitsio) in wcsutils.cc */
  /* To make a long story short, there is nothing to do about WCS's
     after cfitsio, except if distortion corrections apply in the
     pixel space : the only example we know about is the WCS3 private
     toads keywords */
  
#ifdef OLD_CODE
  GtransfoCub cubCorr;
  if (WCS3TransfoFromHeader(In, cubCorr, false))
    {
      CopyWCS(In, Out);
      GtransfoLinShift shift(Dx, Dy);
      cubCorr = cubCorr * shift;
      cubCorr = shift.invert() * cubCorr;
      WCS3Transfo2Header(Out, cubCorr);
    }
#endif
  // just check whether we are in trouble
  if (In.HasKey("WCS3DX"))
    {
      cerr << "Serious problem for WCS handling : we have an old WCS " << endl
	   << " in the header which is no longer handled properly " << endl;
      cerr << "For file =" <<  In.FileName() << endl;
      cerr <<" Aborting. Please contact developers " << endl;
      abort();
    }
      
	

}

bool SubImage::cut_in_fitsimage(string (ReducedImage::*GetFitsFileName)() const,
				bool (ReducedImage::*MakeInFits)())
{
  /* important notice :
     if you think it is a good idea to change the disk representation of 
     pixels here, refrain! because for weights, it is extremely important 
     that 0 remain 0. So do *NOT* alter BITPIX unless you very well know 
     what you are doing.
  */
  if (MakeInFits) {}; // to avoid a warning
  string toProduce = (this->*GetFitsFileName)();
  if (FileExists(toProduce)) return true;
  int dx=0, dy =0;
  string frame=frame_to_string(largeFrame*subFrame, dx, dy);
  ReducedImage large(largeImageName);
  
  // use the fact that cfitsio interprets correctly filename[i:j,k:l]
  string largeFileName = (large.*GetFitsFileName)();
  string cfitsio_subimage_name = largeFileName + frame;
  cout << " extracting " << toProduce << " as :" << endl
       << "    " << cfitsio_subimage_name << endl;
  FitsImage pseudoSubImage(cfitsio_subimage_name);
  if (!pseudoSubImage.IsValid()) return false;
  
  FitsHeader &head = pseudoSubImage;
  Image &pixels = pseudoSubImage;

  // create the output image
  FitsImage subImage(toProduce, head, pixels);
  // handle WCS
  ShiftWCSOrigin(head, subImage, dx, dy);
  // handle Frame describing the "usable part" frame
  if (HasClippedFrame(subImage))
    {
      Frame clippedFrame(subImage,ClippedSizeFrame);
      clippedFrame *= subFrame; // take intersection with
      clippedFrame = ApplyTransfo(clippedFrame,GtransfoLinShift(double(dx), double (dy)));
      clippedFrame.WriteInHeader(subImage);
    }
  // write the origin and subimage coordinates in large image:
  subImage.AddOrModKey("LARGEIMA",large.Name(), 
		       " DbImage from which this one was extracted");
  subImage.AddOrModKey("EXTRASEC", frame, 
		       string(" section that was extracted from" 
			      + largeFileName).c_str());
  /* Has to do something with RA, DEC ?  problem if RA = CRVAL1: CRVAL1 
     should not be changed */
  return true;
}


bool SubImage::MakeFits()
{
  return cut_in_fitsimage( &ReducedImage::FitsName, &ReducedImage::MakeFits);
}

bool SubImage::MakeWeight()
{
  bool cut_ok =  
    cut_in_fitsimage(&ReducedImage::FitsWeightName, &ReducedImage::MakeWeight);
  if (!cut_ok) return cut_ok;
  if (true)
    {// avoid file mode conflict
      FitsHeader head(FitsWeightName());
      // it would be better to generate and use here a type tester for SWARP.
      if (! (head.HasKey("SOFTNAME") && 
	     string(head.KeyVal("SOFTNAME")) == "SWarp")) return true;
    }
  // large image is from Swarp. Normalize weights so that sigma(im*sqrt(w) = 1
  // we need the image
  if (!MakeFits())
    {
      cerr << "cannot normalize weights from SWarp w/o corresponding image"
	   << endl;
      return false;
    }
  FitsImage im(FitsName());
  FitsImage w(FitsWeightName(), RW);
  double normFact = ImageAndWeightError(im ,w);
  normFact = 1./(normFact*normFact);
  w *= normFact;
  w.AddOrModKey("NORMFACT", normFact, 
		"norm factor applied to SWarp weights");
  return true;
}

bool SubImage::MakeCosmic()
{
  return cut_in_fitsimage(&ReducedImage::FitsCosmicName, &ReducedImage::MakeWeight);
}

bool SubImage::MakeSatellite()
{
  return cut_in_fitsimage(&ReducedImage::FitsSatelliteName, &ReducedImage::MakeWeight);
}

bool SubImage::MakeSatur()
{
  return cut_in_fitsimage(&ReducedImage::FitsSaturName, &ReducedImage::MakeSatur);
}

bool SubImage::MakeDead()
{
  return cut_in_fitsimage(&ReducedImage::FitsDeadName, &ReducedImage::MakeDead);
}



#include "sestar.h"

bool SubImage::MakeCatalog()
{
  if (HasCatalog()) return true;
  ReducedImage large(largeImageName);
  if (!large.HasCatalog())
    {
      cout << " Original large image " << largeImageName 
	   << "  has no catalog " << endl
	   << " Making Catalog of the subimage " << endl;
      return ReducedImage::MakeCatalog();
    }
  // else : extract
  SEStarList inList(large.CatalogName());

  for (SEStarIterator i = inList.begin(); i != inList.end();)
    {
      SEStar &s = **i;
      if (!subFrame.InFrame(s))
	i = inList.erase(i);
      else ++i;
    }

  GtransfoLinShift shift(-subFrame.xMin, -subFrame.yMin);
  inList.ApplyTransfo(shift);

  inList.write(CatalogName());
  return true;
}
	  
