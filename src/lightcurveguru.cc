#include <fstream>
#include <algorithm>
#include <functional>
#include "sestar.h"
#include "imagematch.h"
#include "lightcurveguru.h"
#include "transformedimage.h"
#include "fitsimage.h"
#include "reducedutils.h"
#include "fileutils.h"
#include "daophotutils.h"
#include "usnoutils.h"
#include "astroutils.h"
#include "wcsutils.h"

static vector<string> find_bands(const ReducedImageList &Images)
{
  // find all the bands
  vector<string> bands;
  for (ReducedImageCIterator it = Images.begin(); it != Images.end(); ++it) 
    {
      string curband = (*it)->Band();
      if (find(bands.begin(),bands.end(),curband) == bands.end()) 
	{
	  cout << " Found Band " << curband << endl;
	  bands.push_back(curband);
	}
    }  
  return bands;
}

struct IncreasingPointDist : public less<LightCurve> {
  Point pt;
  bool operator () (const LightCurve &s1, const LightCurve &s2)
  {return (s1.refstar.Distance(pt) < s2.refstar.Distance(pt));}
};

static void matchBuilders(const LightCurveBuilder &B1, const LightCurveBuilder &B2, 
			  BaseStarList &B2ListInB1)
{
  Gtransfo *b1Tob2 = NULL;
  Gtransfo *b2Tob1 = NULL;
  ImageListMatch(*B1.Reference, *B2.Reference, b1Tob2, b2Tob1);
  for (LightCurveCIterator fid = B2.begin(); fid != B2.end(); ++fid)
    B2ListInB1.push_back(new BaseStar(b2Tob1->apply(fid->refstar), fid->refstar.flux));
  delete b1Tob2;
  delete b2Tob1;
}

static const ReducedImage* CheckAndRegister(ReducedImageList &Images, ReducedImageList &AlignedImages)
{
  cout << " Checking for " << Images.size() << " images overlap and reduced status " << endl;
  ReducedImageList overlapImages(false);
  FilterWithOverlap(Images, *Images.front(), overlapImages);
  
  cout << " Select a geometric reference " << endl;
  const ReducedImage *geomRef = BestResolutionReference(overlapImages);
  if (!geomRef) return NULL;
   //UsnoProcess(*geomRef,true);
  
  cout << " Register images " << endl;
  ImagesAlign(overlapImages, *geomRef, AlignedImages, DoFits | DoCatalog | DoSatur | DoWeight);
  return geomRef;
}


//----------------------------------------------------------//

LightCurveGuru::~LightCurveGuru()
{
  if (coordRef) delete coordRef; // useless if it becomes ReducedImageRef
}

/*! \page lightfile Syntax of the "lightfile"

A lightfile consists of SN coordinates and images. The DATEZERO
is used to set the lightcurve baseline. It includes all images
for all filters.

\code
DATEZERO
# Here you put the dates min and max within which you are sure
# there is supernova light (DD MM YYYY  DD MM YYYY )
20 12 2000 20 12 2001

COORDREF
# The name of the DbImage where SN coordinates are expressed
MyCoordImage

COORD
# The relative coordinates of the supernova in the COORDREF mentioned above
1061.45 1048.98

IMAGES
# The list of all DbImages in all filters
# you want to build a lightcurve from
Image1_R
Image2_R
Image1_I
Image2_I
\endcode

The COORDREF can be any reduced image, which does not need to be processed for
lightcurve building (typically a subtraction). 

All images in this file should: 
\arg be flatfielded
\arg optionally background subtracted (but preferred)
\arg with a GAIN of 1 (TOADGAIN key)
\arg with proper saturation level (SATURLEV key)
\arg with a catalog produced typically by the make_catalog application
\arg in the DbImage system
\arg optionally with a WCS transfo (but preferred)
*/

bool LightCurveGuru::read(const string &FileName)
{
  cout << " Reading lightcurve file " << endl;
  ifstream rd(FileName.c_str());
  if (!rd) 
    {
      cerr << " LightCurveGuru: file " << FileName << " does not exist" << endl;
      return false;
    }
  string line;
  bool datin=false, coordin=false, imin=false, coordref=false;

  while (!rd.eof())
    {
      getline(rd,line);

      // read a comment
      if (line[0] == '#' || line.length() == 0) {continue;}

      if (line.find("DATEZERO") != line.npos) 
	{datin=true;  coordref=false; imin=false;  coordin=false; continue;}

      if (line.find("COORDREF") != line.npos) 
	{datin=false; coordref=true;  imin=false;  coordin=false; continue;}

      if (line.find("IMAGES") != line.npos)   
	{datin=false; coordref=false; imin=true;   coordin=false; continue;}

      if (line.find("COORD") != line.npos)    
	{datin=false; coordref=false; imin=false;  coordin=true;  continue;}


      // read dates min and max
      if (datin)
	{
	  vector<string> datestring;
	  DecomposeString(datestring,line);
	  if (datestring.size() != 6) 
	    {
	      cerr << " LightCurveGuru::read bad date format " << endl;
	      return false;
	    }
	  jd_inf = JulianDay(atoi(datestring[0].c_str()),
			     atoi(datestring[1].c_str()),
			     atoi(datestring[2].c_str()));

	  jd_sup = JulianDay(atoi(datestring[3].c_str()),
			     atoi(datestring[4].c_str()),
			     atoi(datestring[5].c_str()));
	  // cout << jd_inf << " " << jd_sup << endl;
	  continue;
	}

      // read coordinate reference image
      if (coordref)
	{
	  if (coordRef) delete coordRef;
	  RemovePattern(line," ");
	  coordRef = new ReducedImage(line);
	  if (!coordRef->ActuallyReduced()) 
	    {
	      cerr << " WARNING: coordinate reference " << line << " is not reduced" << endl;
	    }
	  continue;
	}

      // read coordinates of each supernova
      if (coordin)
	{
	  vector<string> coordstring;
	  DecomposeString(coordstring,line);
	  if (coordstring.size() != 2) 
	    {
	      cerr << " LightCurveGuru::read Bad coordinates format " << endl;
	      return false;
	    }
	  double xc = atof(coordstring[0].c_str());
	  double yc = atof(coordstring[1].c_str());
	  //cout << " Coordinates : " << xc << ' ' << yc << endl;
	  if (xc < 0 || yc < 0)
	    {
	      cerr << " Negative coordinates !!! "<< xc << " " << yc << endl;
	      return false;
	    }
	  Supernovae.push_back(Point(xc,yc));
	  continue;
	}

      // get images
      if (imin)
	{
	  RemovePattern(line," ");
	  ReducedImage *red = new ReducedImage(line);
	  if (red->ActuallyReduced()) AllImages.push_back(red);	  
	  else
	    {
	      cerr << " '" << line << "' is not reduced " << endl;
	      delete red;
	      return false;
	    }
	}
    }  
  rd.close();

  if (AllImages.size() == 0) 
    {
      cerr << " WARNING: could not figure out any image name in " <<  FileName << endl;
      return false;
   }

  cout << " Read " << Supernovae.size() << " stars and " << AllImages.size() << " images" << endl;
  return true;
}

bool LightCurveGuru::check_coordinates() const
{
  // a quick check on coordinates, without changing them

  cout << " Quick check for candidates coordinates " << endl;  

  SEStarList coordList(coordRef->CatalogName());
  bool status = true;
  for (unsigned i=0; i<Supernovae.size(); ++i)
    {
      SEStar *closest = coordList.FindClosest(Supernovae[i]);
      double dist = closest->Distance(Supernovae[i]);
      cout << " Found candidate on coordinate reference " << coordRef->Name() 
	   << " x=" << closest->x  << " y=" << closest->y << endl;
      if (dist > 5) 
	{
	  cerr << " WARNING: Are you sure of this object coordinates : " << Supernovae[i] << endl
	       << " I found it " << dist << " pixels away from what you gave me " << endl;
	}
      if (dist > 20) 
	{
	  cerr << " FAILURE: I found it " << dist 
	       << " pixels away from what you gave me " << endl;
	  status = false;
	}
    }

  return status;
}

vector<Point> LightCurveGuru::match_candidates(const ReducedImage &WithImage) const
{
  cout << " Find the candidate stars on " << WithImage.Name() << endl;

  Gtransfo* one2two = NULL; 
  Gtransfo* two2one = NULL;
  ImageListMatch(*coordRef, WithImage, one2two, two2one);
  vector<Point> found;
  for (unsigned i=0; i<Supernovae.size(); ++i)
    {  
      //cout << " Coordinates on " << coordRef->Name() << " :" << Supernova << endl;
      found.push_back(one2two->apply(Supernovae[i]));      
      //cout << " Coordinates on " << WithImage.Name() << " :" << found << endl;
    }
  delete one2two;
  delete two2one;

  return found;
}

bool LightCurveGuru::MakeNights()
{
  vector<string> bands = find_bands(AllImages);
  size_t nbands = bands.size();
  if (nbands==0) 
    {
      cerr << " LightCurveGuru did not find any bands" << endl;
      return false;
    }
  
  for (size_t band=0; band<nbands; ++band)// loop on bands
    {
      cout << " Doing band " << bands[band] << endl;
      string outstr = "nights_"+bands[band]+".dat";

      ReducedImageList bandImages(false);
      for (ReducedImageCIterator it = AllImages.begin(); it!=AllImages.end(); ++it)
	{      
	  const ReducedImage *rim = *it;
	  if (rim->Band() == bands[band]) bandImages.push_back(rim);
	}
      ReducedImageList alignedImages;
      const ReducedImage *geomRef = CheckAndRegister(bandImages, alignedImages);

      cout << " Preparing night images " << endl;
      sort(alignedImages.begin(), alignedImages.end(), IncreasingSeeing);
      LightCurveBuilder current(alignedImages);
  
      cout << " Initialize nights and produce PSFs" << endl;
      FitsHeader geomRefHeader(geomRef->FitsName());
      double equinox = geomRefHeader.KeyVal("TOADEQUI");

      for (NightIterator it = current.Nights.begin(); it != current.Nights.end(); ++it)
	{
	  Night *night = *it;
	  MakeDaoPsfCat(*night, true, false, false, false);
	  FitsHeader nightHeader(night->FitsName(),RW);
	  CopyWCS(geomRefHeader, nightHeader);
	  UpdateRaDec(nightHeader);
	  nightHeader.AddOrModKey("TOADEQUI",equinox,"Equinox updated with RA and DEC");
	  night->init();
	  night->IsZeroRef = (night->julianDate <= jd_inf) || 
	                     (night->julianDate >= jd_sup);
	  night->IsPhotometric = false; // to be changed
	}
      sort(current.Nights.begin(), current.Nights.end(), IncreasingNightSeeing);
      current.Reference = new Night(*current.Nights.front());

      cout << " Initialize lightcurves" << endl;
      bool status;
      if(StdCatalogue.size()==0)
	status = current.Init(match_candidates(*current.Reference));
      else
	status = current.Init(match_candidates(*current.Reference),StdCatalogue);
      if (!status)
	{
	  cerr << " LightCurveBuilder::Init failure for band " << bands[band] << endl;
	  return false;
	}
      BandBuilders.push_back(current);

      cout << " Nights are built for band " << bands[band] << endl;
      current.Nights.dump();
      ofstream out(outstr.c_str());
      out << current.Reference->Name() << endl;
      current.Nights.write(out);
      out.close();

    }

  return true;
}


bool LightCurveGuru::SelectFidRef()
{
  cout << " Looking for fiducials reference " << endl;
  fidRef = &BandBuilders[0];
  for (size_t i=0; i<BandBuilders.size();++i)
    {
      if (BandBuilders[i].size() > fidRef->size())
	{
	  fidRef = &BandBuilders[i];
	}
    }

  for (size_t i=0; i<fidRef->nsn;++i) Supernovae[i] = (*fidRef)[i].refstar;
  if (coordRef) delete coordRef;
  coordRef = fidRef->Reference->Clone();
  cout << " Found " << fidRef->Reference->Name() << " as reference for fiducials" << endl; 
  return true;
}

bool LightCurveGuru::MatchFiducials()
{
  if (fidRef->size() == 0) return false;

  for (size_t i=0; i<BandBuilders.size(); ++i)
    {
      LightCurveBuilder &current = BandBuilders[i];

      // match fidRef with current builder reference
      Gtransfo *curToref = NULL;
      Gtransfo *refTocur = NULL;
      ImageListMatch(*current.Reference, *fidRef->Reference, curToref, refTocur);
      delete curToref;

      // remove objects near the edges
      Frame commonFrame = current.Nights.CommonFrame();
      sort(current.Nights.begin(), current.Nights.end(), IncreasingSeeing);
      double mindist = current.Nights.back()->seeing*2;
      cout << " Band " << current.BandName << " has " 
	 << current.size() << " fiducials " << endl;

      // loop over all fiducials except SN
      LightCurveIterator it = fidRef->begin(); ++it;
      while (it != fidRef->end())
	{
	  Point curStar = refTocur->apply(it->refstar);
	  if (commonFrame.MinDistToEdges(curStar) > mindist) ++it;
	  else it = fidRef->erase(it);
	}
      delete refTocur;
      cout << " After matching bands " << fidRef->BandName << " and " << current.BandName 
	   << " : kept " << fidRef->size() << " fiducials "<< endl;
    }
  
  SelectGoodFiducials();
  
  //size_t nfid = min(lcParams.maxFiducials+Supernovae.size(),fidRef->size());
  //fidRef->resize(nfid);
  cout << " Final number of fiducials : " << fidRef->size()-Supernovae.size() << endl;
  
  return true;
}

bool LightCurveGuru::SelectGoodFiducials() {
  
  // pointer to the first "real" fiducial in fidRef (others are supernovae)
  LightCurveIterator firstfid = fidRef->begin();
  for(unsigned int i=0;i<Supernovae.size();i++)
    ++firstfid;


  // let's first check if no star is saturated in the reference image
  //=================================================================
  SEStarList curList(fidRef->Reference->CatalogName());
  LightCurveIterator fid = firstfid;
  double SaturLev = fidRef->Reference->Saturation();
  while (fid != fidRef->end()) {
    SEStar *star = curList.FindClosest(fid->refstar);
    bool isok = (star->FlagBad() == 0);
    isok &= (star->Flag() < 3);
    isok &= ( star->Fluxmax()+star->Fond() < SaturLev );
    if(isok)
      ++fid;
    else
      fid = fidRef->erase(fid);
  }
  cout << " After removing saturated stars, keep : " << fidRef->size()-Supernovae.size() << " fiducials "<< endl;
  

  // Here we try to select fiducials which are close to supernovae.
  // We want about lcParams.maxFiducials fiducials per supernova
  // So, we try to find a cut in distance to supernova to match
  // this number ...
  //============================================================
  cout << "Now keep <= " << lcParams.maxFiducials << " fiducials per supernova" << endl;
  
  // first guess
  Frame refframe = fidRef->Reference->UsablePart();
  float area = refframe.Nx()*refframe.Ny(); // in pixels2
  float fiddensity = (fidRef->size()-Supernovae.size())/area;
  float distance2cut = lcParams.maxFiducials/fiddensity/M_PI; //pixel2
  float distance2cutmax = refframe.Nx()*refframe.Nx()+refframe.Ny()*refframe.Ny();
  
 
  
  // increase distance2cut if less than lcParams.maxFiducials fiducials for one supernova
  unsigned int nsn = 0;
  LightCurveIterator sn = fidRef->begin();
  unsigned int ok = 0;
  while (sn != fidRef->end() && nsn <Supernovae.size()) { // loop on supernovae
    Point supernova_pos = sn->refstar;
    nsn ++;
    
    ok = 0;
    while(ok< lcParams.maxFiducials || distance2cut>distance2cutmax) {
      ok = 0;
      for( LightCurveIterator fid = firstfid ; fid != fidRef->end(); ++fid) {
	if(supernova_pos.Dist2(fid->refstar)<distance2cut) 
	  ok ++;
      }
      if(ok>=lcParams.maxFiducials)
	break;
      
      if(ok>0)
	distance2cut *= ((float)lcParams.maxFiducials)/ok;
      else
	distance2cut *= 2;
    }
    if(ok<lcParams.maxFiducials) 
      cout << "WARNING in LightCurveGuru::SelectGoodFiducials, only " << ok << "fiducials for supernova " << nsn << endl;
  }
  
  // now remove other fiducials 
  fid = firstfid;
  while(fid != fidRef->end()) {
    Point fid_pos = fid->refstar;
    
    sn = fidRef->begin();
    bool removeit = true;
    for (unsigned int nsn = 0 ; nsn<Supernovae.size(); nsn++) {
      if(fid_pos.Dist2(sn->refstar)<distance2cut) {
	removeit = false;
	break;
      } 
      ++sn;
    }
    if(removeit)
      fid = fidRef->erase(fid);
    else
      ++fid;
  }
  cout << " After selecting good fiducials, keep : " << fidRef->size()-Supernovae.size() << " fiducials "<< endl;
  cout << " Distance cut is " << sqrt(distance2cut) << " pixels" << endl;
  return true;
}

bool LightCurveGuru::MonopolizeCPU()
{
  /*  if (!check_coordinates()) 
    {
      cerr << " LightCurveGuru::MonopolizeCPU : FAILURE checking coordinates " << endl;
      return false;
    }
  */
  if (!MakeNights()) 
    {
      cerr << " LightCurveGuru::MonopolizeCPU : FAILURE building nights " << endl;
      return false;
    }

  if (!SelectFidRef())
    {
      cerr << " LightCurveGuru::MonopolizeCPU : FAILURE selecting fiducial reference image " << endl;
      return false;
    }

  if (!MatchFiducials())
    {
      cerr << " LightCurveGuru::MonopolizeCPU : FAILURE matching fiducials " << endl;
      return false;
    }

  for (size_t i=0; i<BandBuilders.size(); ++i)
    {
      LightCurveBuilder &current = BandBuilders[i];
      cout << " Building lightcurves for band " << current.BandName << endl;

      BaseStarList currentList;
      matchBuilders(current, *fidRef, currentList);
      current.erase(current.begin(), current.end());
      current.InitFiducials(currentList);
      current.write("cat");
      current.BuildKernelsAndSubs();      
      current.SubAperPhotometry(1);
      current.write("sub");
      MakePrecisePsf(*current.Reference); // make sure you do it after the kernel stuff
      if (FileExists("fiducials_simfit_"+current.BandName+".list")) continue;
      current.SimFitPhotometry();
      current.write("simfit");
    }

  return true;
}

