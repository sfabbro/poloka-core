#include <iomanip>

#include "reducedimage.h"
#include "reducedutils.h"
#include "starmatch.h"
#include "imagesum.h"
#include "imagesubtraction.h"
#include "night.h"
#include "sestar.h"
#include "fileutils.h"

//*****************************Night*******************************************

Night::Night() : 
  saturation(0), photomRatio(1), kernelChi2(0),julianDate(0), exposure(0), zeroPoint(0), 
  seeing(0), sigmaX(0), sigmaY(0), thetaXY(0), skyLevel(0), backLevel(0),
  sigmaSky(0), gain(0), flatError(0), readoutNoise(0), profileError(0), 
  IsZeroRef(false), IsPhotometric(false) {}


void Night::init()
{
  photomRatio = 1;
  saturation = kernelChi2 = julianDate = 0;
  exposure = zeroPoint = 0;
  seeing = skyLevel = backLevel = 0;
  sigmaSky = gain = flatError = 0;
  readoutNoise = profileError = 0;
  if (FileExists(FitsName()))
    {
      saturation = Saturation();
      julianDate = JulianDate();
      exposure = Exposure();
      zeroPoint = ZeroPoint();
      seeing = Seeing();
      GetPsfShapeParams(sigmaX, sigmaY, thetaXY);
      skyLevel  = OriginalSkyLevel();
      backLevel = BackLevel();
      sigmaSky = SigmaBack();
      gain = Gain();
      flatError = FlatFieldNoise();
      readoutNoise = ReadoutNoise();
      profileError = ProfileError();
    }

  IsZeroRef = IsPhotometric = false;
}

Night::Night(const string &ImageName) : ReducedImage(ImageName)
{
  init();
}


ReducedImage* Night::Clone() const
{
  return new Night(*this);
}

bool Night::IsOK() const
{
  return
    (
     (seeing > 0) &&
     (saturation > 0) &&
     (skyLevel > 0) && 
     (backLevel < saturation) &&
     (exposure > 0) &&
     (julianDate > 0) &&
     // We comment this line because, with the new stacking, the gain of the stacked image is not necessary 1
     //(gain == 1) &&
     (photomRatio > 0)
     );
}

static void night_not_ok_header()
{
  cerr << "NAME  SEEING SKY BACK SATUR EXPO JD GAIN RATIO" << endl;
  cerr << " -----------------------------------------------\n";
}

static void night_not_ok(const Night *night)
{
  cerr << night->Name() << " " 
       << night->seeing << " "
       << night->skyLevel << " "
       << night->backLevel << " " 
       << night->saturation << " "
       << night->exposure << " "
       << night->julianDate << " "
       << night->gain << " "
       << night->photomRatio << endl;
}


void Night::dump_header(ostream &Stream) const
{
  Stream << "     JD        Expo   Seeing   Sigma   Ratio  Chi2     Image" << endl
	 << "   ------------------------------------------------------------"
	 << endl;
}

void Night::dump(ostream &Stream) const
{
  Stream << setiosflags(ios::fixed);
  Stream << setw(11) << setprecision(2) << julianDate 
	 << setw(8)  << setprecision(0) << exposure 
	 << setw(7)  << setprecision(3) << seeing 
	 << setw(8)  << setprecision(2) << sigmaSky
	 << setw(7)  << setprecision(3) << photomRatio 
	 << setw(8)  << setprecision(2) << kernelChi2
	 << " " << Name() ;
}

//***************************** NightList *******************************************

NightList::NightList(const string &FileName)
{
  read(FileName);
}


NightList::NightList(const ReducedImageList &Images)
{
  cout << " Will stack the images for you" << endl;
  vector<NightSet> allnights;
  ArrangeByNight(Images, allnights);
  unsigned nNights = allnights.size();
  for (unsigned i=0; i<nNights ; ++i)
    { 
      NightSet &nightset = allnights[i];
      string nightName = nightset.GenericName();
      Night *current = new Night(nightName);
      if (!current->ActuallyReduced())
	{
	  ImageSum sum(nightName, nightset);
	  sum.Execute(DoFits | DoCatalog | DoSatur | DoWeight);
	  delete current;
	  current = new Night(nightName);
	}
      push_back(current);
    }
  cout << " Your " << size() << " nights are stacked" << endl;
}

void NightList::InitPhotomRatio(const ReducedImage &Reference)
{
  SEStarList reflist(Reference.CatalogName());
  GtransfoIdentity *ident = new GtransfoIdentity();
  for (NightIterator ni=begin(); ni != end(); ++ni)
    {
      Night *current = *ni;
      SEStarList curlist(current->CatalogName());
      double error;
      current->photomRatio = QuickPhotomRatio(reflist, curlist, error, ident);
    }
  delete ident;
}

Frame NightList::CommonFrame()
{
  sort(begin(), end(), DecreasingNightArea);
  Frame frame_common(front()->UsablePart());
  for (NightCIterator ni=begin(); ni != end(); ++ni)
    {
      const ReducedImage *current = *ni;
      frame_common *= Frame(current->UsablePart());    
    }
  cout << " Common frame for NightList: " << frame_common;
  return frame_common;
}

bool NightList::Init()
{
  cout << " Checking nights for a good start" << endl;
  bool status = true;
  for (NightIterator ni = begin(); ni != end(); ++ni)
    {
      Night *night = *ni;
      if (!night->IsOK())
	{
	  if (status) 
	    {
	      cerr << " WARNING! NightList::Init something wrong with following nights" << endl;
	      night_not_ok_header();
	      status = false;
	    }
	  night_not_ok(night);
	}
    }
  return status;
}

void NightList::dump(ostream &Stream) const
{
  front()->dump_header(Stream);
  for (NightCIterator ni = begin(); ni != end(); ++ni) 
    {
      (*ni)->dump(Stream);
      Stream << endl;
    }

}

void NightList::write(ostream &Stream) const
{
  for (NightCIterator ni = begin(); ni != end(); ++ni) 
    {
      Stream << (*ni)->Name() << " ";
      if ((*ni)->IsPhotometric) Stream << " p ";
      if ((*ni)->IsZeroRef) Stream << " z ";
      Stream << endl;
    }
}

void NightList::write(const string &FileName) const
{ 
  ofstream out(FileName.c_str()); 
  write(out);
  out.close();
}

bool NightList::read(istream &Stream)
{
  while (!Stream.eof())
    {
      string line;
      Stream >> line;
      if (line.length() > 1) 
	{
	  cout << " Found built night " << line << endl;
	  Night *night = new Night(line);
	  bool nextone = true;
	  char c;
	  while (nextone && !Stream.eof())
	    {
	      Stream >> c;
	      switch (c)
		{
		case 'z': night->IsZeroRef = true; break;
		case 'p': night->IsPhotometric = true; break;
		default: nextone = false;
		}
	    }
	  Stream.putback(c);
	  push_back(night);
	}
    }
  return true;
}

bool NightList::read(const string &FileName)
{
  ifstream Stream(FileName.c_str());
  if (!Stream)
    {
      cerr << " NightList cannot open :" << FileName << endl;
      return false;
    }
  return read(Stream);
}

void NightList::FilterAllObjects(BaseStarList &RefList, const Frame &RefFrame, const double MinAssDist) const
{
  // Cleanup all lists
  /* loop over nights on a subset of the best seeing fiducial list
     and make sure the fiducial is present and not saturated on all nights */ 
  cout << " Starting with " << RefList.size() << " objects" << endl;

   // removing doublons (0.1<dist<1 pixels)
   for (BaseStarIterator it = RefList.begin(); it != RefList.end(); )
     {
	BaseStar *refstar = *it;
	BaseStar *curstar = RefList.ClosestNeighbor(*refstar);
	if (refstar->Distance(*curstar) < 1)
	  {
	     it = RefList.erase(it);
	  }
	else ++it;
	
     }
   
  // loop on images
  for (NightCIterator ni=begin(); ni!=end(); ++ni)
    {
      const Night *night = *ni;
      SEStarList curlist(night->CatalogName());
      double minEdgesDist = night->seeing*2;

      // loop on objects
      for (BaseStarIterator it = RefList.begin(); it != RefList.end(); )
	{
	  BaseStar *refstar = *it;
	  SEStar *curstar = curlist.FindClosest(*refstar);
	  double dist = curstar->Distance(*refstar);
	  
	  if ((RefFrame.InFrame(*curstar))
	      && (RefFrame.MinDistToEdges(*curstar) > minEdgesDist)
	      && (dist < MinAssDist)
	     )
	    {
	      ++it;
	    }
	  else 
	    {
	      it = RefList.erase(it);
	    }
	}
     cout << " Left with " << RefList.size() 
	 << " after " << night->Name() << endl;  
    }
}

bool IncreasingNightSeeing(const Night *one, const Night *two)
{ 
  return (one->seeing < two->seeing);
}

bool DecreasingNightRatio(const Night *one, const Night *two)
{ 
  return (one->photomRatio > two->photomRatio);
}

bool IncreasingNightJulian(const Night *one, const Night *two)
{ 
  return (one->julianDate < two->julianDate);
}

bool DecreasingNightSaturation(const Night *one, const Night *two)
{ 
  return ((one->saturation - one->backLevel) < (two->saturation - two->backLevel));
}

bool DecreasingNightArea(const Night *one, const Night *two)
{ 
  return (one->UsablePart().Area() < two->UsablePart().Area());
}

bool SameSeeing(const Night &Reference, const Night &Current, const double &Tol)
{
  return (fabs(Reference.seeing - Current.seeing) < Tol*Reference.seeing);
}
