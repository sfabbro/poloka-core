#include <iomanip>
#include <fstream>
#include "lightcurve.h"

//************************************************
//***************  LightCurve  *******************
//************************************************
bool ByIncreasingFlux(const FiducialStar *one, const FiducialStar *two)
{
  return ( (one->flux * one->night->photomRatio) < (two->flux * two->night->photomRatio));
}
bool ByIncreasingSignalToNoise(const FiducialStar *one, const FiducialStar *two)
{
  return ( (one->flux / sqrt(one->varflux)) < (two->flux / sqrt(two->varflux)) );
}

LightCurve::LightCurve(const NightList &allNights, const PhotStar &OneStar)
  : refstar(OneStar)
{ 
  for (NightCIterator it=allNights.begin(); it != allNights.end(); ++it)
    {
      const Night *night = *it;
      FiducialStar *star = new FiducialStar(OneStar, night);
      star->flux *= night->photomRatio;
      push_back(star);
    }
  IsStandard = false;
}

LightCurve::LightCurve(istream &Stream)
{
  read(Stream);
}

void LightCurve::AveragePosition()
{ 
  refstar.x = refstar.y = 0;
  double sumwx = 0 , sumwy = 0.;

  // loop over nights and weight average position
  // should really do correlation between x and y as well
  for (FiducialStarCIterator it = begin(); it != end(); ++it)
    {
      const FiducialStar *fs = *it;
      double wx = 1./ fs->varx;
      double wy = 1./ fs->vary;
      refstar.x += fs->x * wx;
      refstar.y += fs->y * wy;
      sumwx += wx;
      sumwy += wy;
    }  
  refstar.x /= sumwx;
  refstar.y /= sumwy;
  refstar.varx = 1. / sumwx;
  refstar.vary = 1. / sumwy;
}

void LightCurve::RelativeFluxes()
{
  // do the zero flux by a weighted average
  // and set all fluxes on same photometric scale
  // note we assume zero error on photometric ratio
  double zeroflux = 0;
  double sumw = 0;

  for (FiducialStarCIterator it = begin(); it != end(); ++it)
    {
      const FiducialStar *fs = *it;
      if (fs->night->IsZeroRef) 
	{
	  double weight = fs->night->photomRatio;
	  if (fs->varflux > 1e-10) weight /= fs->varflux;
	  zeroflux += fs->flux * weight;
	  sumw += weight * fs->night->photomRatio;
	}
    }

  double varzeroflux = 1./sumw;
  zeroflux *= varzeroflux;
  
  // subtract each instrumental flux with the scaled back zeroflux
  // we should compute the covariance between fluxes as well
  for (FiducialStarIterator it = begin(); it != end(); ++it) 
    {
      FiducialStar *fs = *it;
      fs->flux -= zeroflux * fs->night->photomRatio;
      fs->varflux += varzeroflux * fs->night->photomRatio * fs->night->photomRatio;
    }
}

void LightCurve::AbsoluteFluxes()
{
  
  // add the refstar flux phometrically scaled to each lightcurve point
  for (FiducialStarIterator it = begin(); it != end(); ++it) 
    {
      FiducialStar *fs = *it;
      fs->flux += refstar.flux * fs->night->photomRatio;
      fs->varflux += refstar.varflux * fs->night->photomRatio * fs->night->photomRatio;
    }

  refstar.flux = 0;
  refstar.varflux = 0;
}

//************************
// read and write routines
//************************

void LightCurve::dump(ostream& Stream) const
{
  int old = Stream.precision();
  Stream << "# julian_date " << endl
	 << "# flux " << endl
	 << "# eflux " << endl
	 << "# image_name " << endl
	 << "# end " << endl;
  Stream << setiosflags(ios::fixed);
  for (FiducialStarCIterator it = begin(); it != end(); ++it)
    {      
      const FiducialStar *fs = *it;
      if (fs->night) Stream << setw(14) << setprecision(2) << fs->night->julianDate;
      Stream << setw(15) << setprecision(3) << fs->flux 
	     << setw(15) << setprecision(3) << sqrt(fs->varflux);
      if (fs->night) Stream << "  " << fs->night->Name();
      Stream << endl;
    }
  Stream << resetiosflags(ios::fixed);
  Stream << setprecision(old);
}

void LightCurve::write(const string &FileName) const
{
  ofstream pr(FileName.c_str());
  write(pr);
  pr.close();
}

void LightCurve::write_header(ostream &Stream) const
{
  Stream << "# reference star" << endl;
  front()->WriteHeader(Stream);
}

void LightCurve::write(ostream &Stream) const
{
  Stream << setiosflags(ios::fixed);
  refstar.writen(Stream);
  for (FiducialStarCIterator it=begin(); it != end(); ++it) 
    {
      (*it)->writen(Stream); 
      Stream << endl;
    }
  Stream  << resetiosflags(ios::fixed);
}

void LightCurve::read(istream &Stream)
{
  string str;
  while (str != "end") Stream >> str;
  refstar.read(Stream,0);
  char c ;
  char buff[400];
  char *format = 0;
  while( Stream >> c )
    {
      Stream.unget() ;
      if ( (c == '#') )
        {
	  Stream.getline(buff,400);
	  char *p = strstr(buff,"format");
	  if (p) format = p + strlen("format");
        }
      else
	{
	  PhotStar *s = PhotStar::read(Stream, format);
	  if (!s) return;
	  FiducialStar *fs = new FiducialStar(*s);
	  delete s;
	  push_back(fs);
	}
    }
}
