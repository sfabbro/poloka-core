#include <iomanip>
#include <fstream>
#include <reducedimage.h>

#include "daophotoption.h"

//***************************** Stream Utilities  *******************
ostream& operator << (ostream &stream, const DaoOption& opt)
{
  size_t oldp = stream.precision();
  stream << setiosflags(ios::fixed);  
  stream << string(27-opt.name.length(),' ')  << opt.name
	 << " =" << setw(9) << setprecision(2) << opt.value;
  stream << resetiosflags(ios::fixed) << setprecision(oldp);
  return stream;
}

ostream& writeDaoOption(const DaoOption& opt, ostream &stream)
{
  stream << opt.ShortName() << "=" << opt.Value() << endl;
  return stream;
}

ostream& operator << (ostream &stream, const DaophotOptions& Options) 
{
  vector<DaoOption>::const_iterator opend = Options.end();
  for (vector<DaoOption>::const_iterator it = Options.begin(); it != opend; it +=2) 
    {
      stream << *it << "    " << *(it+1) << endl;
    }
  return stream;
}

//***************************** DaoOption ******************************
DaoOption::DaoOption(const string& Name, const float& Min, 
		     const float& Max, const float& Default)
  : name(Name), omin(Min), omax(Max), odefault(Default), value(Default) {}

bool DaoOption::SetValue(const float& Value)
{
  if (Value >= omin && Value <= omax) { value=Value; return true;}
  cout << " DaoOption::SetValue() : Replacing " << ShortName() << "=" << Value
       << " by default value = " << odefault << endl;
  value = odefault;
  return false;
}

string DaoOption::ShortName() const 
{
  return name.substr(0,2);
}

//***************************** DaophotOptions ****************************
DaophotOptions::DaophotOptions()
{
  init_daophot();
  read();
}

DaophotOptions::DaophotOptions(const ReducedImage& Rim)
{

  cout << " DaophotOptions::DaophotOptions() : Initialize options \n";

  init_daophot();

  cout << " DaophotOptions::DaophotOptions() : settings from " << Rim.Name() << endl;

  (*this)[ReadNoise].SetValue(Rim.ReadoutNoise());
  (*this)[Gain].SetValue(Rim.Gain());
  (*this)[AduHighDatum].SetValue(Rim.Saturation());
  (*this)[PercError].SetValue(Rim.FlatFieldNoise());
  (*this)[ProfError].SetValue(Rim.ProfileError());

  const float fwhm = Rim.Seeing() * 2.3548;
  (*this)[Fwhm].SetValue(fwhm);
  (*this)[FitRadius].SetValue(fwhm);
  (*this)[PsfRadius].SetValue(4.0 * fwhm);
  (*this)[InnerSky].SetValue(1.5 * fwhm);
  (*this)[OutterSky].SetValue(7.0 * fwhm);
}

DaoOption& DaophotOptions::operator [] (const DaoOptionEnum OptEnum)
{ 
  return vector<DaoOption>::operator[](static_cast<size_t>(OptEnum));
}

const DaoOption& DaophotOptions::operator [] (const DaoOptionEnum OptEnum) const
{ 
  return vector<DaoOption>::operator[](static_cast<const size_t>(OptEnum));
}

void DaophotOptions::init_daophot()
{
  clear();
  resize(26);
  
  //DAOPHOT
  (*this)[ReadNoise]      = DaoOption("READ NOISE (ADU; 1 frame)" , 1e-20, 1e20, 2.   );
  (*this)[Gain]           = DaoOption("GAIN (e-/ADU; 1 frame)"    , 1e-20, 1e20, 1.   );
  (*this)[SigmaLowDatum]  = DaoOption("LOW GOOD DATUM (in sigmas)", 0.   , 1e20, 7.   );
  (*this)[AduHighDatum]   = DaoOption("HIGH GOOD DATUM (in ADU)"  , 0.   , 1e20, 1.5e5);
  (*this)[Fwhm]           = DaoOption("FWHM OF OBJECT"            , 0.2  , 20. , 2.5  );
  (*this)[SigmaThreshold] = DaoOption("THRESHOLD (in sigmas)"     , 0.   , 1e20, 4.   );
  (*this)[LowSharpness]   = DaoOption("LS (LOW SHARPNESS CUTOFF)" , 0.   , 1.  , 0.2  );
  (*this)[HighSharpness]  = DaoOption("HS (HIGH SHARPNESS CUTOFF)", 0.6  , 2.  , 1.   );
  (*this)[LowRoundness]   = DaoOption("LR (LOW ROUNDNESS CUTOFF)" , -2.  , 0.  , -1.  );
  (*this)[HighRoundness]  = DaoOption("HR (HIGH ROUNDNESS CUTOFF)", 0.   , 2.  , 1.   );
  (*this)[WatchProgress]  = DaoOption("WATCH PROGRESS"            , -2.  , 2.  , 0.   );
  (*this)[FitRadius]      = DaoOption("FITTING RADIUS"            , 1.   , 30. , 2.5  );
  (*this)[PsfRadius]      = DaoOption("PSF RADIUS"                , 5.   , 51. , 11.  );
  (*this)[VariablePsf]    = DaoOption("VARIABLE PSF"              , -1.5 , 3.5 , -1.  );
  (*this)[SkyEstimator]   = DaoOption("SKY ESTIMATOR"             , -0.5 , 3.5 , 0.   );
  (*this)[AnalyticPsf]    = DaoOption("ANALYTIC MODEL PSF"        , -6.5 , 6.5 , 3.   );
  (*this)[ExtraCleanPsf]  = DaoOption("EXTRA PSF CLEANING PASSES" , 0.   , 9.5 , 5.   );
  (*this)[UseSaturated]   = DaoOption("USE SATURATED PSF STARS"   , 0.   , 1.  , 0.   );
  (*this)[PercError]      = DaoOption("PERCENT ERROR (in %)"      , 0.   , 100., 0.75 );
  (*this)[ProfError]      = DaoOption("PROFILE ERROR (in %)"      , 0.   , 100., 3.   );
  (*this)[ClipExponent]   = DaoOption("CE (CLIPPING EXPONENT)"    , 0.   , 8.  , 6.   );
  (*this)[Recentroid]     = DaoOption("REDETERMINE CENTROIDS"     , 0.   , 1.  , 1.   );
  (*this)[ClipRange]      = DaoOption("CR (CLIPPING RANGE)"       , 0.   , 10. , 2.5  );
  (*this)[MaxGroup]       = DaoOption("MAXIMUM GROUP SIZE"        , 1.   , 100., 70.  );
  (*this)[InnerSky]       = DaoOption("IS (INNER SKY RADIUS)"     , 0.   , 35. , 0.   );
  (*this)[OutterSky]      = DaoOption("OS (OUTER SKY RADIUS)"     , 0.   , 100., 0.   );
}

float* DaophotOptions::GetArray(const int NoptDaophot) const 
{
  float *opt = new float[NoptDaophot];
  for (int i=0; i<NoptDaophot; ++i) 
    {      
      const DaoOption &tmp = (*this)[static_cast<DaoOptionEnum>(i)];
      opt[i] = tmp.Value();
    }
  return opt;
}

bool DaophotOptions::write(const string& FileName) const 
{
  ofstream stream(FileName.c_str());
  if (!stream) 
    {
      cerr << "DaophotOptions::write() : Error opening " << FileName << endl;
      return false;
    }
  for(vector<DaoOption>::const_iterator it = begin(); it != end(); ++it) 
    writeDaoOption(*it, stream);
  return true;
}

bool DaophotOptions::read(const string& FileName) 
{
  ifstream stream(FileName.c_str());
  if (!stream) 
    {
      cerr << "DaophotOptions::read() : Error opening " << FileName << endl;
      return false;
    }
  string line;
  if (empty()) init_daophot();
  while (getline(stream, line)) 
    {
      vector<string> s;
      DecomposeString(s,line,'=');
      if (s.size() != 2)
	{
	  cerr << "DaophotOptions::read() : Error reading line '" << line 
	       << "' in " << FileName << endl;
	  continue;
	}

      vector<DaoOption>::iterator it = find(begin(), end(), s[0].substr(0,2));
      it->SetValue(atof(s[1].c_str()));
    }
  return true;
}

//***************************** Aperture options ******************************  
bool WriteDaoAperOpt(vector<double> &Radius, double& InnerSky, 
		     double& OutterSky, const string& FileName)
{
  ofstream stream(FileName.c_str());
  if (!stream) 
    {
      cerr << "WriteDaoAperOpt() : Error opening file " << FileName << endl;
      return false;
    }

  if (InnerSky > OutterSky)
    {
      cerr << "WriteDaoAperOpt() : Warning: inner sky bigger than outter sky. Swapping them. " << endl;
      swap(InnerSky, OutterSky);
    }
  else if (InnerSky == OutterSky)
    {
      cerr << "WriteDaoAperOpt() : Error: inner sky = outter sky. " << endl;
      return false;
    }

  size_t nrad = Radius.size();
  if (nrad == 0) 
    {
      cerr << "WriteDaoAperOpt() : Error: no radius specified " << endl;
      return false;
    }

  // sort vector by increasing radius
  sort(Radius.begin(), Radius.end());

  const string sRadius = "123456789ABC";
  size_t nsrad = sRadius.length();
  if (nrad > nsrad) 
    {
      cerr << "WriteDaoAperOpt() : Warning: max number of radii is " << nsrad << endl;
      cerr << "WriteDaoAperOpt() : Warning: will not do the others. \n";      
      nrad = nsrad;
    }

  for (size_t i=0; i<nsrad; ++i)
    stream << "A" << sRadius[i] << "=" << Radius[i] << endl;
  stream << "IS=" << InnerSky << endl
	 << "OS=" << OutterSky << endl ;

  return true;
}

