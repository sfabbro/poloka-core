// This may look like C code, but it is really -*- C++ -*-
//  daophotoption.h
//
// Last change: Time-stamp: <07 Dec 2003 at 19:52:07 by Sebastien Fabbro>
//
//
#ifndef DAOPHOTOPTION__H
#define DAOPHOTOPTION__H

#include <vector>

//! the option enum to index the vector of options
typedef enum {ReadNoise=0,
	      Gain,
	      SigmaLowDatum,
	      AduHighDatum,
	      Fwhm,
	      SigmaThreshold,
	      LowSharpness,
	      HighSharpness,
	      LowRoundness,
	      HighRoundness,
	      WatchProgress,
	      FitRadius,
	      PsfRadius,
	      VariablePsf,
	      SkyEstimator,
	      AnalyticPsf,
	      ExtraCleanPsf,
	      UseSaturated,
	      PercError,
	      ProfError,
	      ClipExponent,
	      Recentroid,
	      ClipRange,
	      MaxGroup,
	      InnerSky,
	      OutterSky} DaoOptionEnum;


//! a small class to hold options name, limits and value
class DaoOption { 

  string name;     // option name, daophot style
  float  omin;     // min value
  float  omax;     // max value 
  float  odefault; // default value
  float  value;    // current value

public:

  DaoOption() {}

  //! fill up the option
  DaoOption(const string& Name, const float& Min, 
	    const float&  Max,  const float& Default);
 
  //! use this function to set the option value, it checks against limits
  bool   SetValue(const float& Value);

  //! the current value of the option
  float& Value() { return value;}

  //! the current value of the option
  const float&  Value() const { return value;}

  //! the first two letter of the option name
  string ShortName() const;

  //! check 2 letters against option short name 
  bool operator == (const string& OpName) const { return ShortName() == OpName;}

  //! enable cout << DaoOption << endl;
  friend ostream& operator << (ostream &stream, const DaoOption& opt);
};


//! a vector of DaoOption. Useful in Daophot
class DaophotOptions : public vector<DaoOption> {

  void init_daophot(); // fill with default values

public:

  //! empty constructor fill up with default values
  DaophotOptions();

  //! initialize all options from a ReducedImage
  DaophotOptions(const ReducedImage& Rim);

  //! enable DaophotOptions[ReadoutNoise].Value()
  DaoOption& operator [] (const DaoOptionEnum OptEnum);

  const DaoOption& operator [] (const DaoOptionEnum OptEnum) const;

  //! transforms the vector into an array to be read by daophot fortran routines
  float *GetArray(const int NoptDaophot=20) const;

  //! write on disk daophot style
  bool write(const string FileName="daophot.opt") const;

  //! read the daophot style option file
  bool read(const string FileName="daophot.opt");

  //! enable dumping of all options, ala daophot 
  friend ostream& operator << (ostream &stream, const DaophotOptions& Options);
};


//! a simple routine to write a daophot aperture option file to be read from Daophot aperture photometry
bool WriteDaoAperOpt(vector<double> &Radius, double& InnerSky, 
		     double& OutterSky , const string FileName="photo.opt");


#endif // DAOPHOTOPTION__H

