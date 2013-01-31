// -*- C++ -*- 
// 
// \file des.cc
// telinst for DES
// 
#ifdef VIRTUAL_INSTRUMENTS

#include <assert.h>

class DES : public VirtualInstrument {  /* TYPE_SELECTOR TYPE_TESTER */
  
public:
  string TelInstName() const {return "Blanco/DECam";}
  string TelName () const {return "Blanco";}
  string InstName () const { return "DECam";}
  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { 
    if (StringToUpper(Head.KeyVal("INSTRUME")) == "DECAM") 
      {
	return new DES; 
      }
  else return NULL;
  }
  
  // translators
  SIMPLE_TRANSLATOR(TOADRA,  "TELRA");
  SIMPLE_TRANSLATOR(TOADDEC, "TELDEC");  
  SIMPLE_TRANSLATOR(TOADPIXS,"PIXSCAL1");
  //  SIMPLE_TRANSLATOR(TOADRDON,"RDNOISE");
  SIMPLE_TRANSLATOR(TOADUTIM,"TIME-OBS"); // ALSO DATE-OBS available 
  SIMPLE_TRANSLATOR(TOADTYPE,"OBSTYPE");
  SIMPLE_TRANSLATOR(TOADMJD,"MJD-OBS");

  FitsKey TOADMMJD(const FitsHeader &Head, const bool Warn) const
  { // MJD(2003/01/01) = 52640.
    return FitsKey("TOADMMJD", double(Head.KeyVal("MJD-OBS",Warn))-52640.);
  }

  //SIMPLE_TRANSLATOR(TOADCHIP,"EXTVER");

  TRANSLATOR_DEC(TOADCHIP)
  {
    string chip = Head.KeyVal("CCDNUM");
    int chip_value = atoi(chip.c_str());
    if (chip_value<10) chip = "0"+chip;
    return FitsKey("TOADCHIP",chip);
  }
  
  FitsKey TOADBAND(const FitsHeader &Head, const bool Warn) const
  { 
    string name = Head.KeyVal("FILTER",Warn);
    return FitsKey("TOADBAND",name.substr(0,1));
  }

  RETURN_A_VALUE(TOADNAMP,2);

  // in principle Illu is correct

  // Add new key used to compute magnitude of the object
  
  double AmpGain(const FitsHeader &Head, const int Iamp) const;
  double AmpRdNoise(const FitsHeader &Head, const int Iamp) const;
  FitsKey TOADGAIN(const FitsHeader &Head, const bool Warn) const;
  FitsKey TOADRDON(const FitsHeader &Head, const bool Warn) const;
  Frame  OverscanRegion(const FitsHeader &Head, const int Iamp) const;
  Frame  TotalIlluRegion(const FitsHeader &Head) const;
  Frame  IlluRegion(const FitsHeader &Head, const int Iamp) const;
  Frame  AmpRegion(const FitsHeader &Head, const int Iamp) const;
  bool   IsTrimmed (const FitsHeader &Head) const;
  Frame  extract_frame(const FitsHeader &Head, const char *KeyGenName,int Iamp=0) const;
};

bool DES::IsTrimmed(const FitsHeader &Head) const {
  
  int nx = Head.KeyVal("NAXIS1");
  int ny = Head.KeyVal("NAXIS2");
  if (nx < 2048) // one amplifier
    return (nx==1024);
  return (nx==2048);
}

double DES::AmpGain(const FitsHeader &Head, int Iamp) const
{
  string keyname;
  if (Iamp==1) keyname = "GAINA";
  else if (Iamp==2) keyname = "GAINB";
  double gain = Head.KeyVal(keyname, true);
  return gain;
}

double DES::AmpRdNoise(const FitsHeader &Head, const int Iamp) const
{
  string keyname;
  if(Iamp==1) keyname = "RDNOISEA";
  else if(Iamp==2) keyname = "RDNOISEB";
  double rdnoise = Head.KeyVal(keyname, true);
  return rdnoise;
}

FitsKey DES::TOADGAIN(const FitsHeader &Head, const bool Warn) const
{
  return FitsKey("TOADGAIN", 0.5 * (AmpGain(Head, 1) + AmpGain(Head, 2)));
}

FitsKey DES::TOADRDON(const FitsHeader &Head, const bool Warn) const
{
  return FitsKey("TOADRDON", 0.5 * (AmpRdNoise(Head, 1) + AmpRdNoise(Head, 2)));
}

Frame DES::extract_frame(const FitsHeader &Head, const char *KeyGenName, int Iamp) const
{
  char keyname[16];

  if (Iamp==0)
    {
      sprintf(keyname,"%s",KeyGenName);
    }
  else
    {
      int namp = TOADNAMP(Head, false);
      if(namp == 2)
	{
	  char ampChar = 'A'+(Iamp-1);
	  sprintf(keyname,"%s%c",KeyGenName,ampChar);
	  if (!Head.HasKey(keyname)) sprintf(keyname,"%s",KeyGenName);
	}
      else
	{
	  assert(false);
	  //	  sprintf(keyname, "%s", KeyGenName);
	}
    }
  string keyval = Head.KeyVal(keyname, true);
  Frame result;
  fits_imregion_to_frame(keyval, result);
  
  return result;
}

Frame DES::OverscanRegion(const FitsHeader &Head, const int Iamp) const
{
  if(IsTrimmed(Head)) // if trimmed no more OverscanRegion
    return Frame();
  Frame ret = extract_frame(Head,"BIASSEC",Iamp);
  return ret;
}

Frame DES::IlluRegion(const FitsHeader &Head, const int Iamp) const
{
  return extract_frame(Head,"DATASEC", Iamp);
}

Frame DES::TotalIlluRegion(const FitsHeader &Head) const
{
  if(IsTrimmed(Head))
    return Frame(Head,WholeSizeFrame);
  Frame result;
  
  fits_imregion_to_frame(Head.KeyVal("DATASEC"), result);
  return result;
}

Frame DES::AmpRegion(const FitsHeader &Head, const int Iamp) const
{
  int namp = TOADNAMP(Head, false);
  if(namp == 2)
    return extract_frame(Head,"CCDSEC",Iamp);
  else
    {
      return extract_frame(Head,"CCDSEC",Iamp);
    }
}

#endif

