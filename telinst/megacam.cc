#include "wcsutils.h"

class Megacam: public VirtualInstrument {  /* TYPE_SELECTOR TYPE_TESTER */

public:
  string TelInstName() const {return "CFHT/Megacam";}
  string TelName () const {return "CFHT";}
  string InstName () const { return "Megacam";}
  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { if (StringToUpper(Head.KeyVal("DETECTOR")) == "MEGACAM") return new Megacam; 
  else return NULL;}
  
  bool GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const;

  // translators
  SIMPLE_TRANSLATOR(TOADPIXS,"PIXSCAL1");
  SIMPLE_TRANSLATOR(TOADRDON,"RDNOISE");
  SIMPLE_TRANSLATOR(TOADUTIM,"UTC-OBS");
  SIMPLE_TRANSLATOR(TOADTYPE,"OBSTYPE");
  //SIMPLE_TRANSLATOR(TOADCHIP,"EXTVER");

  TRANSLATOR_DEC(TOADCHIP)
  {
  string chip = Head.KeyVal("EXTVER");
  int chip_value = atoi(chip.c_str());
  if (chip_value<10) chip = "0"+chip;
  return FitsKey("TOADCHIP",chip);
  };

  FitsKey TOADBAND(const FitsHeader &Head, const bool Warn) const
  { string name = Head.KeyVal("FILTER",Warn);
  // extract the first character of .e.g :  'i.MP9701' 
  return FitsKey("TOADBAND",name.substr(0,1));
  }

  RETURN_A_VALUE(TOADNAMP,2);

  // in principle Illu is correct

  // Add new key used to compute magnitude of the object

  double AmpGain(const FitsHeader &Head, const int Iamp) const;
  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const;
  Frame TotalIlluRegion(const FitsHeader &Head) const;
  Frame IlluRegion(const FitsHeader &Head, const int Iamp) const;
  Frame AmpRegion(const FitsHeader &Head, const int Iamp) const;
  bool   IsTrimmed (const FitsHeader &Head) const;
  Frame  extract_frame(const FitsHeader &Head, const char *KeyGenName,int Iamp) const;
};

bool Megacam::IsTrimmed(const FitsHeader &Head) const {
  
  int nx = Head.KeyVal("NAXIS1");
  if(nx>2048)
    return false;
  return true; 
}

bool Megacam::GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const
{  
  if (HasLinWCS(Head)) return TanLinWCSFromHeader(Head,Guess);
  else return false;
}

double Megacam::AmpGain(const FitsHeader &Head, int Iamp) const
{
  string keyname;
  if (Iamp==1) keyname = "GAINA";
  else if (Iamp==2) keyname = "GAINB";
  double gain = Head.KeyVal(keyname, true);
  return gain;
}

Frame Megacam::extract_frame(const FitsHeader &Head, const char *KeyGenName, int Iamp) const
{
  char ampChar = 'A'+(Iamp-1);
  char keyname[16];
  sprintf(keyname,"%s%c",KeyGenName,ampChar);
  string keyval = Head.KeyVal(keyname, true);
  Frame result;
  fits_imregion_to_frame(keyval, result);
  //! Megacam conventions : xMax and YMax are included to the frame
  result.yMax++;
  result.xMax++;
  return result;
}
  

Frame Megacam::OverscanRegion(const FitsHeader &Head, const int Iamp) const
{
  if(IsTrimmed(Head)) // if trimmed no more OverscanRegion
    return Frame();
  return extract_frame(Head,"BSEC",Iamp);
}

Frame Megacam::IlluRegion(const FitsHeader &Head, const int Iamp) const
{
  if(!IsTrimmed(Head))
    return extract_frame(Head,"DSEC",Iamp);
  
  // if trimmed we do something rather manual
  
  // first check if the geometry is what we expect
  string DSECA = Head.KeyVal("DSECA");
  if(DSECA!="[33:1056,1:4612]") {
    cerr << "warning in Megacam::IlluRegion, DSECA is not as expected [33:1056,1:4612]!=" << DSECA<< endl;
  }
  string DSECB = Head.KeyVal("DSECB");
  if(DSECB!="[1057:2080,1:4612]") {
    cerr << "warning in Megacam::IlluRegion, DSECB is not as expected [1057:2080,1:4612]!=" << DSECB << endl;
  }
  Frame result;
  if(Iamp==1)
    fits_imregion_to_frame("[1:1024,1:4612]", result);
  if(Iamp==2)
    fits_imregion_to_frame("[1025:2048,1:4612]", result);
  //! Megacam conventions : xMax and YMax are included to the frame
  result.yMax++;
  result.xMax++;
  return result;
}

Frame Megacam::TotalIlluRegion(const FitsHeader &Head) const
{
  if(IsTrimmed(Head))
    return Frame(Head,WholeSizeFrame);
  
  Frame result;
  if (Head.HasKey("DATASEC")) {
    string datasec = Head.KeyVal("DATASEC");
    if(datasec=="[0:0,0:0]") {
      datasec="[33:2080,1:4612]";
      cout << "Megacam::TotalIlluRegion warning: DATASEC=[0:0,0:0], use instead " << datasec << endl;     
    }
    fits_imregion_to_frame(datasec, result);
  }else{
    result = Frame(Head,WholeSizeFrame);
  }
  return result;
}



Frame Megacam::AmpRegion(const FitsHeader &Head, const int Iamp) const
{
  return extract_frame(Head,"CSEC",Iamp);
}


