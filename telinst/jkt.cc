class Jkt : public VirtualInstrument {
public:
  static VirtualInstrument *Acceptor(const FitsHeader &Head);// at the end for prototyping
  string TelName() const {return "Jkt";}

  SIMPLE_TRANSLATOR(TOADFILT,"JAGFBAND");
  SIMPLE_TRANSLATOR(TOADEQUI,"EQUINOX");
  SIMPLE_TRANSLATOR(TOADUTIM,"UTSTART");
  SIMPLE_TRANSLATOR(TOADEXPO,"EXPTIME");
  bool GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const;
};

bool Jkt::GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const
{
// Ra,Dec at image center, North = down, East = left
  return ComputeLinWCS(Head,Head.ImageCenter(),
		       RotationFlip(Down,Left),
		       Guess);
}

// it is obviously useless to have 2 derived classes if they do not overload
// anything. But it seems that apart form the bias and ill regions the 2 
// kind of headers were decoded the same way. 
// will become usefull if we ever flatfield this data.

class JktTek1 : public Jkt { /* TYPE_SELECTOR */
public :
  string InstName() const {return "Tek1";}
  string TelInstName() const {return "JktTek1";}

};

class JktSit1 : public Jkt { /* TYPE_SELECTOR */
public :
  string InstName() const {return "Sit1";}
  string TelInstName() const {return "JktSit1";}

  //taken from http://www.ing.iac.es/Engineering/detectors/ultra_site2.htm
  //SITe2 chip
  TRANSLATOR_DEC(TOADPIXS)
    {
      return FitsKey("TOADPIXS",0.33);
    }
  TRANSLATOR_DEC(TOADRDON)
    {
      string keyval = Head.KeyVal("CCDSPEED");
      if (strstr(keyval.c_str(), "FAST")) return FitsKey("TOADRDON",8.); 
      return FitsKey("TOADRDON",7.); 
    }


  //  TRANSLATOR_DEC(TOADGAIN)
  //    {
  //     string keyval = Head.KeyVal("CCDSPEED");
  //     if (strstr(keyval.c_str(), "FAST")) return FitsKey("TOADGAIN",2.58); 
  //     return FitsKey("TOADGAIN",1.95); 
  //   }

  SIMPLE_TRANSLATOR(TOADGAIN,"GAIN")

  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  {
    string oscan;
    oscan = "[2071:2088,1:2120]";    
    Frame result;
    fits_imregion_to_frame(oscan, result);
    return result;
  }
  Frame IlluRegion(const FitsHeader &Head, const int Iamp) const
  {
    return TotalIlluRegion(Head);
  }

  Frame TotalIlluRegion(const FitsHeader &Head) const
  {
    int nx = Head.KeyVal("NAXIS1");
    int ny = Head.KeyVal("NAXIS2");
    string illu;
    if ((nx==1700)&&(ny==1700)) illu ="[1:1700,1:1700]";
      else  illu ="[401:1700,401:1700]";
    // keep only the non vigneted part of the image
    //   illu ="[1:2069,1:2120]";
    Frame result;
    fits_imregion_to_frame(illu, result);
    return result;
  }

 };

class JktSit2 : public Jkt { /* TYPE_SELECTOR */
public :
  string InstName() const {return "Sit2";}
  string TelInstName() const {return "JktSit2";}

  //taken from http://www.ing.iac.es/Engineering/detectors/ultra_site2.htm
  //SITe2 chip
  TRANSLATOR_DEC(TOADPIXS)
    {
      return FitsKey("TOADPIXS",0.33);
    }
  TRANSLATOR_DEC(TOADRDON)
    {
      string keyval = Head.KeyVal("CCDSPEED");
      if (strstr(keyval.c_str(), "FAST")) return FitsKey("TOADRDON",7.); 
      return FitsKey("TOADRDON",6.); 
    }


  //  TRANSLATOR_DEC(TOADGAIN)
  //    {
  //     string keyval = Head.KeyVal("CCDSPEED");
  //     if (strstr(keyval.c_str(), "FAST")) return FitsKey("TOADGAIN",2.58); 
  //     return FitsKey("TOADGAIN",1.95); 
  //   }

  SIMPLE_TRANSLATOR(TOADGAIN,"GAIN")

  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  {
    string oscan;
    oscan = "[2071:2088,1:2120]";    
    Frame result;
    fits_imregion_to_frame(oscan, result);
    return result;
  }
  Frame IlluRegion(const FitsHeader &Head, const int Iamp) const
  {
    return TotalIlluRegion(Head);
  }

  Frame TotalIlluRegion(const FitsHeader &Head) const
  {
    int nx = Head.KeyVal("NAXIS1");
    int ny = Head.KeyVal("NAXIS2");
    string illu;
    if ((nx==1700)&&(ny==1700)) illu ="[1:1700,1:1700]";
      else  illu ="[401:1700,401:1700]";
    // keep only the non vigneted part of the image
    //   illu ="[1:2069,1:2120]";
    Frame result;
    fits_imregion_to_frame(illu, result);
    return result;
  }

 };

VirtualInstrument *Jkt::Acceptor(const FitsHeader &Head)
{ 
  if (CheckKey(Head,"DETECTOR","SIT1")) return new JktSit1; 
  if (CheckKey(Head,"DETECTOR","SIT2")) return new JktSit2; 
  if (CheckKey(Head,"DETECTOR","TEK1")) return new JktTek1;
  return NULL;
}
