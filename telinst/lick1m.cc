// This stuff is only valid for Lick 1m reduced image coming 
// from the LBL database

#include "fitsimage.h" // for KeyVal()

class Lick1m : public VirtualInstrument {
public:
  string TelName() const {return "Lick1m";}

  SIMPLE_TRANSLATOR(TOADEXPO,"EXPOSURE");
  SIMPLE_TRANSLATOR(TOADUTIM,"TIME");
  SIMPLE_TRANSLATOR(TOADRDON,"RONOISE");
  RETURN_A_VALUE(TOADTYPE,"Unknown"); 
  bool GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const;
};

bool Lick1m::GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const
{
// Ra,Dec at image center, North = down, East = left
  return ComputeLinWCS(Head,Head.ImageCenter(),
		       RotationFlip(Down,Left),
		       Guess);
}

class Lick1mDewar5 : public Lick1m { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "Lick1mDewar5";};
  string InstName() const { return "Dewar5";}
  SIMPLE_TRANSLATOR(TOADEXPO,"EXPOSURE");
  static VirtualInstrument *Acceptor(const FitsHeader &Head)
    { 
     if (CheckKey(Head,"DET_SPEC", "li1m_dewar5")) 
       return new Lick1mDewar5;
     return NULL;
    }
};

class Lick1mDewar5_large : public Lick1m { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "Lick1mDewar5_large";};
  string InstName() const { return "Dewar5_large";}

  static VirtualInstrument *Acceptor(const FitsHeader &Head)
    { 
      double ro = (double) Head.KeyVal("RONOISE");
     if (CheckKey(Head,"DET_SPEC", "li1m_dewar5") && ro > 1.6 && ro < 1.8) 
       return new Lick1mDewar5_large;
     return NULL;
    }
};

class Lick1mDewar2 : public Lick1m { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "Lick1mDewar2";};
  string InstName() const { return "Dewar2";}

  static VirtualInstrument *Acceptor(const FitsHeader &Head)
    { 
     if (CheckKey(Head,"DET_SPEC", "li1m_dewar2") || CheckKey(Head,"DET_SPEC", "li1m_dewar2_small")) 
       return new Lick1mDewar2;
     return NULL;
    }
};


#ifdef OLD_FITSTOAD
FitsKey FitsHeader::Lick1mDewar2Format(const string& KeyName) const 
{
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.18); // check this
  if (KeyName == "TOADFILT") return KeyVal("FILTER");
  if (KeyName == "TOADEXPO") return KeyVal("EXPOSURE");
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADRDON") return KeyVal("RONOISE");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADEQUI") return KeyVal("EQUINOX") ;
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("TIME") ;
  if (KeyName == "TOADDATE") {
    int dd, mm, yy ;
    string kval = KeyVal("DATE-OBS") ;
    char buf[64] ;
    if(sscanf(kval.c_str(), "%d-%d-%dT%s", &yy, &mm, &dd, buf)) {
      sprintf(buf, "%02d/%02d/%04d", dd, mm, yy) ;
      return FitsKey(KeyName, string(buf)) ;
    }
  }
  if (KeyName == "TOADSCAN")
    return FitsKey(KeyName,string("[0,0;0,0]")) ;
  if (KeyName == "TOADILLU")
    return FitsKey(KeyName,string("[0,1024;0,1024]")) ;
  if (KeyName == "TOADTYPE") return FitsKey(KeyName, "UNKNOWN");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return FitsKey(KeyName,0); // begining @ 0
  if (KeyName == "TOADBAND") {
    string kval = KeyVal("FILTER") ;
    return FitsKey(KeyName, ToadBand(StringToUpper(kval)));
  }
  if (KeyName == "TOADNAMP") return FitsKey(KeyName, 1) ;
  return FitsKey(KeyName,NOVAL);
}
#endif
