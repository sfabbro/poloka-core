// 
// Site2K_6 CCD mounted on the CTIO 1.5m telescope. 
// Intensively used to follow-up the Spring '99 nearby SNIa's
// 2 remarks: 
//      (*) this is a 4 amplifier instrument, with 4 different gains. 
//          What we use here, is the mean, computed when these images
//          were loaded into the database. If we want this to work with 
//          raw images, this should be recomputed on the fly from the GAIN[1-4]
//          keys. Same pb with the RDNOISE
//      (*) this camera has 2 filter wheels. In practice only the 
//          2d wheel has been used ==> we use the FILTER2 key. 
//          This is not robust at all. 
// 

class Ctio1p5 : public VirtualInstrument {
public:
  string TelName() const {return "Ctio1p5";}
  
  SIMPLE_TRANSLATOR(TOADPIXS,"XPIXSIZE");
  SIMPLE_TRANSLATOR(TOADRDON,"RONOISE"); // CHECK THIS !!!
  SIMPLE_TRANSLATOR(TOADFILT,"FILTER2");
  SIMPLE_TRANSLATOR(TOADEQUI,"EPOCH");
  bool GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const;
  
};

bool Ctio1p5::GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const
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

class Ctio1p5mSite2K6 : public Ctio1p5 { /* TYPE_SELECTOR */
public :
  string InstName() const {return "Site2K6";}
  string TelInstName() const {return "Ctio1p5mSite2K6";}

  static VirtualInstrument *Acceptor(const FitsHeader &Head)
    { if (CheckKey(Head,"DETECTOR","Site2K_6")) return new Ctio1p5mSite2K6; return NULL;}
    
};
    
class Ctio1p5mTek1K2 : public Ctio1p5 { /* TYPE_SELECTOR */
 public:
  string InstName() const {return "Tek1K2";}
  string TelInstName() const {return "Ctio1p5mTek1K2";}

  static VirtualInstrument *Acceptor(const FitsHeader &Head)
    { if (CheckKey(Head,"DETECTOR","Tek1K_2")) return new Ctio1p5mSite2K6; return NULL;}

};



#ifdef OLD_TOADS
FitsKey FitsHeader::Ctio1p5mSite2K6Format(const string &KeyName) const
{
  if (KeyName == "TOADINST") return KeyVal("DETECTOR");
  if (KeyName == "TOADPIXS") return KeyVal("XPIXSIZE");
  if (KeyName == "TOADFILT") { 
      string kval = KeyVal("FILTER2");
      kval = StringToUpper(kval) ;
      return FitsKey(KeyName, kval) ;
  }
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADRASC") return KeyVal("RA") ;
  if (KeyName == "TOADDECL") return KeyVal("DEC") ;
  if (KeyName == "TOADRDON") return KeyVal("RONOISE"); // check this!
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADEQUI") return KeyVal("EPOCH") ;
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UT") ;
  if (KeyName == "TOADDATE") {
      // from dd/mm/yy to dd-mm-yyyy
      char buff[64] ;
      int dd, mm, yy ;
      string kval = KeyVal("DATE-OBS") ;
      if( sscanf(kval.c_str(), "%d/%d/%d",&dd,&mm,&yy) ) {
	  if(yy<50) yy += 2000 ; else yy += 1900 ; // not very convincing, but...
	  sprintf(buff, "%02d/%02d/%04d", dd, mm, yy) ;
	  return FitsKey(KeyName, string(buff)) ;
      }
  }
  // I don't know exactly what to put here. The CTIO images
  // I am currently using have been clipped during the reduction
  // process. -- I hope putting [0,0;0,0] won't make toads crash. --nrl
  if (KeyName == "TOADSCAN") return FitsKey(KeyName,string("[0,0;0,0]")) ;
  if (KeyName == "TOADILLU") return FitsKey(KeyName,string("[0,2048;0;2048]")) ;
  if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return FitsKey(KeyName,0); // begining @ 0
  if (KeyName == "TOADBAND") {
      string kval = KeyVal("FILTER2");
      kval = StringToUpper(kval) ;
      return FitsKey(KeyName, kval) ;
  } // same as TOADSFILT...
  if (KeyName == "TOADNAMP") {
      int nax, nay, namp ;
      string kval = KeyVal("NAMPSXY") ;
      sscanf(kval.c_str(), "'%d%d'", &nax, &nay) ;
      namp = nax * nay ; if(namp <=0) namp = 1 ;
      return FitsKey(KeyName, namp) ;
  }
  return FitsKey(KeyName,NOVAL);      
}


FitsKey FitsHeader::Ctio1p5mTek1K2Format(const string &KeyName) const
{
  if (KeyName == "TOADINST") return KeyVal("DETECTOR");
  if (KeyName == "TOADPIXS") return KeyVal("XPIXSIZE");
  if (KeyName == "TOADFILT") { 
      string kval = KeyVal("FILTER2");
      kval = StringToUpper(kval) ;
      return FitsKey(KeyName, kval) ;
  }
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADRASC") return KeyVal("RA") ;
  if (KeyName == "TOADDECL") return KeyVal("DEC") ;
  if (KeyName == "TOADRDON") return KeyVal("RONOISE"); // check this!
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADEQUI") return KeyVal("EPOCH") ;
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UT") ;
  if (KeyName == "TOADDATE") {
      // from dd/mm/yy to dd-mm-yyyy
      char buff[64] ;
      int dd, mm, yyyy ;
      string kval = KeyVal("DATE-OBS") ;
      if( sscanf(kval.c_str(), "%d-%d-%dT%s",&yyyy,&mm,&dd, buff) ) {
	  sprintf(buff, "%02d/%02d/%04d", dd, mm, yyyy) ;
	  return FitsKey(KeyName, string(buff)) ;
      }
  }
  // I don't know exactly what to put here. The CTIO images
  // I am currently using have been clipped during the reduction
  // process. -- I hope putting [0,0;0,0] won't make toads crash. --nrl
  if (KeyName == "TOADSCAN") return FitsKey(KeyName,string("[0,0;0,0]")) ;
  if (KeyName == "TOADILLU") return FitsKey(KeyName,string("[0,1024;0;1024]")) ;
  if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return FitsKey(KeyName,0); // begining @ 0
  if (KeyName == "TOADBAND") {
      string kval = KeyVal("FILTER2");
      kval = StringToUpper(kval) ;
      return FitsKey(KeyName, kval) ;
  } // same as TOADSFILT...
  if (KeyName == "TOADNAMP") {
      int nax, nay, namp ;
      string kval = KeyVal("NAMPSXY") ;
      sscanf(kval.c_str(), "'%d%d'", &nax, &nay) ;
      namp = nax * nay ; if(namp <=0) namp = 1 ;
      return FitsKey(KeyName, namp) ;
  }
  return FitsKey(KeyName,NOVAL);      
}

#endif
