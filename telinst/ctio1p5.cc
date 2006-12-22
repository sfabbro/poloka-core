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
#ifdef USE_WCS
static bool GuessLinWCS_ctio1p5(const FitsHeader &Head, TanPix2RaDec &Guess) {
  return ComputeLinWCS(Head,Head.ImageCenter(),
		       RotationFlip(Down,Left),
		       Guess);
}
static int toto_Ctio1p5mSite2K6 = Add2GuessWCSMap("Ctio1p5mSite2K6",GuessLinWCS_ctio1p5); 
static int toto_Ctio1p5mTek1K2  = Add2GuessWCSMap("Ctio1p5mTek1K2",GuessLinWCS_ctio1p5);
#endif
#ifdef VIRTUAL_INSTRUMENTS
class Ctio1p5 : public VirtualInstrument {
public:
  string TelName() const {return "Ctio1p5";}
  
  SIMPLE_TRANSLATOR(TOADPIXS,"XPIXSIZE");
  SIMPLE_TRANSLATOR(TOADRDON,"RONOISE"); // CHECK THIS !!!
  SIMPLE_TRANSLATOR(TOADFILT,"FILTER2");
  SIMPLE_TRANSLATOR(TOADEQUI,"EPOCH");
};

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
#endif


