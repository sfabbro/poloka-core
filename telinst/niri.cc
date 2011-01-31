/*
  Attempting to add VLT ISAAC as one instrument.

*/
#ifdef VIRTUAL_INSTRUMENTS
/* this code does not work (yet) because
I was not able to find a discriminant for NIRI images */


class Niri : public VirtualInstrument { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "Niri/Gemini";};
  string TelName() const {return "Gemini";}
  string InstName() const { return "Niri";}
  
  static VirtualInstrument *Acceptor(const FitsHeader &Head) {
    if (CheckKeyToUpper(Head,"INSTRUME","NIRI")) 
      return new Niri;
    return NULL;
  }

  TRANSLATOR_DEC(TOADPIXS) 
  {
    /* not sure the WCS is always there
    double cd11 = Head.KeyVal("CD1_1");
    double cd12 = Head.KeyVal("CD1_2");
    double cd21 = Head.KeyVal("CD2_1");
    double cd22 = Head.KeyVal("CD2_2");
    return FitsKey("TOADSPIX",sqrt(fabs(cd11*cd22-cd12*cd21))*3600);
    */
    return FitsKey("TOADPIXS",0.116);
  }


  SIMPLE_TRANSLATOR(TOADUTIM,"TIME-OBS");

  FitsKey TOADBAND(const FitsHeader &Head, const bool Warn) const {
    string filt1 = Head.KeyVal("FILTER1", true);
    if (filt1.find("OPEN")==string::npos) return FitsKey("TOADBAND",filt1);
    return FitsKey("TOADBAND",string(Head.KeyVal("FILTER2",true)));
  }


  FitsKey TOADFILT(const FitsHeader &Head, const bool Warn) const {
    return TOADBAND(Head,Warn);
  }

  TRANSLATOR_DEC(TOADGAIN) {
    double ncomb = Head.KeyVal("NCOMBINE");
    // from an email from Chris Lidman : 
    return FitsKey("TOADGAIN",12.3*ncomb);
  }

  TRANSLATOR_DEC(TOADRDON) {
    // found on 
    // http://www.gemini.edu/sciops/instruments/niri/hodappetal.ps 
    // (section 6.3) but I am not sure. I don't think it is used though.
    return FitsKey("TOADRDON" , 45); 
  }

  SIMPLE_TRANSLATOR(TOADTYPE,"OBSTYPE");


  //  SIMPLE_TRANSLATOR(TOADCHIP,"ESO DET CHIP NO");


  


};
#endif
