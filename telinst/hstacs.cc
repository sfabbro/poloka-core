
class HstACS : public VirtualInstrument { /* TYPE_SELECTOR */

public :
  
  string TelInstName() const {return "HstACS";}
  string InstName() const { return "ACS";}
  string TelName() const {return "HST";}
  
  static VirtualInstrument* Acceptor(const FitsHeader &Head)
  { 
    if (CheckKey(Head,"INSTRUME","ACS")) return new HstACS; return NULL;}
  
  RETURN_A_VALUE(TOADPIXS,0.05);
  
  SIMPLE_TRANSLATOR(TOADEXPO,"EXPTIME");
  // I set the readout noise to be the one for A amplifier (That could be B,C or D).
  // I don't know yet how to determine the amplifier
  SIMPLE_TRANSLATOR(TOADRDON,"READNSEA");
  RETURN_A_VALUE(TOADAIRM,0.0);
  SIMPLE_TRANSLATOR(TOADUTIM,"EXPSTART");
  SIMPLE_TRANSLATOR(TOADOBJE,"TARGNAME");
  // There are 4 different calibrated gain for each amplifier
  // We should use it but I don't yet how to determine the amplifier
  SIMPLE_TRANSLATOR(TOADGAIN,"CCDGAIN");
  
  
//  RETURN_A_VALUE(TOADBAND,ToadBand(StringToUpper("FILTNAM1")));	  
  RETURN_A_VALUE(TOADBAND,ToadBand(StringToUpper(Head.KeyVal("TOADFILT")))); 
  RETURN_A_VALUE(TOADRASC,RaDegToString(Head.KeyVal("RA_TARG")));  
  RETURN_A_VALUE(TOADDECL,DecDegToString(Head.KeyVal("DEC_TARG")));
  SIMPLE_TRANSLATOR(TOADCHIP,"CCDCHIP");
  TRANSLATOR_DEC(TOADPZPT)
  {
    
    string filter = Head.KeyVal("TOADFILT");
    double photflam =  -2.5*log10(double(Head.KeyVal("PHOTFLAM"))) ;
    double zerop = Head.KeyVal("PHOTZPT");
    zerop += photflam;

    // Gives zerop for Cousins RI

    return FitsKey("TOADPZPT",zerop);
  };

  TRANSLATOR_DEC(TOADFILT)
  {
    string filter = Head.KeyVal("FILTER1");
    if (strstr(filter.c_str(),"CLEAR"))
      filter = Head.KeyVal("FILTER2");
    
    return FitsKey("TOADFILT",filter);
  }  

};
