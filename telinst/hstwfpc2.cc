#ifdef VIRTUAL_INSTRUMENTS
class HstWfpc2 : public VirtualInstrument { /* TYPE_SELECTOR */

public :
  
  string TelInstName() const {return "HstWfpc2";}
  string InstName() const { return "Wfpc2";}
  string TelName() const {return "HST";}
  
  static VirtualInstrument* Acceptor(const FitsHeader &Head)
  { if (CheckKey(Head,"INSTRUME","WFPC2")) return new HstWfpc2; return NULL;}
  
  RETURN_A_VALUE(TOADPIXS,0.046);
  
  SIMPLE_TRANSLATOR(TOADEXPO,"EXPTIME");
  RETURN_A_VALUE(TOADRDON,0.72);
  RETURN_A_VALUE(TOADAIRM,0.0);
  SIMPLE_TRANSLATOR(TOADUTIM,"EXPSTART");
  SIMPLE_TRANSLATOR(TOADOBJE,"TARGNAME");
  SIMPLE_TRANSLATOR(TOADGAIN,"ATODGAIN");
  SIMPLE_TRANSLATOR(TOADFILT,"FILTNAM1");  
//  RETURN_A_VALUE(TOADBAND,ToadBand(StringToUpper("FILTNAM1")));	  
  RETURN_A_VALUE(TOADBAND,ToadBand(StringToUpper(Head.KeyVal("FILTNAM1")))); 
  RETURN_A_VALUE(TOADRASC,RaDegToString(Head.KeyVal("RA_TARG")));  
  RETURN_A_VALUE(TOADDECL,DecDegToString(Head.KeyVal("DEC_TARG")));
  SIMPLE_TRANSLATOR(TOADCHIP,"DETECTOR");
  TRANSLATOR_DEC(TOADPZPT)
  {
    
    string filter = Head.KeyVal("FILTNAM1");
    double photflam =  -2.5*log10(double(Head.KeyVal("PHOTFLAM"))) ;
    double zerop = Head.KeyVal("PHOTZPT");
    zerop += photflam;

    // Gives zerop for Cousins RI
    if (filter == "F675W")
	zerop += -0.63;
    if (filter == "F814W")
      zerop += -1.22;

    return FitsKey("TOADPZPT",zerop);
  };

};
#endif
