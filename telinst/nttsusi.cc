class NttSusi :  public VirtualInstrument { /* TYPE_SELECTOR */

public :

  string TelInstName() const {return "NttSusi";}
  string InstName() const { return "Susi";}
  string TelName() const {return "NTT";}

  static VirtualInstrument* Acceptor(const FitsHeader &Head)
  { if (CheckKey(Head,"TELESCOP","ESO-NTT")) return new NttSusi; return NULL;}

  RETURN_A_VALUE(TOADPIXS,0.161 );
  SIMPLE_TRANSLATOR(TOADEXPO,"EXPTIME");
  RETURN_A_VALUE(TOADRDON,4.55);
  SIMPLE_TRANSLATOR(TOADFILT,"FILT1NAM");  
  
  RETURN_A_VALUE(TOADBAND,ToadBand(StringToUpper(Head.KeyVal("FILT1NAM")))); 
  SIMPLE_TRANSLATOR(TOADGAIN,"OUT1GAIN");
  SIMPLE_TRANSLATOR(TOADAIRM,"AIRMST");
  SIMPLE_TRANSLATOR(TOADTYPE,"DPRCATG");


};
