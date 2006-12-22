/*
did no check anything (miss images !)

 */
#ifdef VIRTUAL_INSTRUMENTS
class WiynS2kb : public VirtualInstrument {  /* TYPE_SELECTOR */

public :

  string TelInstName() const {return "WiynS2kb";}
  string InstName() const { return "S2kb";};
  string TelName() const { return "WIYN";};

  static VirtualInstrument* Acceptor(const FitsHeader &Head)
  { if (CheckKey(Head,"DETECTOR","s2kb")) return new WiynS2kb; return NULL;}						  
  RETURN_A_VALUE(TOADPIXS,0.197);
  SIMPLE_TRANSLATOR(TOADFILT,"FILTERS");
  TRANSLATOR_DEC(TOADRASC) // see Seb's comments below
  {
    string rastr = Head.KeyVal("RA");
    double radeg = RaStringToDeg(rastr);
    while (radeg > 360) radeg -= 360; // !!!
    return FitsKey("TOADRASC",RaDegToString(radeg));
  }
  RETURN_A_VALUE(TOADBAND,ToadBand(StringToUpper(string(Head.KeyVal("FILTERS")))));
  TRANSLATOR_DEC(TOADEQUI) // EQUINOX is the observation date.. (old fits standard)
  {
    if (CheckKey(Head,"RADECSYS","FK5")) return FitsKey("TOADEQUI",2000.);
    if (CheckKey(Head,"RADECSYS","FK4")) return FitsKey("TOADEQUI",1950.);
    return Head.KeyVal("EQUINOX");
  }
  
};
#endif
