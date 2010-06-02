/*
  Attempting to add VLT ISAAC as one instrument.

*/
#ifdef VIRTUAL_INSTRUMENTS
class VltHawk : public VirtualInstrument { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "VltHawk";};
  string TelName() const {return "Vlt";}
  string InstName() const { return "Hawk";}
  
  static VirtualInstrument *Acceptor(const FitsHeader &Head) {
    if (CheckKeyToUpper(Head,"INSTRUME","HAWKI")) 
      return new VltHawk;
    return NULL;
  }

  SIMPLE_TRANSLATOR(TOADPIXS,"HIERARCH ESO INS PIXSCALE");
  SIMPLE_TRANSLATOR(TOADUTIM,"UTC");

  FitsKey TOADBAND(const FitsHeader &Head, const bool Warn) const {
    string filt1 = Head.KeyVal("FILTER1", true);
    if (filt1.find("OPEN")==string::npos) return FitsKey("TOADBAND",filt1);
    return FitsKey("TOADBAND",string(Head.KeyVal("FILTER2",true)));
  }


  FitsKey TOADFILT(const FitsHeader &Head, const bool Warn) const {
    return TOADBAND(Head,Warn);
  }

  TRANSLATOR_DEC(TOADGAIN) {

    cout << " Could not find a gain for HAWKI" << endl;
    return  FitsKey("TOADGAIN" , 1);
    /*
    if ( string(Head.KeyVal("HIERARCH ESO DET WIN TYPE")) == "SW")
      return FitsKey("TOADGAIN" , 1);
    if ( string(Head.KeyVal("HIERARCH ESO DET WIN TYPE")) != "LW")
      return FitsKey("TOADGAIN" , 1);
    return  FitsKey("TOADGAIN" , 1);
    */
  }

  TRANSLATOR_DEC(TOADRDON) {
    // found on 
    // http://www.eso.org/observing/etc/doc/ut4/hawki/helphawki.html
    return FitsKey("TOADRDON" , 5); 
  }

  /*
    Convert the RA and DEC from degrees, to hours/degrees, minute and
    seconds.
   */
  /*
    There are two values for the air mass. One from the beginning of
    the observation and one for the end of it. The mean is returned.
   */
  /*
    Observation type.
   */
  SIMPLE_TRANSLATOR(TOADTYPE,"HIERARCH ESO DPR TYPE");

};
#endif
