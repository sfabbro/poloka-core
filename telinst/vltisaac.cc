/*
  Attempting to add VLT ISAAC as one instrument.

*/
 
class VltIsaac : public VirtualInstrument { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "VltIsaac";};
  string TelName() const {return "Vlt";}
  string InstName() const { return "Isaac";}
  
  static VirtualInstrument *Acceptor(const FitsHeader &Head) {
    if (CheckKeyToUpper(Head,"INSTRUME","ISAAC")) 
      return new VltIsaac;
    return NULL;
  }
  
  SIMPLE_TRANSLATOR(TOADPIXS,"HIERARCH ESO INS PIXSCALE");
  SIMPLE_TRANSLATOR(TOADUTIM,"UTC");

  TRANSLATOR_DEC(TOADFILT) {
    string tmp;
    if ( string(Head.KeyVal("HIERARCH ESO INS FILT1 TYPE")) == "FILTER" )
      tmp =  string(Head.KeyVal("HIERARCH ESO INS FILT1 NAME"));
    if ( string(Head.KeyVal("HIERARCH ESO INS FILT2 FILTER")) == "FILTER" )
      tmp =  tmp + " " + string(Head.KeyVal("HIERARCH ESO INS FILT2 NAME"));
    return FitsKey("TOADFILT", tmp);
  }

  TRANSLATOR_DEC(TOADGAIN) {
    if ( string(Head.KeyVal("HIERARCH ESO DET WIN TYPE")) == "SW")
      //      return FitsKey("TOADGAIN" , 4.5);
      return FitsKey("TOADGAIN" , 1);
    if ( string(Head.KeyVal("HIERARCH ESO DET WIN TYPE")) != "LW")
      //      return FitsKey("TOADGAIN", 7.8);
      return FitsKey("TOADGAIN" , 1);
    return  FitsKey("TOADGAIN" , 1);
  }

  TRANSLATOR_DEC(TOADRDON) {
    if ( string(Head.KeyVal("HIERARCH ESO DET WIN TYPE")) == "SW")
      return FitsKey("TOADRDON" , 11);
    if ( string(Head.KeyVal("HIERARCH ESO DET WIN TYPE")) != "LW")
      return FitsKey("TOADRDON", 40);

    return FitsKey("TOADGAIN" , 0);
  }

  /*
    Convert the RA and DEC from degrees, to hours/degrees, minute and
    seconds.
   */
  TRANSLATOR_DEC(TOADRASC) {
    double raDeg;
    string raStr = Head.KeyVal("RA");

    sscanf(raStr.c_str(),"%lf", &raDeg);
    return FitsKey("TOADRASC", RaDegToString(raDeg));
  }
  TRANSLATOR_DEC(TOADDECL) {
    double decDeg;
    string decStr = Head.KeyVal("DEC");

    sscanf(decStr.c_str(),"%lf", &decDeg);
    return FitsKey("TOADDECL", DecDegToString(decDeg));
  }

  /*
    There are two values for the air mass. One from the beginning of
    the observation and one for the end of it. The mean is returned.
   */
  TRANSLATOR_DEC(TOADAIRM) {
    double startAirm = Head.KeyVal("HIERARCH ESO TEL AIRM START");
    double endAirm = Head.KeyVal("HIERARCH ESO TEL AIRM END");
    return FitsKey("TOADAIRM",(startAirm+endAirm)/2);
  }
  
  /*
    Observation type.
   */
  SIMPLE_TRANSLATOR(TOADTYPE,"HIERARCH ESO DPR TYPE");

};
