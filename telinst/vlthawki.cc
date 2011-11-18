/*
  Attempting to add VLT HAWK-I as one instrument.
*/
#ifdef VIRTUAL_INSTRUMENTS
class VltHawki : public VirtualInstrument { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "VltHawki";};
  string TelName() const {return "Vlt";}
  string InstName() const { return "Hawki";}
  
  static VirtualInstrument *Acceptor(const FitsHeader &Head) {
    if (CheckKeyToUpper(Head,"INSTRUME","HAWKI")) 
      return new VltHawki;
    return NULL;
  }

  SIMPLE_TRANSLATOR(TOADPZPT,"ZP");

  SIMPLE_TRANSLATOR(TOADPIXS,"HIERARCH ESO INS PIXSCALE");
  // SIMPLE_TRANSLATOR(TOADPIXS,"PIXSCALE");
  SIMPLE_TRANSLATOR(TOADUTIM,"UTC");
  // SIMPLE_TRANSLATOR(TOADRDON,"HIERARCH ESO DET CHIP RON");
  // SIMPLE_TRANSLATOR(TOADGAIN,"HIERARCH ESO DET CHIP GAIN");

  TRANSLATOR_DEC(TOADAIRM) {
    float am = 0.5*(float(Head.KeyVal("HIERARCH ESO TEL AIRM START")) +
		    float(Head.KeyVal("HIERARCH ESO TEL AIRM END")));
    return FitsKey("TOADAIRM",am);
  }

  TRANSLATOR_DEC(TOADEXPO) {
    float exptime = int(Head.KeyVal("HIERARCH ESO DET NDIT")) * 
      float(Head.KeyVal("HIERARCH ESO DET DIT"));
    return FitsKey("TOADEXPO",exptime);
  }


  TRANSLATOR_DEC(TOADFILT) {
    string tmp;

    /*
    if ( string(Head.KeyVal("FILTER1")) != "OPEN" )
      tmp =  string(Head.KeyVal("FILTER1"));
    if ( string(Head.KeyVal("FILTER2")) != "OPEN" )
      tmp =  tmp + " " + string(Head.KeyVal("FILTER2"));
    */
    if ( string(Head.KeyVal("HIERARCH ESO INS FILT1 NAME")) != "OPEN" ) {
      tmp =  string(Head.KeyVal("HIERARCH ESO INS FILT1 NAME"));
    }
    if ( string(Head.KeyVal("HIERARCH ESO INS FILT2 NAME")) != "OPEN" ) {
      tmp =  tmp + " " + string(Head.KeyVal("HIERARCH ESO INS FILT2 NAME"));
    }

    return FitsKey("TOADFILT", tmp);
  }

  //  TRANSLATOR_DEC(TOADGAIN) {
  //    if ( string(Head.KeyVal("HIERARCH ESO DET WIN TYPE")) == "SW")
  //      return FitsKey("TOADGAIN" , 4.5);
  //      return FitsKey("TOADGAIN" , 1);
  //    if ( string(Head.KeyVal("HIERARCH ESO DET WIN TYPE")) != "LW")
  //      return FitsKey("TOADGAIN", 7.8);
  //      return FitsKey("TOADGAIN" , 1);
  //    return  FitsKey("TOADGAIN" , 1);
  //  }

//  TRANSLATOR_DEC(TOADRDON) {
//    if ( string(Head.KeyVal("HIERARCH ESO DET WIN TYPE")) == "SW")
//      return FitsKey("TOADRDON" , 11);
//    if ( string(Head.KeyVal("HIERARCH ESO DET WIN TYPE")) != "LW")
//      return FitsKey("TOADRDON", 40);
//
//    return FitsKey("TOADGAIN" , 0);
//  }

  /*
    Convert the RA and DEC from degrees, to hours/degrees, minute and
    seconds.
   */
//  TRANSLATOR_DEC(TOADRASC) {
//    double raDeg;
//    string raStr = Head.KeyVal("RA");
//
//    sscanf(raStr.c_str(),"%lf", &raDeg);
//    return FitsKey("TOADRASC", RaDegToString(raDeg));
//  }
//  TRANSLATOR_DEC(TOADDECL) {
//    double decDeg;
//    string decStr = Head.KeyVal("DEC");
//
//    sscanf(decStr.c_str(),"%lf", &decDeg);
//    return FitsKey("TOADDECL", DecDegToString(decDeg));
//  }

  /*
    There are two values for the air mass. One from the beginning of
    the observation and one for the end of it. The mean is returned.
   */
//  TRANSLATOR_DEC(TOADAIRM) {
//    double startAirm = Head.KeyVal("HIERARCH ESO TEL AIRM START");
//    double endAirm = Head.KeyVal("HIERARCH ESO TEL AIRM END");
//    return FitsKey("TOADAIRM",(startAirm+endAirm)/2);
//  }
  
  /*
    Observation type.
   */
//  SIMPLE_TRANSLATOR(TOADTYPE,"HIERARCH ESO DPR TYPE");

};
#endif
