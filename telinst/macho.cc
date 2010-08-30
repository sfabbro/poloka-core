#ifdef VIRTUAL_INSTRUMENTS
class Macho : VirtualInstrument { /* TYPE_SELECTOR */

public :
  string TelInstName() const {return "MtStromlo50in/MACHO";}
  string TelName () const {return "MtStromlo50in";}
  string InstName () const {return "MACHO";}

  static VirtualInstrument *Acceptor(const FitsHeader &Head) {
    if (DeleteWhiteSpaces(Head.KeyVal("SITE")) == "Mt.StromloObservatory" && 
	DeleteWhiteSpaces(Head.KeyVal("TELESCOP")) == "50inch")
      return new Macho;
    else 
      return NULL;
  }

  SIMPLE_TRANSLATOR(TOADPIXS,"SCALE");
  SIMPLE_TRANSLATOR(TOADOBJE,"FIELD");
  SIMPLE_TRANSLATOR(TOADRDON,"RD-NOISE");
  SIMPLE_TRANSLATOR(TOADMJD,"MJD-OBS");
  SIMPLE_TRANSLATOR(TOADTYPE,"OBSTYPE");
  
  TRANSLATOR_DEC(TOADCHIP)
  {
    string ampid = Head.KeyVal("AMP-ID");
    return FitsKey("TOADCHIP", atoi(ampid.erase(1,1).c_str()));
  }

  TRANSLATOR_DEC(TOADRASC)
  {
    if (!Head.HasKey("FIELDRA") || !Head.HasKey("FIELDDEC"))
      return FitsKey("TOADRASC", string(Head.KeyVal("RA")));

    string raString  = Head.KeyVal("FIELDRA");
    string decString = Head.KeyVal("FIELDDEC");
    double ra  = RaStringToDeg(raString);
    double dec = DecStringToDeg(decString);
    double cosfactor = cos(dec*M_PI/180) * 60.;
    string ampid = Head.KeyVal("AMP-ID");

    if (ampid == "0_0")
      return FitsKey("TOADRASC", RaDegToString(ra + 20.4/cosfactor));
    else if (ampid == "0_1")
      return FitsKey("TOADRASC", RaDegToString(ra +  7.0/cosfactor));
    else if (ampid == "1_0")
      return FitsKey("TOADRASC", RaDegToString(ra -  9.6/cosfactor));
    else if (ampid == "1_1")
      return FitsKey("TOADRASC", RaDegToString(ra -  9.7/cosfactor));
    else if (ampid == "2_0")
      return FitsKey("TOADRASC", RaDegToString(ra -  9.8/cosfactor));
    else if (ampid == "2_1")
      return FitsKey("TOADRASC", RaDegToString(ra -  9.8/cosfactor));
    else if (ampid == "3_0")
      return FitsKey("TOADRASC", RaDegToString(ra + 11.3/cosfactor));
    else if (ampid == "3_1")
      return FitsKey("TOADRASC", RaDegToString(ra + 11.2/cosfactor));
    else if (ampid == "4_0")
      return FitsKey("TOADRASC", RaDegToString(ra + 17.2/cosfactor));
    else if (ampid == "4_1")
      return FitsKey("TOADRASC", RaDegToString(ra +  7.2/cosfactor));
    else if (ampid == "5_0")
      return FitsKey("TOADRASC", RaDegToString(ra -  5.4/cosfactor));
    else if (ampid == "5_1")
      return FitsKey("TOADRASC", RaDegToString(ra - 15.5/cosfactor));
    else if (ampid == "6_0")
      return FitsKey("TOADRASC", RaDegToString(ra - 15.3/cosfactor));
    else if (ampid == "6_1")
      return FitsKey("TOADRASC", RaDegToString(ra -  5.3/cosfactor));
    else if (ampid == "7_0")
      return FitsKey("TOADRASC", RaDegToString(ra + 11.6/cosfactor));
    else if (ampid == "7_1")
      return FitsKey("TOADRASC", RaDegToString(ra + 11.7/cosfactor));
    else return FitsKey("TOADRASC", raString);
  }
  
  TRANSLATOR_DEC(TOADDECL)
  {
    if (!Head.HasKey("FIELDDEC"))
      return FitsKey("TOADDECL", string(Head.KeyVal("DEC")));
    string decString = Head.KeyVal("FIELDDEC");
    double dec = DecStringToDeg(decString);
    string ampid = Head.KeyVal("AMP-ID");

    if (ampid == "0_0")
      return FitsKey("TOADDECL", DecDegToString(dec + 20.4/60));
    else if (ampid == "0_1")
      return FitsKey("TOADDECL", DecDegToString(dec + 10.4/60));
    else if (ampid == "1_0")
      return FitsKey("TOADDECL", DecDegToString(dec + 15.6/60));
    else if (ampid == "1_1")
      return FitsKey("TOADDECL", DecDegToString(dec +  5.8/60));
    else if (ampid == "2_0")
      return FitsKey("TOADDECL", DecDegToString(dec -  6.5/60));
    else if (ampid == "2_1")
      return FitsKey("TOADDECL", DecDegToString(dec - 16.4/60));
    else if (ampid == "3_0")
      return FitsKey("TOADDECL", DecDegToString(dec - 16.4/60));
    else if (ampid == "3_1")
      return FitsKey("TOADDECL", DecDegToString(dec -  6.5/60));
    else if (ampid == "4_0")
      return FitsKey("TOADDECL", DecDegToString(dec - 10.8/60));
    else if (ampid == "4_1")
      return FitsKey("TOADDECL", DecDegToString(dec - 10.9/60));
    else if (ampid == "5_0")
      return FitsKey("TOADDECL", DecDegToString(dec - 10.9/60));
    else if (ampid == "5_1")
      return FitsKey("TOADDECL", DecDegToString(dec - 10.9/60));
    else if (ampid == "6_0")
      return FitsKey("TOADDECL", DecDegToString(dec + 10.5/60));
    else if (ampid == "6_1")
      return FitsKey("TOADDECL", DecDegToString(dec + 10.4/60));
    else if (ampid == "7_0")
      return FitsKey("TOADDECL", DecDegToString(dec + 16.1/60));
    else if (ampid == "7_1")
      return FitsKey("TOADDECL", DecDegToString(dec +  6.1/60));
    else return FitsKey("TOADDECL", decString);
  }

  FitsKey TOADUTIM(const FitsHeader &Head, const bool warn) const
  {
    double jd = Head.KeyVal("EXPSTART");    
    string ut = DateFromJulianDay(jd+2400000.5);
    return FitsKey("TOADUTIM", ut);
  }
  RETURN_A_VALUE(TOADNAMP, 1);
};
#endif
