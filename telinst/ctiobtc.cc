
#ifdef VIRTUAL_INSTRUMENTS
class CtioBtc : public VirtualInstrument { /* TYPE_SELECTOR */
public:
  string TelInstName() const { return "CtioBtc";};
  string InstName() const { return "Btc";};
  string TelName() const { return "CTIO";};


  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { if (CheckKeyToUpper(Head,"INSTRUME","BTC")) return new CtioBtc; return NULL;}

  RETURN_A_VALUE(TOADPIXS,0.425);
  SIMPLE_TRANSLATOR(TOADRDON, "RON");

  TRANSLATOR_DEC(TOADRASC)
  {
      // if LBL changed the RA in the header!!! recognize the fools with the welldept key
      // actually sometimes they didnot, depending on who did the reduction. Usually
      // they did it with image name suffix cln and it seems that they have the key CCDPROC
      // 01/09/00: not always. how to get out of it?? OVERSCAN seems to be ok
      //      if (HasKey("WELLDEPT") && HasKey("OVERSCAN")) return KeyVal("RA");
      int chip = Head.KeyVal("AMPID");
      if ((chip == 2) || (chip == 3) )
	{
	  string raString = Head.KeyVal("RA");
	  double ra = RaStringToDeg(raString);
	  string decString = Head.KeyVal("DEC");
	  double dec = DecStringToDeg(decString);
	  double cosfactor = cos(dec*M_PI/180)*60.;
       	  return FitsKey("TOADRASC",RaDegToString(ra + 20.4/cosfactor));
	}
      else return Head.KeyVal("RA");
    };

  TRANSLATOR_DEC(TOADDECL)
  {
      // same comment as RA
      if (Head.HasKey("WELLDEPT") && Head.HasKey("OVERSCAN")) return Head.KeyVal("DEC");
      int chip = Head.KeyVal("AMPID");
      if ((chip == 3) || (chip == 4))
	{
	  string decString = Head.KeyVal("DEC");
	  double dec = DecStringToDeg(decString);
	  return FitsKey("TOADDECL",DecDegToString(dec - 20.4/60.));
	}
      else  return Head.KeyVal("DEC");
  }
  SIMPLE_TRANSLATOR(TOADEQUI,"EPOCH");
  FitsKey TOADOBJE(const FitsHeader &Head, const bool Warn) const
  {
    string obje =  Head.KeyVal("TITLE");
    RemovePattern(obje," ");
    if (obje.length() == 0) obje = static_cast<string>(Head.KeyVal("OBJECT"));
    return FitsKey("TOADOBJE",obje);
  }
  SIMPLE_TRANSLATOR(TOADCHIP,"AMPID");

};
#endif
