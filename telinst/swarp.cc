#ifdef VIRTUAL_INSTRUMENTS
class Swarp: public VirtualInstrument {  /* TYPE_SELECTOR */

public :
  string TelInstName() const {return "Swarp";}
  string TelName () const {return "no telescope";}
  string InstName () const { return "no instrument";}

  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { if  (string(Head.KeyVal("SOFTNAME")) == "SWarp") return new Swarp; 
  else return NULL;}


  TRANSLATOR_DEC(TOADMJD)
  {
    if (Head.HasKey("MJDDATE",false))
      {
	double mjd = Head.KeyVal("MJDDATE");
	return FitsKey("MJDDATE",mjd);
      }
    if (Head.HasKey("MJD_OBS",false))
      {
	double mjd = Head.KeyVal("MJD_OBS");
	return FitsKey("MJD_OBS",mjd);
      }
    if (Head.HasKey("MJD-OBS",false))
      {
	double mjd = Head.KeyVal("MJD-OBS");
	return FitsKey("MJD-OBS",mjd);
      }
    return VirtualInstrument::TOADMJD(Head,Warn);
  }
    

  
  // translators
  TRANSLATOR_DEC(TOADPIXS)
  {
    if (Head.HasKey("CDELT2", false))
      {
	double cdelt2 = Head.KeyVal("CDELT2");
	return FitsKey("TOADPIXS",fabs(cdelt2*3600.));
      }
    else if (Head.HasKey("CD1_1",false))
      {
	double cd1_1 = Head.KeyVal("CD1_1");
	double cd1_2 = Head.KeyVal("CD1_2");
	double cd2_2 = Head.KeyVal("CD2_2");
	double cd2_1 = Head.KeyVal("CD2_1");
	return FitsKey("TOADPIXS",sqrt(fabs(cd1_1*cd2_2-cd1_2*cd2_1))*3600.);
      }	
    else return FitsKey("TOADPIXS", 0.);
  }

  // Add new key used to compute magnitude of the object

  SIMPLE_TRANSLATOR(TOADPZPT,"ZEROPT");

  TRANSLATOR_DEC(TOADRASC)
  {
    if (Head.HasKey("CRVAL1"))
      {
	double radeg = Head.KeyVal("CRVAL1");
	return FitsKey("TOADRASC",RaDegToString(radeg));
      }
    else return VirtualInstrument::TOADRASC(Head,Warn);
  }


  TRANSLATOR_DEC(TOADDECL)
  {
    if (Head.HasKey("CRVAL2"))
      {
	double radeg = Head.KeyVal("CRVAL2");
	return FitsKey("TOADDECL",DecDegToString(radeg));
      }
    else return VirtualInstrument::TOADDECL(Head,Warn);
  }

};
#endif

