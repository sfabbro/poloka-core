#ifdef VIRTUAL_INSTRUMENTS
class Andor : public VirtualInstrument { /* TYPE_SELECTOR */

public :

  static VirtualInstrument *Acceptor(const FitsHeader &Head);
  string InstName () const {return "Andor";}

  SIMPLE_TRANSLATOR(TOADEXPO,"EXPOSURE");
  RETURN_A_VALUE(TOADNAMP, 1);

  TRANSLATOR_DEC(TOADPIXS) {
    if (Head.HasKey("CDELT2", false)) {
      double cdelt2 = Head.KeyVal("CDELT2");
      return FitsKey("TOADPIXS",fabs(cdelt2*3600.));
    } else if (Head.HasKey("CD1_1",false)) {
      double cd1_1 = Head.KeyVal("CD1_1");
      double cd1_2 = Head.KeyVal("CD1_2");
      double cd2_2 = Head.KeyVal("CD2_2");
      double cd2_1 = Head.KeyVal("CD2_1");
      return FitsKey("TOADPIXS",sqrt(fabs(cd1_1*cd2_2-cd1_2*cd2_1))*3600.);
    } else 
      return FitsKey("TOADPIXS", 0.);
  }

  TRANSLATOR_DEC(TOADUTIM) {
    double jd = Head.KeyVal("MJD-OBS");    
    string ut = DateFromJulianDay(jd+2400000.5);
    return FitsKey("TOADUTIM", ut);
  }

  // TRANSLATOR_DEC(AIRMASS) {
  // write code to translate airmass from FITS header
  // }

};

class ZadkoAndor : public Andor { 

  string TelInstName() const {return "Zadko/Andor";}
  string TelName () const {return "Zadko";}

  RETURN_A_VALUE(TOADRDON, 6.5);
  RETURN_A_VALUE(TOADGAIN, 2);

};

class TarotAndor : public Andor {
  string TelInstName() const {return "TAROT/Andor";}
  string TelName () const {return "TAROT";}

  RETURN_A_VALUE(TOADRDON, 8.5);
  RETURN_A_VALUE(TOADGAIN, 2);

};

VirtualInstrument* Andor::Acceptor(const FitsHeader &Head) {
  if (DeleteWhiteSpaces(Head.KeyVal("INSTRUME")) == "CAMERAANDORDW436S") {
    if (DeleteWhiteSpaces(Head.KeyVal("TELESCOP")) == "TAROTCHILI")
      return new TarotAndor;
    if (DeleteWhiteSpaces(Head.KeyVal("TELESCOP")) == "ZadkoAustralia")
      return new ZadkoAndor;
  }
  return NULL;
}


#endif
