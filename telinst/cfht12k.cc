#ifdef USE_WCS
static bool GuessLinWCS_cfh12k(const FitsHeader &Head, TanPix2RaDec &Guess) {
  if (HasLinWCS(Head)) return TanLinWCSFromHeader(Head,Guess);
  else
    {// untested code : not clear that the offsets are correct
      // not clear either that this code is of any use.
      double ra,dec;
      RaDec2000(Head,ra,dec);
      //double cosdec = cos(dec*M_PI/180.);
      int chip = Head.KeyVal("TOADCHIP");
      double pixs = Head.KeyVal("TOADPIXS");
      double crpix1 = ((chip%6) - 3) * 0.116 * 3600 /pixs;
      double crpix2 = - (chip/6) * 0.225 * 3600/ pixs;
      ComputeLinWCS(Head, Point(crpix1,crpix2), 
		    GtransfoIdentity(), Guess);
      return true;
    }
}
static int toto_cfh12k = Add2GuessWCSMap("Cfht12K",GuessLinWCS_cfh12k);
#endif

/*
  TO BE DONE : implement the necessary swaps due to a bug in the BIASSEC : done in the default decoding
*/
#ifdef VIRTUAL_INSTRUMENTS
class Cfht12K: public VirtualInstrument {  /* TYPE_SELECTOR TYPE_TESTER */

public :
  string TelInstName() const {return "Cfht12K";}
  string TelName () const {return "CFHT";}
  string InstName () const { return "12K";}

  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { if  (StringToUpper(Head.KeyVal("DETECTOR")) == "CFH12K") return new Cfht12K; else return NULL;}

  // translators
  SIMPLE_TRANSLATOR(TOADPIXS,"PIXSCAL1");
  virtual RETURN_A_VALUE(TOADRDON,5);
  SIMPLE_TRANSLATOR(TOADUTIM,"UTC-OBS");
  SIMPLE_TRANSLATOR(TOADTYPE,"OBSTYPE");
  SIMPLE_TRANSLATOR(TOADCHIP,"IMAGEID");

  // in principle Illu and overscan defaults are correct

  // Add new key used to compute magnitude of the object

  TRANSLATOR_DEC(TOADPZPT)
  {
    double ccdqer[12] = {86.,77.,72.,95.,85.,90.,71.,76.,76.,79.,71.,84.};
    double ccdqei[12] = {57.,55.,52.,80.,70.,70.,52.,54.,51.,53.,54.,51.};
    double zpccd9r=26.52;
    double zpccd9i=26.22;
    string band = Head.KeyVal("TOADBAND");
    double zerop=0;
    int chip = Head.KeyVal("IMAGEID");
    double AirMass = Head.KeyVal("AIRMASS");
    if (band=="I") 
      {
	double a=0.05;
	zerop=zpccd9i+2.5*log10(ccdqei[chip]/ccdqei[9])-a*AirMass;
      }
    if (band=="R") 
      {
	double a=0.06;
	zerop=zpccd9r+2.5*log10(ccdqer[chip]/ccdqer[9])-a*AirMass;
      }
    double offset= 0;
    if (Head.HasKey("INIZEROP"))
      offset = Head.KeyVal("INIZEROP");
    else
      offset = 2.5*log10(double(Head.KeyVal("TOADEXPO")));
    return FitsKey("TOADPZPT",zerop+offset);
  };

};
#endif
