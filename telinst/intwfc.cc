
/* 3 leaf classes : 
   Int1Wfc (1 chip)
   Int4Wfc (4chips)
   IntWfcNewDaq (4 chips after daq upgrade in 99)
*/

/* CHECK : 1<=(TOADCHIP,IMAGEID)<=4 */
/* Illu overscan and GuessLinWcs are implemented for the 4-chip 
setup as it was in 2000.
*/
#ifdef USE_WCS
static bool GuessLinWCS_int1wfc(const FitsHeader &Head, TanPix2RaDec &Guess) {
  cerr << " GuessWCS still to be impemented for Int1Wfc " << endl;
  return false;
}
static int toto_Int1Wfc = Add2GuessWCSMap("Int1Wfc",GuessLinWCS_int1wfc);

static bool GuessLinWCS_int4wfc(const FitsHeader &Head, TanPix2RaDec &Guess) {
  int id = Head.KeyVal("TOADCHIP");
  switch (id)
    {
    case 1: ComputeLinWCS(Head,Point(-1120 ,1980.), RotationFlip(Left,Down), Guess); break;
    case 2: ComputeLinWCS(Head,Point( 4200 , 980.), RotationFlip(Down,Right),Guess); break;
    case 3: ComputeLinWCS(Head,Point( 3054.,1978.), RotationFlip(Left,Down), Guess); break;
    case 4: ComputeLinWCS(Head,Point(  993, 1995.), RotationFlip(Left,Down), Guess); break;
    default: ComputeLinWCS(Head,Head.ImageCenter(),GtransfoIdentity(), Guess) ; break;
    }
  return true;
}
static int toto_Int4Wfc      = Add2GuessWCSMap("Int4Wfc",GuessLinWCS_int4wfc);
static int toto_IntWfcNewDaq = Add2GuessWCSMap("IntWfcNewDaq",GuessLinWCS_int4wfc);
#endif

#ifdef VIRTUAL_INSTRUMENTS
class IntWfc : public VirtualInstrument {    /* TYPE_SELECTOR */
public :
  string TelName() const {return "INT";};
  static VirtualInstrument *Acceptor(const FitsHeader &Head); // at the end


  SIMPLE_TRANSLATOR(TOADPIXS,"SECPPIX")
  TRANSLATOR_DEC(TOADFILT)
    {return FitsKey("TOADFILT", string(Head.KeyVal("WFFBAND"))+string(Head.KeyVal("WFFPSYS")));}
  virtual SIMPLE_TRANSLATOR(TOADRDON,"READNOIS") // overloaded for IntWfcNewDaq
  SIMPLE_TRANSLATOR(TOADRASC,"CAT-RA")
  SIMPLE_TRANSLATOR(TOADDECL,"CAT-DEC")
  TRANSLATOR_DEC(TOADEQUI)
    {
      string keyval = Head.KeyVal("CAT-EQUI");
      if (strstr(keyval.c_str(), "B1950")) return FitsKey("TOADEQUI",1950.0); 
      return FitsKey("TOADEQUI",2000.0); 
    }
  SIMPLE_TRANSLATOR(TOADUTIM,"UTSTART");
  TRANSLATOR_DEC(TOADBAND)
    { 
      string filter = Head.KeyVal("WFFBAND");
      return FitsKey("TOADBAND",ToadBand(StringToUpper(filter)));	  
    }
#ifdef INSTRUMENTWCS
  bool GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const
  {  cerr << " GuessWCS still to be impemented for IntWfc " << endl;return false; }
#endif
};
class Int4Wfc : public IntWfc  /* TYPE_TESTER */
{
  string TelInstName() const { return "Int4Wfc";}
  string InstName() const { return "Wfc(4chips)";}


  virtual FitsKey TOADCHIP(const FitsHeader &Head, const bool Warn) const
  {
    int chip = Head.KeyVal("DASCHAN", Warn);
    if (chip == 5) chip = 2;
    return FitsKey("TOADCHIP",chip);
  };
};

class Int1Wfc : public IntWfc
{
  string TelInstName() const {return "Int1Wfc";}
  string InstName () const { return "Wfc1chip";}
};  


class IntWfcNewDaq: public Int4Wfc {      /*  TYPE_TESTER */
public :
  string TelInstName() const {return "IntWfcNewDaq";}
  string InstName () const { return "WfcNewDaq";}
  SIMPLE_TRANSLATOR(TOADCHIP,"IMAGEID");

  FitsKey TOADRDON(const FitsHeader &Head, const bool Warn) const
  {
    double noise;
    if (Head.HasKey("RDNOISE") && ( (noise = Head.KeyVal("RDNOISE")) > 0. ))
	{
	    return FitsKey("TOADRDON",noise);
	}
    else 
      { 
	static double ronoise[4] = {11,13,9,11}; 
	return FitsKey("TOADRDON",ronoise[int(Head.KeyVal("IMAGEID"))-1]);
      }
    }

  Frame TotalIlluRegion(const FitsHeader &Head) const
  {Frame result; 
  fits_imregion_to_frame(string(Head.KeyVal("TRIMSEC")), result);
  return result;
  };

  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  {
    string oscan;
    oscan = "[2110:2150,1:2048]";    
    Frame result;
    fits_imregion_to_frame(oscan, result);
    return result;
  }
  

};



VirtualInstrument* IntWfc::Acceptor(const FitsHeader &Head)
{
  string KEYINST = StringToUpper(string(Head.KeyVal("INSTRUME")));
  if (KEYINST != "WFC") return NULL;
  double jd = Head.KeyVal("JD");
  if (Head.HasKey("JD") &&  jd > 2451410.5822505) /* roughly Aug 20th 1999 when the Daq was changed */
    { return new(IntWfcNewDaq);}
  if (Head.HasKey("CCDTYPE"))
    {
      string keyccd  = Head.KeyVal("CCDTYPE");
      keyccd = StringToLower(keyccd);
      if (keyccd=="loral") {return new Int1Wfc;}
      if (keyccd=="eev42") {return new Int4Wfc;}
    }
  // WARNING ! THIS IS THE ONLY WAY WE FOUND TO SEPARATE 
  // Int1Wfc from Int4Wfc using main header keywords..
  else if (Head.HasKey("JD"))
    {
      if (jd > 2451000.0) return new Int4Wfc; 
      else return new Int1Wfc; 
    }

return NULL;
}
#endif
