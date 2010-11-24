#ifdef USE_WCS


#ifdef STORAGE
// Some raw images do not have any CRPIX/CD?_? cards (ex: focus images)
// here are values to cook them up, if needed.

struct WCSTab {
  double crpix1;
  double crpix2;
  double cd1_1;
  double cd1_2;
  double cd2_1;
  double cd2_2;
};


static WCSTab wcscards[] = {
  { -7300.4, 9635.9,  5.194E-05,    0.0,          0.0,         -5.194E-05   },
  { -5220.2, 9647.5,  5.194E-05,    0.0,          0.0,         -5.194E-05   },
  { -3118.8, 9640.5,  5.194E-05,    0.0,          0.0,         -5.194E-05   },
  { -1023.5, 9653.2,  5.185E-05,   -6.20280E-08,  7.31168E-09, -5.17978E-05 },
  {  1083.9, 9639.9,  5.19917E-05, -2.12536E-08, -1.36452E-07, -5.16196E-05 },
  {  3191.2, 9639.9,  5.16979E-05, -7.47275E-08, -5.09956E-08, -5.16072E-05 },
  {  5312.4, 9609.0,  5.194E-05,    0.0,          0.0,         -5.194E-05 },
  {  7394.6, 9587.7,  5.15103E-05, -3.13944E-07,  3.25835E-07, -5.15464E-05 },
  {  9486.6, 9576.0,  5.194E-05,    0.0,          0.0,         -5.194E-05 },
  { -7346.7, 4642.5,  5.194E-05,    0.0,          0.0,         -5.194E-05 },
  { -5257.2, 4649.7,  5.18648E-05,  2.40623E-08, -9.13086E-08, -5.18950E-05 },
  { -3156.8, 4643.3,  5.18075E-05,  2.13197E-08, -5.96327E-08, -5.18790E-05 },
  { -1049.3, 4638.0,  5.19400E-05,  5.21987E-10,  4.09338E-10, -5.19400E-05 },
  {  1066.2, 4633.7,  5.194E-05,    0.0,          0.0,         -5.194E-05 },
  {  3180.8, 4624.1,  5.19627E-05, -1.27804E-07, -1.46688E-07, -5.19433E-05 },
  {  5308.3, 4613.1,  5.194E-05,    0.0,          0.0,         -5.194E-05 },
  {  7393.9, 4610.4,  5.17317E-05, -1.22648E-07, -3.02587E-08, -5.19467E-05 },
  {  9487.0, 4592.2,  5.194E-05,    0.0,          0.0,         -5.194E-05 },
  {  9496.3, 4587.2, -5.194E-05,    0.0,          0.0,          5.194E-05 },
  {  7402.2, 4606.1, -5.16748E-05,  4.50615E-08,  3.14292E-08,  5.18040E-05 },
  {  5295.4, 4613.0, -5.17704E-05,  1.12236E-07, -8.12493E-08,  5.19440E-05 },
  {  3179.9, 4625.3, -5.194E-05,    0.0,          0.0,          5.194E-05 },
  {  1066.2, 4633.7, -5.194E-05,    0.0,          0.0,          5.194E-05 },
  { -1053.7, 4637.7, -5.19321E-05, -1.11346E-08,  2.92162E-08,  5.19464E-05 },
  { -3168.1, 4649.4, -5.17087E-05, -4.59707E-09,  9.82726E-08,  5.18915E-05 },
  { -5233.7, 4632.5, -5.194E-05,    0.0,          0.0,          5.194E-05 },
  { -7348.7, 4644.8, -5.194E-05,    0.0,          0.0,          5.194E-05 },
  {  9498.2, 9574.1, -5.194E-05,    0.0,          0.0,          5.194E-05 },
  {  7405.4, 9596.6, -5.194E-05,    0.0,          0.0,          5.194E-05 },
  {  5303.0, 9612.3, -5.194E-05,    0.0,          0.0,          5.194E-05 },
  {  3194.3, 9631.1, -5.194E-05,    0.0,          0.0,          5.194E-05 },
  {  1091.0, 9634.0, -5.17860E-05,  3.56207E-08, -4.20932E-07,  5.16339E-05 },
  { -1017.7, 9643.7, -5.194E-05,    0.0,          0.0,          5.194E-05 },
  { -3121.3, 9643.4, -5.194E-05,    0.0,          0.0,          5.194E-05 },
  { -5230.7, 9651.0, -5.194E-05,    0.0,          0.0,          5.194E-05 },
  { -7308.7, 9639.7, -5.194E-05,    0.0,          0.0,          5.194E-05 }
};


static bool GuessLinWCS_skymapper(const FitsHeader &Head, TanPix2RaDec &Guess) {
  return false;
}

static int toto_megacam = Add2GuessWCSMap("CFHT/SkyMapper",GuessLinWCS_skymapper);
#endif /* STORAGE */

#endif

#ifdef VIRTUAL_INSTRUMENTS
class SkyMapper: public VirtualInstrument {  /* TYPE_SELECTOR */

public:
  string TelInstName() const {return "SkyMapper/SkyMapper";}
  string TelName () const {return "SkyMapper";}
  string InstName () const { return "SkyMapper";}

  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { 
    if (Head.HasKey("CONTROLR")) return new SkyMapper();
  else return NULL;
  }
  
  
  // translators
  FitsKey TOADPIXS(const FitsHeader &Head, const bool Warn) const
  {
    string ccdsum = Head.KeyVal("CCDSUM");
    if (ccdsum == "1 1") return FitsKey("TOADPIXS",0.5);
    if (ccdsum == "2 2") return FitsKey("TOADPIXS",1);
    cout << " SkyMapper::TOADPIXS : bizarre CCDSUM=" << ccdsum << " dont know what to do " << endl;
}
    
  FitsKey TOADRDON(const FitsHeader &Head, const bool Warn) const
  {
    double rdnoise = Head.KeyVal("RDNOISE");
    if (rdnoise <= 0) return FitsKey("TOADRDON",5);// To Be checked
  }

  SIMPLE_TRANSLATOR(TOADFILT,"FILTNAME");
  SIMPLE_TRANSLATOR(TOADBAND,"FILTNAME");

  RETURN_A_VALUE(TOADEQUI,2000);

  SIMPLE_TRANSLATOR(TOADTYPE,"IMAGETYP");

#ifdef STORAGE
  SIMPLE_TRANSLATOR(TOADUTIM,"UTC-OBS");

  SIMPLE_TRANSLATOR(TOADMJD,"MJD-OBS");

  FitsKey TOADMMJD(const FitsHeader &Head, const bool Warn) const
  { // MJD(2003/01/01) = 52640.
    return FitsKey("TOADMMJD", double(Head.KeyVal("MJD-OBS",Warn))-52640.);
  }
#endif

  //SIMPLE_TRANSLATOR(TOADCHIP,"EXTVER");

  TRANSLATOR_DEC(TOADCHIP)
  {
  string chip = Head.KeyVal("EXTVER");
  int chip_value = atoi(chip.c_str());
  if (chip_value<10) chip = "0"+chip;
  return FitsKey("TOADCHIP",chip);
  };

  FitsKey TOADNAMP( const FitsHeader &Head, const bool Warn) const
  {
    if( Head.HasKey("AMPLIST") )
      return FitsKey("TOADNAMP", 2);
    else
      return FitsKey("TOADNAMP", 1);
  }

  //  RETURN_A_VALUE(TOADNAMP,2);

  // in principle Illu is correct

  // Add new key used to compute magnitude of the object

  double AmpGain(const FitsHeader &Head, const int Iamp) const;
  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const;
  Frame TotalIlluRegion(const FitsHeader &Head) const;
  Frame IlluRegion(const FitsHeader &Head, const int Iamp) const;
  Frame AmpRegion(const FitsHeader &Head, const int Iamp) const;
  bool   IsTrimmed (const FitsHeader &Head) const;
};


bool SkyMapper::IsTrimmed(const FitsHeader &Head) const {
  
  int nx = Head.KeyVal("NAXIS1");
  int ny = Head.KeyVal("NAXIS2");
  return (nx==1024 || nx == 2048);
}

double SkyMapper::AmpGain(const FitsHeader &Head, int Iamp) const
{
  string keyname;
  if (Iamp==1) keyname = "GAINA";
  else if (Iamp==2) keyname = "GAINB";
  double gain = Head.KeyVal(keyname, true);
  return gain;
}


Frame SkyMapper::OverscanRegion(const FitsHeader &Head, const int Iamp) const
{
  if(IsTrimmed(Head)) // if trimmed no more OverscanRegion
    return Frame();
  
  int nx = Head.KeyVal("NAXIS1");
  if (nx>=2048)
    throw(PolokaException(" ERROR : SkyMapper::OverscanRegion not implemented for full CCDs"));

  if (Head.HasKey("BIASSEC")) 
    {
      string bias= Head.KeyVal("BIASSEC");
      Frame biasFrame;
      if (!fits_imregion_to_frame(bias,biasFrame))
	throw(PolokaException("SkyMapper::OverscanRegion cannot decode a frame in "+bias));
      return biasFrame;
    }
  throw(PolokaException("SkyMapper::OverscanRegion finds no BIASSEC in "+Head.FileName()));

}

Frame SkyMapper::TotalIlluRegion(const FitsHeader &Head) const
{
  // to be fixed if we have full CCDs images
  return IlluRegion(Head,1);
}


Frame SkyMapper::IlluRegion(const FitsHeader &Head, const int Iamp) const
{
  int nx = Head.KeyVal("NAXIS1");
  if (nx>=2048)
    throw(PolokaException(" ERROR : SkyMapper::IlluRegion not implemented for full CCDs"));


  if (IsTrimmed(Head)) return Frame(Head);
  else
    {
      // wrong DATASEC (nov 2010)
      Frame illu;
	fits_imregion_to_frame("[51:1074,1:4096]", illu);
      return illu;
    }
}


Frame SkyMapper::AmpRegion(const FitsHeader &Head, const int Iamp) const
{
  int nx=Head.KeyVal("NAXIS1");
  if (nx>=2048)
    throw(PolokaException(" ERROR : SkyMapper::AmpRegion not implemented for full CCDs"));
  return Frame(Head);
}
#endif


