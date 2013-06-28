// -*- C++ -*- 
// 
// \file skymapper.cc
//
// Open questions:
//  (*) what about TOADPIXS ? 0.5 ? 1. ? 
//  (*) Amp region: amp illu region or amp region + overscans ?
//  (*) c'est quoi ce contructor de Frame avec des XNEWBEG ? Ce sont des mots clefs standard ?
// 
//  (*) RON missing 
//  (*) GAINS missing -- for the moment, I return 1.
//  (*) FILTER missing in the bias images 
//  (*) TOADEQUI ?
//  (*) for the OVERSCAN, we use BIASSEC. Is it correct ? 
//  (*) for the ILLUREGION, we use TRIMSEC. Is it correct ? 
//  (*) how come DATASEC includes the overscan ?
//  (*) UTC or MJD in the headers would help a lot. 
// 
#ifdef VIRTUAL_INSTRUMENTS


class SkyMapper : public VirtualInstrument {  /* TYPE_SELECTOR TYPE_TESTER */
public:
  string TelInstName() const {return "SkyMapper/SkyMapper";}
  string TelName () const {return "SkyMapper";}
  string InstName () const { return "SkyMapper";}
  
  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  {
    // should also check NAXIS1
    if (Head.HasKey("CONHWV") && (int)Head.KeyVal("NAXIS1")<2048)
      return new SkyMapper();
    return 0;
  }
  
  // translators
  FitsKey TOADPIXS(const FitsHeader &Head, const bool Warn) const
  {
    string ccdsum = Head.KeyVal("CCDSUM");
    if (ccdsum == "1 1") return FitsKey("TOADPIXS", 0.4);
    if (ccdsum == "2 2") return FitsKey("TOADPIXS", 0.8);
    cout << "SkyMapper::TOADPIXS : bizarre CCDSUM=" << ccdsum << " dont know what to do " << endl;
    return FitsKey("TOADPIXS", 0.6);
  }
    
  FitsKey TOADRDON(const FitsHeader &Head, const bool Warn) const
  {
    double rdnoise = Head.KeyVal("RDNOISE");
    if(rdnoise<=0) return FitsKey("TOADRDON", 5);// To Be checked
    return FitsKey("TOADRDON",rdnoise);
  }
  
  // for the filter, we have a problem: key not always present 
  SIMPLE_TRANSLATOR(TOADFILT,"FILTNAME");
  SIMPLE_TRANSLATOR(TOADBAND,"FILTNAME");
  RETURN_A_VALUE(TOADEQUI,2000);
  SIMPLE_TRANSLATOR(TOADTYPE,"IMAGETYP");

#ifdef STORAGE
  SIMPLE_TRANSLATOR(TOADUTIM,"UTC-OBS");
  SIMPLE_TRANSLATOR(TOADMJD,"MJD-OBS");
  FitsKey TOADMMJD(const FitsHeader &Head, const bool Warn) const
  { 
    // MJD(2003/01/01) = 52640.
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
  
  // how can I decide that I have one amp only ? 
  // just with the 
  FitsKey TOADNAMP( const FitsHeader &Head, const bool Warn) const
  {
    if( Head.HasKey("AMPLIST") )
      return FitsKey("TOADNAMP", 2);
    return FitsKey("TOADNAMP", 1);
  }

  //  RETURN_A_VALUE(TOADNAMP,2);

  // in principle Illu is correct

  // Add new key used to compute magnitude of the object
  
  double AmpGain(const FitsHeader &Head, const int Iamp) const;
  Frame  PrescanRegion(const FitsHeader &Head, const int Iamp) const;
  Frame  OverscanRegion(const FitsHeader &Head, const int Iamp) const;
  Frame  TotalIlluRegion(const FitsHeader &Head) const;
  Frame  IlluRegion(const FitsHeader &Head, const int Iamp) const;
  Frame  AmpRegion(const FitsHeader &Head, const int Iamp) const;
  bool   IsTrimmed (const FitsHeader &Head) const;
};

bool SkyMapper::IsTrimmed(const FitsHeader &Head) const 
{
  int nx = Head.KeyVal("NAXIS1");
  int ny = Head.KeyVal("NAXIS2");
  return (nx==1024 || nx == 2048);
}

double SkyMapper::AmpGain(const FitsHeader &Head, int Iamp) const
{
  return 1.;
  // string keyname;
  // if (Iamp==1) keyname = "GAINA";
  // else if (Iamp==2) keyname = "GAINB";
  // double gain = Head.KeyVal(keyname, true);
  // return gain;
}

Frame SkyMapper::PrescanRegion(const FitsHeader &Head, const int Iamp) const
{
  if(IsTrimmed(Head)) // if trimmed, no more PrescanRegion
    return Frame();
  int ny = Head.KeyVal("NAXIS2");
  return Frame(0, 0, 49, ny-1);
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
  //  return IlluRegion(Head,1);
  return IlluRegion(Head, 0);
}

Frame SkyMapper::IlluRegion(const FitsHeader &Head, const int Iamp) const
{
  int nx = Head.KeyVal("NAXIS1");
  if (nx>=2048)
    throw(PolokaException("ERROR : SkyMapper::IlluRegion not implemented for full CCDs"));
  
  if (IsTrimmed(Head)) 
    return Frame(Head);
  else
    {
      int ny = Head.KeyVal("NAXIS2");
      Frame illu(50, 0, 1074, ny-1);
      //      string trimsec = (string)Head.KeyVal("TRIMSEC");
      //      fits_imregion_to_frame(trimsec, illu);
      //      fits_imregion_to_frame("[51:1074,1:4096]", illu);
      return illu;
    }
}

Frame SkyMapper::AmpRegion(const FitsHeader &Head, const int Iamp) const
{
  int nx=Head.KeyVal("NAXIS1");
  if (nx>=2048)
    throw(PolokaException("ERROR : SkyMapper::AmpRegion not implemented for full CCDs"));
  return Frame(Head);
}

#endif


