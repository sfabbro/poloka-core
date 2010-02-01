#ifdef USE_WCS
static bool GuessLinWCS_snifs(const FitsHeader &Head, TanPix2RaDec &Guess)
{
  /* (1000,0) is roughly the location of the ra,dec point in rebinned 
     images containing both CCDs. The WCS should be guessed a little bit more seriously. 
     These images are however useless for astrometry. They should be cut in 
     2 pieces (or patched)
   */
  if (TanLinWCSFromHeader(Head,Guess, /* warn = */ false)) return true;
  /* These acquisition images should be cropped (because both halves are inconsistent)
     In case one only kept the guiding CCD : */
  double crpix1 = Head.HasKey("CRPIX1") ? double(Head.KeyVal("CRPIX1")) : 0.;
  double crpix2 = Head.HasKey("CRPIX2") ? double(Head.KeyVal("CRPIX2")) : 0.;
  // handle here rebinned images :
  double scale = double(Head.KeyVal("SECPIX1"))/0.136;
  return ComputeLinWCS(Head, Point(crpix1+625*scale,crpix2+35*scale), RotationFlip(Down,Left), Guess);
 
}


static int toto_snifs = Add2GuessWCSMap("UHSnifs",GuessLinWCS_snifs);
#endif

#ifdef VIRTUAL_INSTRUMENTS
class SnifsPhoto : public VirtualInstrument { /* TYPE_SELECTOR */
public:  
  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  {if (CheckKey(Head,"INSTRUME","SNIFS")) 
  return new SnifsPhoto;
  return NULL;}

  string TelName() const {return "UH";}
  string InstName() const {return "SnifsPhoto";}
  string TelInstName() const {return "UHSnifs";}

  SIMPLE_TRANSLATOR(TOADPIXS,"SECPIX1");
  RETURN_A_VALUE(TOADGAIN,1);

  //fix TOADUTIM for 2006 post software changes
  FitsKey TOADUTIM(const FitsHeader &Head, const bool Warn) const
  {
    //buggy headers, I only want the last part of the UTC string (time)
    if(Head.HasKey("UTC")){
      string key=Head.KeyVal("UTC");
      return FitsKey("TOADUTIM",key.substr(11));
    }
    //default
    return FitsKey("TOADUTIM",string(Head.KeyVal("UT")));
   }

  //use efftime or opentime for exposure time
  FitsKey TOADEXPO(const FitsHeader &Head, const bool Warn) const
  {
    if(Head.HasKey("EFFTIME")){
      return FitsKey("TOADEXPO",string(Head.KeyVal("EFFTIME")));
    }
    if(Head.HasKey("OPENIME")){
      return FitsKey("TOADEXPO",string(Head.KeyVal("OPENTIME")));
    }
    //default
    return FitsKey("TOADEXPO",string(Head.KeyVal("EXPTIME")));
   }

  //Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  //{
  //  cerr <<" warning not defined for " << TelInstName() << endl;
  //  return Frame;
  //}

};
#endif
