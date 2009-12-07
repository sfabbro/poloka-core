#ifdef USE_WCS
static bool GuessLinWCS_snifs(const FitsHeader &Head, TanPix2RaDec &Guess)
{
  return ComputeLinWCS(Head, Point(1000,0), RotationFlip(Down,Left), Guess);
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
