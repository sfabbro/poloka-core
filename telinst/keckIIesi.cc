#ifdef USE_WCS
static bool GuessLinWCS_KeckIIesi(const FitsHeader &Head, TanPix2RaDec &Guess)
{
  //double angle =double(Head.KeyVal("ROTDEST"))*M_PI/180.0 + (M_PI/2.0);
  double angle =double(Head.KeyVal("ROTDEST"));
  double xc,yc;
  
  if (angle>45.0) 
    {
      angle = M_PI;
      // The values depend on the illuregion selected.
      xc = 310;
      yc = 975;
    }
  else
    {
      angle = M_PI/2.0;
      xc = 300;
      yc = 650;
    }
  
  Point center(xc,yc);
  double c = cos(angle);
  double s = sin(angle);
  GtransfoLin rotate(0,0,c,-s,s,c);
  
  return ComputeLinWCS (Head, center, rotate, Guess);
}
static int toto_KeckIIesi = Add2GuessWCSMap("KeckIIesi",GuessLinWCS_KeckIIesi);
#endif

#ifdef VIRTUAL_INSTRUMENTS
class KeckIIesi : public VirtualInstrument { /* TYPE_SELECTOR */

public :
  
  string TelInstName() const {return "KeckIIesi";}
  string InstName() const { return "esi";}
  string TelName() const {return "KECKII";}
  
  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { 
    const char *p =StringToUpper(Head.KeyVal("INSTRUME")).c_str(); 
    if  (strstr(p,"ESI")==p) return new KeckIIesi;
    else return NULL;
  }
  
  RETURN_A_VALUE(TOADPIXS,0.15);
  RETURN_A_VALUE(TOADNAMP,1);
  TRANSLATOR_DEC(TOADGAIN)
  {
    if (Head.HasKey("GAIN")) return FitsKey("TOADGAIN",string(Head.KeyVal("GAIN"))); 
    //    else return FitsKey("TOADGAIN",(1.97+2.1)/2.0);
    else return FitsKey("TOADGAIN",(1.97+2.1)/2.0);
  }
  RETURN_A_VALUE(TOADRDON,2.5); 
  SIMPLE_TRANSLATOR(TOADEXPO,"EXPOSURE");
  SIMPLE_TRANSLATOR(TOADTYPE,"OBJECT");
  SIMPLE_TRANSLATOR(TOADFILT,"DWFILNAM");
  SIMPLE_TRANSLATOR(TOADBAND,"DWFILNAM");
  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  {
    string keyval = "[3:12,1:1550]";
    Frame result;
    fits_imregion_to_frame(keyval, result);
    return result;
  }

  Frame TotalIlluRegion(const FitsHeader &Head) const
  {
    string illu ="[201:800,101:1400]";
    // This illu region has an impact on the GuessLinWCS routine

    Frame result;
    fits_imregion_to_frame(illu, result);
    return result;
  }
};
#endif
