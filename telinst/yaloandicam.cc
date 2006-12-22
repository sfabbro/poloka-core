#ifdef VIRTUAL_INSTRUMENTS
static double jd_ref1 = JulianDay(1, 11, 1999); // change in gain and readout
static double jd_ref2 = JulianDay(8,  9, 2000); // a major date of increase in YALO DAQ

class YaloAndicam : public VirtualInstrument { /* TYPE_SELECTOR */
public :
  
  static VirtualInstrument *Acceptor(const FitsHeader &Head); // at the end for prototyping
  string TelName() const { return "Yalo";}
  SIMPLE_TRANSLATOR(TOADAIRM,"SECZ");
  bool IsBinned(const FitsHeader &Head) const 
  { 
    int binx = Head.KeyVal("CCDXBIN");
    int biny = Head.KeyVal("CCDYBIN");
    return (( binx== 2) && (biny == 2));
  }
  TRANSLATOR_DEC(TOADPIXS)
  {
    if (!IsBinned(Head)) return FitsKey("TOADPIXS",0.298);
    return FitsKey("TOADPIXS",0.596);
  };
};

// YALO ANDICAM setup before 01/11/1999
class YaloAndicam1999 : public YaloAndicam 
{
public:
  string InstName() const { return "Andicam1999";}
  string TelInstName() const { return "YaloAndicam1999";}
  RETURN_A_VALUE(TOADGAIN,4.0);
  RETURN_A_VALUE(TOADRDON,12.0);
  RETURN_A_VALUE(TOADNAMP,2);
  SIMPLE_TRANSLATOR(TOADFILT,"FILTERID");
  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  {
    string oscan;
    if (!IsBinned(Head)) oscan = "[2066:2079,1:2048]";
    else oscan = "[1033:1040,1:1024]";
    
    Frame result;
    fits_imregion_to_frame(oscan, result);
    return result;
  }
  
  Frame IlluRegion(const FitsHeader &Head, const int Iamp) const
  {      
    return TotalIlluRegion(Head);
  }
  
  Frame TotalIlluRegion(const FitsHeader &Head) const
  {
    string illu;
    if (!IsBinned(Head)) illu ="[1040:2063,1:2048]";
    else illu ="[520:1031,1:1024]";
    Frame result;
    fits_imregion_to_frame(illu, result);
    return result;
  }
};

// YALO ANDICAM setup before 08/09/2000
class YaloAndicam2000 : public YaloAndicam1999
{
  string InstName() const { return "Andicam2000";}
  string TelInstName() const { return "YaloAndicam2000";}
  RETURN_A_VALUE(TOADGAIN,8.0);
  RETURN_A_VALUE(TOADRDON,20.0);
  SIMPLE_TRANSLATOR(TOADFILT,"CCDFLTID");
};

// YALO ANDICAM setup after 08/09/2000
class YaloAndicam2001 : public YaloAndicam
{
  string InstName() const { return "Andicam2001";}
  string TelInstName() const { return "YaloAndicam2001";}
  RETURN_A_VALUE(TOADGAIN,3.6);
  RETURN_A_VALUE(TOADRDON,12.0);
  RETURN_A_VALUE(TOADNAMP,2);
  SIMPLE_TRANSLATOR(TOADFILT,"CCDFLTID");
  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  {
    string oscan;
    if (Iamp==1) oscan = "[8:32,1:2048]";
    else oscan = "[2113:2142,1:2048]";
    Frame result;
    fits_imregion_to_frame(oscan, result);
    return result;
  }
  
  Frame IlluRegion(const FitsHeader &Head, const int Iamp) const
  {      
    string illu;
    if (Iamp == 1) illu ="[49:1072,1:2048]";
    else illu = "[1073:2094,1:2048]";
    Frame result;
    fits_imregion_to_frame(illu, result);
    return result;
  }
  
  Frame TotalIlluRegion(const FitsHeader &Head) const
  {
    string illu = "[49:2094,1:2048]";
    Frame result;
    fits_imregion_to_frame(illu, result);
    return result;
  }
};

// YALO ANDICAM infrared array
class YaloAndicamIR2001 : public YaloAndicam
{
  string TelInstName() const {return "YaloAndicamIR";}
  string InstName() const { return "AndicamIR";}
  RETURN_A_VALUE(TOADNAMP,2);
  RETURN_A_VALUE(TOADPIXS,0.222);
  RETURN_A_VALUE(TOADGAIN,6.6);
  RETURN_A_VALUE(TOADRDON,14.1);
  SIMPLE_TRANSLATOR(TOADFILT,"IRFLTID");
  
  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  {
    return Frame(0,0,0,0);
  }
  
  Frame IlluRegion(const FitsHeader &Head, const int Iamp) const
  {      
    return TotalIlluRegion(Head);
  }
  
  Frame TotalIlluRegion(const FitsHeader &Head) const
  {
    string illu ="[1:1024,1:1024]";
    Frame result;
    fits_imregion_to_frame(illu, result);
    return result;
  }
};  

VirtualInstrument *YaloAndicam::Acceptor(const FitsHeader &Head)
{ 
  string keyinst = StringToUpper(Head.KeyVal("INSTRUME"));
  RemovePattern(keyinst," ");
  if (keyinst == "ANDYCAMIR") return new YaloAndicamIR2001; // probably add a new DAQ before jd_ref2  
  if (keyinst == "ANDICAM")  // watch out the Y and the I on Andi!!
    {
      double jd = Head.KeyVal("JD"); // do not use JulianDay function because it uses TOADKEYS.
      if (jd < jd_ref1)  return new YaloAndicam1999;
      if (jd < jd_ref2)  return new YaloAndicam2000;
      return new YaloAndicam2001; 
    }
  return NULL;
}
#endif
