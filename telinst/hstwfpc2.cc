
class HstWfpc2 : public VirtualInstrument { /* TYPE_SELECTOR */

public :
  
  string TelInstName() const {return "HstWfpc2";}
  string InstName() const { return "Wfpc2";}
  string TelName() const {return "HST";}
  
  static VirtualInstrument* Acceptor(const FitsHeader &Head)
  { if (CheckKey(Head,"INSTRUME","WFPC2")) return new HstWfpc2; return NULL;}
  
  RETURN_A_VALUE(TOADPIXS,0.05);
  
  SIMPLE_TRANSLATOR(TOADEXPO,"EXPTIME");
  RETURN_A_VALUE(TOADRDON,0.72);
  RETURN_A_VALUE(TOADAIRM,0.0);
  SIMPLE_TRANSLATOR(TOADUTIM,"EXPSTART");
  SIMPLE_TRANSLATOR(TOADOBJE,"TARGNAME");
  SIMPLE_TRANSLATOR(TOADGAIN,"ATODGAIN");
  SIMPLE_TRANSLATOR(TOADFILT,"FILTNAM1");  
//  RETURN_A_VALUE(TOADBAND,ToadBand(StringToUpper("FILTNAM1")));	  
  RETURN_A_VALUE(TOADBAND,ToadBand(StringToUpper(Head.KeyVal("FILTNAM1")))); 
  RETURN_A_VALUE(TOADRASC,RaDegToString(Head.KeyVal("RA_TARG")));  
  RETURN_A_VALUE(TOADDECL,DecDegToString(Head.KeyVal("DEC_TARG")));
  SIMPLE_TRANSLATOR(TOADCHIP,"DETECTOR");
  TRANSLATOR_DEC(TOADPZPT)
  {
    
    string filter = Head.KeyVal("FILTNAM1");
    double photflam =  -2.5*log10(double(Head.KeyVal("PHOTFLAM"))) ;
    double zerop = Head.KeyVal("PHOTZPT");
    zerop += photflam;

    // Gives zerop for Cousins RI
    if (filter == "F675W")
	zerop += -0.63;
    if (filter == "F814W")
      zerop += -1.22;

    return FitsKey("TOADPZPT",zerop);
  };

};

#ifdef STORAGE
FitsKey FitsHeader::HstWfpc2Format(const string &KeyName) const
{
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.1);
  if (KeyName == "TOADFILT") return KeyVal("FILTNAM1");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("ATODGAIN");
  // check the crap with PHOTFLAM or PHOTZPT for a proper calib.
  if (KeyName == "TOADRDON") return FitsKey(KeyName,0.72);  
  if (KeyName == "TOADRASC") 
      {
	  double rastr=KeyVal("RA_TARG");
	  return FitsKey(KeyName,RaDegToString(rastr));
      }
  if (KeyName == "TOADDECL")       
      {
	  double decstr=KeyVal("DEC_TARG");
	  return FitsKey(KeyName,DecDegToString(decstr));
      }

  if (KeyName == "TOADEQUI") return KeyVal("EQUINOX");
  if (KeyName == "TOADAIRM") return FitsKey(KeyName,0.0);
  if (KeyName == "TOADUTIM") return KeyVal("EXPSTART");
  if (KeyName == "TOADDATE") return KeyVal("DATE-OBS");
  if (KeyName == "TOADSCAN") 
    {
      string sec  = "not yet implemented";
      return FitsKey(KeyName,sec);
    }
  if (KeyName == "TOADILLU") 
    {
      string sec  = "not yet implemented";
      return FitsKey(KeyName,sec);
    }
  if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
  if (KeyName == "TOADOBJE") return KeyVal("TARGNAME");
  if (KeyName == "TOADCHIP") return KeyVal("DETECTOR");
  if (KeyName == "TOADBAND")
    {
      string filter = KeyVal("FILTNAM1");
      return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
    }
  return FitsKey(KeyName,NOVAL);      
}
#endif
