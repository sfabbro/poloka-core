/*
did no check anything (miss images !)

 */

class WiynS2kb : public VirtualInstrument {  /* TYPE_SELECTOR */

public :

  string TelInstName() const {return "WiynS2kb";}
  string InstName() const { return "S2kb";};
  string TelName() const { return "WIYN";};

  static VirtualInstrument* Acceptor(const FitsHeader &Head)
  { if (CheckKey(Head,"DETECTOR","s2kb")) return new WiynS2kb; return NULL;}						  
  RETURN_A_VALUE(TOADPIXS,0.197);
  SIMPLE_TRANSLATOR(TOADFILT,"FILTERS");
  TRANSLATOR_DEC(TOADRASC) // see Seb's comments below
  {
    string rastr = Head.KeyVal("RA");
    double radeg = RaStringToDeg(rastr);
    while (radeg > 360) radeg -= 360; // !!!
    return FitsKey("TOADRASC",RaDegToString(radeg));
  }
  RETURN_A_VALUE(TOADBAND,ToadBand(StringToUpper(string(Head.KeyVal("FILTERS")))));
  TRANSLATOR_DEC(TOADEQUI) // EQUINOX is the observation date.. (old fits standard)
  {
    if (CheckKey(Head,"RADECSYS","FK5")) return FitsKey("TOADEQUI",2000.);
    if (CheckKey(Head,"RADECSYS","FK4")) return FitsKey("TOADEQUI",1950.);
    return Head.KeyVal("EQUINOX");
  }
  
};

#ifdef OLD_FITSTOAD
FitsKey FitsHeader::WiynS2kbFormat(const string &KeyName) const
{
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.197);
  if (KeyName == "TOADINST") return KeyVal("DETECTOR");
  if (KeyName == "TOADFILT") return KeyVal("FILTERS");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADRDON") return KeyVal("RDNOISE");
  /*
    Below is an example of problem associated with the wiyn: watch coordinates and equinox, they don't match.
    We then choose the RA and DEC keys, they seem more exact.
    RA      = '28:55:57.70'        / Telescope RA
    DEC     = '-3:36:57.68'        / Telescope DEC
    RAOFFST = '00:00:0.00'         / Telescope RA offset
    DECOFFST= '00:00:0.00'         / Telescope DEC offset
    RADECSYS= 'FK5'                / coordinate system
    TARGRA  = '04:56:11.60'        / right ascension
    TARGDEC = '-3:41:26.00'        / declination
    EQUINOX =               1998.1 / equinox of position
    EPOCH   =               1998.1 / same as EQUINOX (for back compat.)
  */
  if (KeyName == "TOADRASC") 
    {
      string rastr = KeyVal("RA");
      double radeg = RaStringToDeg(rastr);
      while (radeg > 360) radeg -= 360;
      return FitsKey(KeyName,RaDegToString(radeg));
    }
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADEQUI") /* return KeyVal("EQUINOX"); */
    {
      string radecsys = KeyVal("RADECSYS");
      //should be 2000 but works only with 1950...should really learn the coordinate systems
      if (strstr(radecsys.c_str(),"FK5")) return FitsKey(KeyName,2000.0);
      if (strstr(radecsys.c_str(),"FK4")) return FitsKey(KeyName,2000.0); // never seen that, according to RA and DEC keys
      return FitsKey(KeyName,2000.0);
    } 
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UT");
  if (KeyName == "TOADDATE") 
    {  
      int dd,mm,yy;
      string keyval = KeyVal("DATE-OBS");
      if (sscanf(keyval.c_str(),"%d/%d/%d",&dd,&mm,&yy) == 3)
 	{
 	  yy = yy + 1900;
          char date_string[64];
	  sprintf(date_string,"%d/%d/%d",dd,mm,yy);
 	  return FitsKey(KeyName,string(date_string));
 	}
    }
  if (KeyName == "TOADSCAN") 
    {
      string keyval = KeyVal("BIASSEC");
      int x0,x1,y0,y1,nx,ny;
      if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
	{
	  nx = x1-x0+1;
	  ny = y1-y0+1;
	  char sec_charstar[64];
	  sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	  return FitsKey(KeyName, string(sec_charstar));
	}
    }
  if (KeyName == "TOADILLU") 
    {
      string keyval = KeyVal("DATASEC");
      int x0,x1,y0,y1,nx,ny;
      if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
      {
	nx = x1-x0+1;
	ny = y1-y0+1;
	char sec_charstar[64];
	sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	return FitsKey(KeyName, string(sec_charstar));
      }
    }
  if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return FitsKey(string("TOADCHIP"),"1");
  if (KeyName == "TOADBAND")
    {
      string filter = KeyVal("FILTERS");
      return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
    }
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);

  return FitsKey(KeyName,NOVAL);      
}

#endif
