

class CtioBtc : public VirtualInstrument { /* TYPE_SELECTOR */
public:
  string TelInstName() const { return "CtioBtc";};
  string InstName() const { return "Btc";};
  string TelName() const { return "CTIO";};


  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { if (CheckKeyToUpper(Head,"INSTRUME","BTC")) return new CtioBtc; return NULL;}

  RETURN_A_VALUE(TOADPIXS,0.425);
  SIMPLE_TRANSLATOR(TOADRDON, "RON");

  TRANSLATOR_DEC(TOADRASC)
  {
      // if LBL changed the RA in the header!!! recognize the fools with the welldept key
      // actually sometimes they didnot, depending on who did the reduction. Usually
      // they did it with image name suffix cln and it seems that they have the key CCDPROC
      // 01/09/00: not always. how to get out of it?? OVERSCAN seems to be ok
      //      if (HasKey("WELLDEPT") && HasKey("OVERSCAN")) return KeyVal("RA");
      int chip = Head.KeyVal("AMPID");
      if ((chip == 2) || (chip == 3) )
	{
	  string raString = Head.KeyVal("RA");
	  double ra = RaStringToDeg(raString);
	  string decString = Head.KeyVal("DEC");
	  double dec = DecStringToDeg(decString);
	  double cosfactor = cos(dec*M_PI/180)*60.;
       	  return FitsKey("TOADRASC",RaDegToString(ra + 20.4/cosfactor));
	}
      else return Head.KeyVal("RA");
    };

  TRANSLATOR_DEC(TOADDECL)
  {
      // same comment as RA
      if (Head.HasKey("WELLDEPT") && Head.HasKey("OVERSCAN")) return Head.KeyVal("DEC");
      int chip = Head.KeyVal("AMPID");
      if ((chip == 3) || (chip == 4))
	{
	  string decString = Head.KeyVal("DEC");
	  double dec = DecStringToDeg(decString);
	  return FitsKey("TOADDECL",DecDegToString(dec - 20.4/60.));
	}
      else  return Head.KeyVal("DEC");
  }
  SIMPLE_TRANSLATOR(TOADEQUI,"EPOCH");
  FitsKey TOADOBJE(const FitsHeader &Head, const bool Warn) const
  {
    string obje =  Head.KeyVal("TITLE");
    RemovePattern(obje," ");
    if (obje.length() == 0) obje = static_cast<string>(Head.KeyVal("OBJECT"));
    return FitsKey("TOADOBJE",obje);
  }
  SIMPLE_TRANSLATOR(TOADCHIP,"AMPID");

};
 
#ifdef OLD_FITSTOAD
FitsKey FitsHeader::CtioBtcFormat(const string &KeyName) const
{
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.425);
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADFILT") return KeyVal("FILTER");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADRDON") return KeyVal("RON");
  if (KeyName == "TOADRASC") 
    {
      // if LBL changed the RA in the header!!! recognize the fools with the welldept key
      // actually sometimes they didnot, depending on who did the reduction. Usually
      // they did it with image name suffix cln and it seems that they have the key CCDPROC
      // 01/09/00: not always. how to get out of it?? OVERSCAN seems to be ok
      //      if (HasKey("WELLDEPT") && HasKey("OVERSCAN")) return KeyVal("RA");
      int chip = KeyVal("AMPID");
      if ((chip == 2) || (chip == 3) )
	{
	  string raString = KeyVal("RA");
	  double ra = RaStringToDeg(raString);
	  string decString = KeyVal("DEC");
	  double dec = DecStringToDeg(decString);
	  double cosfactor = cos(dec*M_PI/180)*60.;
       	  return FitsKey(KeyName,RaDegToString(ra + 20.4/cosfactor));
	}
      else return KeyVal("RA");
    }

  if (KeyName == "TOADDECL") 
    {
      // same comment as RA
      if (HasKey("WELLDEPT") && HasKey("OVERSCAN")) return KeyVal("DEC");
      int chip = KeyVal("AMPID");
      if ((chip == 3) || (chip == 4))
	{
	  string decString = KeyVal("DEC");
	  double dec = DecStringToDeg(decString);
	  return FitsKey(KeyName,DecDegToString(dec - 20.4/60.));
	}
      else  return KeyVal("DEC");
    }

  if (KeyName == "TOADEQUI") return KeyVal("EPOCH");
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UT");
  if (KeyName == "TOADDATE") 
    {
      int dd,mm,yy;
      string keyval = KeyVal("DATEOBS");
      if (sscanf(keyval.c_str(),"%d/%d/%d",&mm,&dd,&yy) == 3)
 	{
          char date_string[64];
	  sprintf(date_string,"%d/%d/%d",dd,mm,yy);
 	  return FitsKey(KeyName,string(date_string));
 	}
    };
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
	  return FitsKey(KeyName,string(sec_charstar));
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
  if (KeyName == "TOADOBJE") 
    {
      string obje =  KeyVal("TITLE");
      RemovePattern(obje," ");
      if (obje.length() == 0) obje = KeyVal("OBJECT");
      return FitsKey(KeyName,obje);
    }
  if (KeyName == "TOADCHIP") {return KeyVal("AMPID");}
  if (KeyName == "TOADBAND")
    {
      string Filter = KeyVal("FILTER");       
      return FitsKey(KeyName,ToadBand(StringToUpper(Filter)));	  
    }
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);
  
  return FitsKey(KeyName,NOVAL);
  
}
#endif
