/*
did no check anything (miss images !)

 */

class CtioMosaic2 : public VirtualInstrument {  /* TYPE_SELECTOR */

public :

  string TelInstName() const {return "CtioMosaic2";};
  string InstName() const { return "Mosaic2";};
  string TelName() const { return "CTIO";};


  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { if (CheckKey(Head,"DETECTOR","Mosaic2")) return new CtioMosaic2; return NULL;}						  
  SIMPLE_TRANSLATOR(TOADPIXS,"PIXSCAL1");
  SIMPLE_TRANSLATOR(TOADUTIM,"TIME-OBS");
  //  RETURN_A_VALUE(TOADDATE,toads_date(string(Head.KeyVal("DATEOBS"))));
  SIMPLE_TRANSLATOR(TOADTYPE,"OBSTYPE");
  SIMPLE_TRANSLATOR(TOADCHIP,"CHIPID");
  RETURN_A_VALUE(TOADBAND,string(Head.KeyVal("FILTER")).substr(0,1));
  TRANSLATOR_DEC(TOADRDON)
  {
    string RoNoise = Head.KeyVal("RON");
    //cout << "toto" << endl;
    double RoN = atof(RoNoise.c_str());
    //cout << RoN << endl;
    
    if (RoN>0)
      return FitsKey("TOADRDON",RoNoise);
    
    
    RoNoise = static_cast<string>(Head.KeyVal("RDNOISE"));
    RoN = atof(RoNoise.c_str());
    if (RoN>0)
      return FitsKey("TOADRDON",RoN);

    RoN = 8;

    return FitsKey("TOADRDON",RoN);
  }



};

#ifdef OLD_FITSTOAD
FitsKey FitsHeader::CtioMosaic2Format(const string &KeyName) const
{
  if (KeyName == "TOADPIXS") {
    double xpxsz = KeyVal("PIXSIZE1") ;
    double ypxsz = KeyVal("PIXSIZE2") ;
    return FitsKey(KeyName,0.5*(xpxsz+ypxsz));
  }
  if (KeyName == "TOADINST") return KeyVal("DETECTOR");
  if (KeyName == "TOADFILT") return KeyVal("FILTER");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADRDON") return KeyVal("RDNOISE");
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  
  if (KeyName == "TOADEQUI") return KeyVal("EQUINOX");
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("TIME-OBS");
  if (KeyName == "TOADDATE") 
    {
      int dd,mm,yy;
      string keyval = KeyVal("DATEOBS");
      if (sscanf(keyval.c_str(),"%d-%d-%d",& yy,&mm,&dd) == 3)
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
  if (KeyName == "TOADTYPE") return KeyVal("OBSTYPE");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT") ;
  if (KeyName == "TOADCHIP") {return KeyVal("CHIPID");} // 
  if (KeyName == "TOADBAND")
    {
      string Filter = KeyVal("FILTER");
      return FitsKey(KeyName,Filter.substr(0,1)) ;
    }
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);
  
  return FitsKey(KeyName,NOVAL);
  
}
#endif
