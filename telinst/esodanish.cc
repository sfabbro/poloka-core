


class DanishDfosc : VirtualInstrument { /* TYPE_SELECTOR */

public :
  string TelInstName() const {return "DanishDfosc";}
  string TelName () const {return "EcoDanish";}
  string InstName () const { return "Dfosc";}

  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  {if (CheckKeyToUpper(Head,"INSTRUME","DFOSC")) return new DanishDfosc; return NULL;};

  RETURN_A_VALUE(TOADPIXS,0.39);
  SIMPLE_TRANSLATOR(TOADFILT,"OPTI-2NM"); 
  SIMPLE_TRANSLATOR(TOADRDON,"RONOISE");
  // Do not put any specifit TOADDATE converter ... check if it works for both 1900 and 2000. (p.A)
  SIMPLE_TRANSLATOR(TOADUTIM,"TM-START"); // do you really need to average with Tm-end?
  SIMPLE_TRANSLATOR(TOADTYPE,"EXPOTYPE");
  RETURN_A_VALUE(TOADBAND,StringToUpper(string(Head.KeyVal("OPTI-2NM"))));
};



#ifdef OLD_TOADS



// ESO/Danish 1.54m 
FitsKey FitsHeader::DanishDfoscFormat(const string &KeyName) const
{
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.390);
  if (KeyName == "TOADFILT") return KeyVal("OPTI-2NM");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADRDON") return KeyVal("RONOISE");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADEQUI") return KeyVal("EPOCH") ;
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") {
      float tm_st, tm_end, tm ;
      tm_st = (float)KeyVal("TM_START") ;
      tm_end = (float)KeyVal("TM_END") ;
      tm = 0.5*(tm_st+tm_end)/3600. ;
      string ret = UtDeciToString(tm) ;
      return FitsKey(KeyName, ret) ;
  }
  if (KeyName == "TOADDATE") { 
      string date_str = KeyVal("DATE-OBS") ;
      int dd, mm, yy ;
      char buf[64] ;
      // Date format changed between 1999 and 2000
      // that's scary
      if(date_str.find('/') != string::npos) {
	  if( sscanf(date_str.c_str(), "%d/%d/%d", &dd, &mm, &yy) ) {
	      if(yy<50) yy += 2000 ; else yy += 1900 ;
	      sprintf(buf, "%02d/%02d/%04d", dd, mm, yy) ;
	      return FitsKey(KeyName, string(buf)) ;
	  }
      }
      else {
	  if( sscanf(date_str.c_str(), "%d-%d-%d", &yy, &mm, &dd) ) {
	      sprintf(buf, "%02d/%02d/%04d", dd, mm, yy) ;
	      return FitsKey(KeyName, string(buf)) ;
	  }
      }
  }
  if (KeyName == "TOADSCAN") {
      return FitsKey(KeyName,string("[2035,65;0,2050]")) ;
  }
  if (KeyName == "TOADILLU") {
      return FitsKey(KeyName,string("[0,2034;0;2050]")) ;
  }
  if (KeyName == "TOADTYPE") return KeyVal("EXPOTYPE");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return FitsKey(KeyName,0); // begining @ 0
  if (KeyName == "TOADBAND") {
      string kval = KeyVal("OPTI-2NM") ;
      //      string ret = kval.substr(0,1);
      return FitsKey(KeyName, ToadBand(StringToUpper(kval)));
  }
  if (KeyName == "TOADNAMP") return FitsKey(KeyName, 1) ;
  return FitsKey(KeyName,NOVAL);
}


#endif
