
class CfhtMos : public VirtualInstrument { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "CfhtMos";};
  string TelName() const {return "CFHT";}
  string InstName() const { return "Mos";}


  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { if (CheckKeyToUpper(Head,"INSTRUME","MOS")) 
    return new CfhtMos;
    return NULL;}

  RETURN_A_VALUE(TOADPIXS,0.6); // I doubt it
  SIMPLE_TRANSLATOR(TOADEXPO,"INTTIME");
  RETURN_A_VALUE(TOADBAND,end_of_string(Head.KeyVal("FILTER")));
  SIMPLE_TRANSLATOR(TOADUTIME,"UTIME");
};


#ifdef OLD_FITSTOADS
 // CFHT Mos instrument
FitsKey FitsHeader::CfhtMosFormat(const string &KeyName) const
{
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.6); // check this
  if (KeyName == "TOADFILT") return KeyVal("FILTER");
  if (KeyName == "TOADEXPO") return KeyVal("INTTIME");
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADRDON") return KeyVal("RONOISE");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADEQUI") return KeyVal("EPOCH") ;
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UTIME") ;
  if (KeyName == "TOADDATE") {
    int dd, mm, yy ;
    string kval = KeyVal("DATE-OBS") ;
    char buf[64] ;
    if(sscanf(kval.c_str(), "%d/%d/%d", &dd, &mm, &yy)) {
      sprintf(buf, "%02d/%02d/%04d", dd, mm, yy+1900) ;
      return FitsKey(KeyName, string(buf)) ;
    }
  }
  if (KeyName == "TOADSCAN")
    return FitsKey(KeyName,string("[0,0;0,0]")) ;
  if (KeyName == "TOADILLU")
    return FitsKey(KeyName,string("[0,1400;100,1300]")) ;
  if (KeyName == "TOADTYPE") return KeyVal("OBSTYPE");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return FitsKey(KeyName,0); // begining @ 0
  if (KeyName == "TOADBAND") {
    string kval = KeyVal("FILTER") ;
    return FitsKey(KeyName, ToadBand(StringToUpper(kval)));
  }
  if (KeyName == "TOADNAMP") return FitsKey(KeyName, 1) ;
  return FitsKey(KeyName,NOVAL);
}
#endif
