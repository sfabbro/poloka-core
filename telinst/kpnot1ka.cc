class Kpno2p1mT1ka : public VirtualInstrument { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "Kpno2p1mT1ka";};
  string TelName() const {return "KPNO 2.1m";}
  string InstName() const { return "T1ka";}

  static VirtualInstrument* Acceptor(const FitsHeader &Head)
  { if (CheckKey(Head,"DETECTOR","t1ka")) 
    return new Kpno2p1mT1ka; return NULL;}
  RETURN_A_VALUE(TOADPIXS,0.3);
  SIMPLE_TRANSLATOR(TOADEQUI,"EPOCH");
  
};

#ifdef OLD_FITSTOAD
FitsKey FitsHeader::Kpno2p1mT1kaFormat(const string& KeyName) const 
{
  if (KeyName == "TOADINST") return KeyVal("DETECTOR");
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.3); // check this
  if (KeyName == "TOADFILT") return KeyVal("FILTER");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADRDON") return KeyVal("RDNOISE");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADEQUI") return KeyVal("EPOCH") ;
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UT") ;
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
       return FitsKey(KeyName,string("[0,1024;0,1024]")) ;
   if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
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
 
