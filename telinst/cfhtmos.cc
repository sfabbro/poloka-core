#ifdef VIRTUAL_INSTRUMENTS
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
#endif
