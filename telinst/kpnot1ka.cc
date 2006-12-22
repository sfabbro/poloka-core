#ifdef VIRTUAL_INSTRUMENTS
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
#endif
