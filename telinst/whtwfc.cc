
class WhtWfc : public VirtualInstrument { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "WhtWfc";};
  string TelName() const {return "WHT";};
  string InstName() const { return "WFC";};
  
  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { 
    if (CheckKeyToUpper(Head,"DETECTOR","WHTWFC")) 
      return new WhtWfc;
    return NULL;
  }
  
  SIMPLE_TRANSLATOR(TOADRDON,"READNOIS");
  SIMPLE_TRANSLATOR(TOADFILT,"PFAFILT");
  // as find on the web http://www.ing.iac.es/Astronomy/observing/instruments.html
  RETURN_A_VALUE(TOADPIXS,0.24);
};
