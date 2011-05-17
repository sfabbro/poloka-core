#ifdef VIRTUAL_INSTRUMENTS
class MagellanMegacam : public VirtualInstrument { /* TYPE_SELECTOR */
public:
  static VirtualInstrument *Acceptor(const FitsHeader &Head) { 
    if (CheckKey(Head, "TELESCOP", "clay_f5" ) && CheckKey(Head,"DETECTOR","megacam")) 
      return new MagellanMegacam;
    return NULL;
  }

  string TelName() const { return "Magellan"; }
  string InstName() const { return "Megacam"; }
  string TelInstName() const { return "MagellanMegacam"; }

  SIMPLE_TRANSLATOR(TOADMJD,"MJD");
  SIMPLE_TRANSLATOR(TOADUTIM,"UT");
  SIMPLE_TRANSLATOR(TOADPIXS,"SECPIX2");
  SIMPLE_TRANSLATOR(TOADBAND,"BFILTNIC");
  SIMPLE_TRANSLATOR(TOADFILT,"BFILTID");
  RETURN_A_VALUE(TOADCHIP, 1);

};
#endif

