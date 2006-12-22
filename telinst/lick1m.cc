// This stuff is only valid for Lick 1m reduced image coming 
// from the LBL database
#ifdef VIRTUAL_INSTRUMENTS
class Lick1m : public VirtualInstrument {
public:
  string TelName() const {return "Lick1m";}

  SIMPLE_TRANSLATOR(TOADEXPO,"EXPOSURE");
  SIMPLE_TRANSLATOR(TOADUTIM,"TIME");
  SIMPLE_TRANSLATOR(TOADRDON,"RONOISE");
  RETURN_A_VALUE(TOADTYPE,"Unknown"); 

};


class Lick1mDewar5 : public Lick1m { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "Lick1mDewar5";};
  string InstName() const { return "Dewar5";}
  SIMPLE_TRANSLATOR(TOADEXPO,"EXPOSURE");
  static VirtualInstrument *Acceptor(const FitsHeader &Head)
    { 
     if (CheckKey(Head,"DET_SPEC", "li1m_dewar5")) 
       return new Lick1mDewar5;
     return NULL;
    }
};

class Lick1mDewar5_large : public Lick1m { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "Lick1mDewar5_large";};
  string InstName() const { return "Dewar5_large";}

  static VirtualInstrument *Acceptor(const FitsHeader &Head)
    { 
      double ro = (double) Head.KeyVal("RONOISE");
     if (CheckKey(Head,"DET_SPEC", "li1m_dewar5") && ro > 1.6 && ro < 1.8) 
       return new Lick1mDewar5_large;
     return NULL;
    }
};

class Lick1mDewar2 : public Lick1m { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "Lick1mDewar2";};
  string InstName() const { return "Dewar2";}

  static VirtualInstrument *Acceptor(const FitsHeader &Head)
    { 
     if (CheckKey(Head,"DET_SPEC", "li1m_dewar2") || CheckKey(Head,"DET_SPEC", "li1m_dewar2_small")) 
       return new Lick1mDewar2;
     return NULL;
    }
};
#endif

