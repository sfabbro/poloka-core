#ifdef VIRTUAL_INSTRUMENTS
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
#endif
