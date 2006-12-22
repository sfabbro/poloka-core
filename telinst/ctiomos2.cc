/*
did no check anything (miss images !)

 */
#ifdef VIRTUAL_INSTRUMENTS
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
#endif
