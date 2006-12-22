#ifdef VIRTUAL_INSTRUMENTS
class KeckILris : public VirtualInstrument { /* TYPE_SELECTOR */

public :
  
  string TelInstName() const {return "KeckILris";}
  string InstName() const { return "Lris";}
  string TelName() const {return "KECKI";}
  
  static VirtualInstrument* Acceptor(const FitsHeader &Head)
  { if (CheckKey(Head,"INSTRUME","LRIS")) return new KeckILris; return NULL;}
  
  RETURN_A_VALUE(TOADPIXS,0.215);
  RETURN_A_VALUE(TOADNAMP,2);
  TRANSLATOR_DEC(TOADGAIN)
  {
    if (Head.HasKey("GAIN")) return FitsKey("TOADGAIN",string(Head.KeyVal("GAIN"))); 
    //    else return FitsKey("TOADGAIN",(1.97+2.1)/2.0);
    else return FitsKey("TOADGAIN",(1.97+2.1)/2.0);
  }  
  RETURN_A_VALUE(TOADRDON,(6.3+6.6)/2.0); 
  SIMPLE_TRANSLATOR(TOADEXPO,"EXPOSURE");
  SIMPLE_TRANSLATOR(TOADTYPE,"OBJECT");
  SIMPLE_TRANSLATOR(TOADFILT,"REDFILT");
  SIMPLE_TRANSLATOR(TOADBAND,"REDFILT");
  //  RETURN_A_VALUE(TOADBAND, ToadBand(StringToUpper(Head.KeyVal("REDFILT"))));	  
  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  {
    if (Iamp==1)
      {
	string keyval = "[1:234,1:2048]";
	Frame result;
	fits_imregion_to_frame(keyval, result);
	return result;
      }
    if (Iamp==2)
      {
	string keyval = "[1801:2250,1:2048]";
	Frame result;
	fits_imregion_to_frame(keyval, result);
	return result;
      }
    else
      {
	cerr <<Head.FileName()<<": Iamp = "<<Iamp<<" and be 1 or 2 !!!!!"<< endl;
	return Frame();
      }
  }
  Frame IlluRegion(const FitsHeader &Head, const int Iamp) const
  {
    if (Iamp == 1)
      {
	string illu ="[240:1068,1:2048]";
	Frame result;
	fits_imregion_to_frame(illu, result);
	return result;
      }
    if (Iamp == 2)
      {
	string illu ="[1068:1800,1:2048]";
	Frame result;
	fits_imregion_to_frame(illu, result);
	return result;
      }
    return TotalIlluRegion(Head);
  }



  Frame TotalIlluRegion(const FitsHeader &Head) const
  {
    string illu ="[240:1820,1:2048]";
    Frame result;
    fits_imregion_to_frame(illu, result);
    return result;
  }
};
#endif














