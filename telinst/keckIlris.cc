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

#ifdef FITSTOADS
//*********************************************************************
FitsKey FitsHeader::KeckIILrisFormat(const string &KeyName) const
{
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,2);
  if (KeyName == "TOADINST") return KeyVal("PONAME");
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.215);
  if (KeyName == "TOADFILT") return KeyVal("REDFILT");
  if (KeyName == "TOADEXPO") return KeyVal("EXPOSURE");
  if (KeyName == "TOADGAIN") 
     {
	if (HasKey("GAIN")) return KeyVal("GAIN");
	else return FitsKey(KeyName,(1.97+2.1)/2.0);
     }
  if (KeyName == "TOADRDON") return FitsKey(KeyName,(6.3+6.6)/2.0);
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADEQUI") return KeyVal("EQUINOX");
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UT");
  if (KeyName == "TOADDATE") 
    {
      string keyval;
      if (HasKey("DATE-OBS")) keyval = KeyVal("DATE-OBS");
      else keyval= KeyVal("DATE");
      int dd,mm,yy;
      if (sscanf(keyval.c_str(),"%d/%d/%d",&dd,&mm,&yy) == 3)
 	{ 
	  char date_charstar[64];
	  sprintf(date_charstar,"%d/%d/%d",dd,mm,yy+1900);
	  return FitsKey(KeyName,string(date_charstar));
	}
      if (sscanf(keyval.c_str(),"%d-%d-%d",&yy,&mm,&dd) == 3)
 	{ 
	  char date_charstar[64];
	  sprintf(date_charstar,"%d/%d/%d",dd,mm,yy);
	  return FitsKey(KeyName,string(date_charstar));
	}
    } 

  if (KeyName == "TOADSCAN") 
    {
      //Only considering the right side of overscan
      string sec  = "[2091,80;1,2048]";
      return FitsKey(KeyName,sec);
    }
  if (KeyName == "TOADILLU") 
    {
      string sec  = "[240,1640;1,2048]";
      return FitsKey(KeyName,sec);
    }
  if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return FitsKey(KeyName,1);
  if (KeyName == "TOADBAND")
    {
      string filter = KeyVal("REDFILT");
      return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
    }
  return FitsKey(KeyName,NOVAL);      
}
//*********************************************************************
#endif















