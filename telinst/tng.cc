class TngLrs : public VirtualInstrument { /* TYPE_SELECTOR */
public:
  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  {if (CheckKey(Head,"INSTRUME","LRS")) 
   return new TngLrs;
   return NULL;}

  string TelName() const {return "Tng";}
  string InstName() const {return "Lrs";}
  string TelInstName() const {return "TngLrs";}

  //SIMPLE_TRANSLATOR(TOADFILT,"JAGFBAND");
  SIMPLE_TRANSLATOR(TOADEQUI,"EQX-OBS");
  SIMPLE_TRANSLATOR(TOADUTIM,"EXPSTART");
  SIMPLE_TRANSLATOR(TOADEXPO,"EXPTIME");
  SIMPLE_TRANSLATOR(TOADGAIN,"GAIN_1");
  SIMPLE_TRANSLATOR(TOADRDON,"RMS_1");

  //taken from http://www.tng.iac.es/instruments/lrs/lrs.html
  RETURN_A_VALUE(TOADCHIP,1);
  RETURN_A_VALUE(TOADPIXS,0.275);
   
  RETURN_A_VALUE(TOADDECL,DecDegToString(Head.KeyVal("DEC-DEG")));
  TRANSLATOR_DEC(TOADRASC)
    {
      double ra = Head.KeyVal("RA-DEG");
      string RA = RaDegToString(15.*ra);
      return FitsKey("TOADRASC",RA);
      }

  TRANSLATOR_DEC(TOADFILT)
  {
    string keyval = Head.KeyVal("GRM_ID");
    string keyval2 = Head.KeyVal("FLT_ID");
    keyval = keyval[keyval.length()-1];
    //SLT_ID may not exist
    if (keyval2 == "") return FitsKey("TOADFILT",keyval);
    else return FitsKey("TOADFILT",keyval2);
  }

  TRANSLATOR_DEC(TOADTYPE)
  {
    //is it a spectrum?
    if (string(Head.KeyVal("SLT_ID")) == "NULL" || string(Head.KeyVal("SLT_ID")) == "OPEN"){

      string keyval = Head.KeyVal("GRM_ID");
      string keyval2 = Head.KeyVal("OBS-TYPE");

      //OBS-TYPE does not exist
      if (keyval2 == ""){

	//lamp flats?
	if ((keyval.find("LAMP",0) != string::npos)||(keyval.find("lamp",0) != string::npos)){
	  return FitsKey("TOADTYPE","undefined");
	}
	//flats
	if (keyval.find("FLAT",0) != string::npos || keyval.find("flat",0) != string::npos || keyval.find("Flat",0) != string::npos){
	  return FitsKey("TOADTYPE","flat"); 
	}
	//bias
	if (keyval.find("BIAS",0) != string::npos || keyval.find("bias",0) != string::npos || keyval.find("NULL",0) != string::npos){
	  if (double(Head.KeyVal("EXPTIME")) == 0.0) return FitsKey("TOADTYPE","bias");
	  //bias with exptime > 0s?
	  else return FitsKey("TOADTYPE","bias_error");
	}
	//spectra...
	if (keyval.find("Helium",0) != string::npos || keyval.find("Hz",0) != string::npos){
	  return FitsKey("TOADTYPE","spectrum");
	}
	//all the rest
	return FitsKey("TOADTYPE","object");
      }
      //OBS-TYPE exists
      else{
	if (keyval2.find("BIAS",0) != string::npos && double(Head.KeyVal("EXPTIME")) != 0.0) return FitsKey("TOADTYPE","bias_error");
	else return FitsKey("TOADTYPE",keyval2);
      }
    }
    else return FitsKey("TOADTYPE","spectrum");
  }
  
  TRANSLATOR_DEC(TOADOBJE)
  {
    string keyval = Head.KeyVal("OBJCAT");
    string keyval2 = Head.KeyVal("OBJECT");
    string keyval3 = Head.KeyVal("GRM_ID");

    //    if (keyval == "")
    //      {
    //	if (keyval2 == "") return FitsKey("TOADOBJE",keyval3);
    //	else return FitsKey("TOADOBJE",keyval2);
    //     }
    //   else return FitsKey("TOADOBJE",keyval);
    return FitsKey("TOADOBJE",keyval+" "+keyval2+" "+keyval3);
  }

  
  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  {
    string oscan;
    oscan = "[2071:2088,1:2050]";    
    Frame result;
    fits_imregion_to_frame(oscan, result);
    return result;
  }
  
  Frame IlluRegion(const FitsHeader &Head, const int Iamp) const
  {
    return TotalIlluRegion(Head);
  }

  Frame TotalIlluRegion(const FitsHeader &Head) const
  {
    int nx = Head.KeyVal("NAXIS1");
    int ny = Head.KeyVal("NAXIS2");
    string illu;
    if ((nx==2000)&&(ny==2000)) illu ="[1:2000,1:2000]";
      else  illu ="[21:2020,21:2020]";
    // keep only the non vigneted part of the image
    //   illu ="[1:2069,1:2120]";
    Frame result;
    fits_imregion_to_frame(illu, result);
    return result;
  }

  bool GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const;

};

bool TngLrs::GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const
{
// Ra,Dec at image center, North = down, East = left
  return ComputeLinWCS(Head,Head.ImageCenter(),
		       RotationFlip(Down,Left),
		       Guess);
}


