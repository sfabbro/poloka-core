#ifdef VIRTUAL_INSTRUMENTS
class NotMosca : public VirtualInstrument { /* TYPE_SELECTOR */
public:  
  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  {if ( CheckKey(Head,"INSTRUME","MOSCA") )
    return new NotMosca;
  return NULL;}

  string TelName() const {return "Not";}
  string InstName() const {return "Mosca";}
  string TelInstName() const {return "NotMosca";}

  SIMPLE_TRANSLATOR(TOADRDON,"RDNOISE");
  SIMPLE_TRANSLATOR(TOADGAIN,"GAIN");
  SIMPLE_TRANSLATOR(TOADPIXS,"PIXSCALE");
  RETURN_A_VALUE(TOADRASC,RaDegToString(Head.KeyVal("RA")));
  RETURN_A_VALUE(TOADDECL,DecDegToString(Head.KeyVal("DEC")));

  TRANSLATOR_DEC(TOADUTIM)
  {
    float ut = Head.KeyVal("UT");
    int hh = (int) floor(ut);
    int mm = (int) floor((ut-hh)*60.0);
    float ss = ((ut-hh)*60.0-mm)*60.0;
    char utim[15];
    sprintf(utim,"%02d:%02d:%4.1f",hh,mm,ss);
    return FitsKey("TOADUTIM",utim);
  }

  TRANSLATOR_DEC(TOADFILT)
  {
    // assumed that only one filter is used at a certain time
    // else only the last 'if' will be returned
    string keyval = Head.KeyVal("FILTER");
    keyval = static_cast<string>(Head.KeyVal("FAFLTNM"));
    if (!strstr(keyval.c_str(),"Open")) return FitsKey("TOADFILT",keyval);
    keyval = static_cast<string>(Head.KeyVal("FBFLTNM"));
    if (!strstr(keyval.c_str(),"Open")) return FitsKey("TOADFILT",keyval);
    return FitsKey("TOADFILT","undefined");
  }

  TRANSLATOR_DEC(TOADOBJE)
  {
    string keyval = Head.KeyVal("OBJECT");
    string keyval2 = Head.KeyVal("TCSTGT");
    if (keyval.length()!=0) return FitsKey("TOADOBJE",keyval);
    if (keyval2.length()!=0) return FitsKey("TOADOBJE",keyval2);
    return FitsKey("TOADOBJE","undefined");
  }
  
  TRANSLATOR_DEC(TOADTYPE)
  {
    string keyval = Head.KeyVal("OBJECT");
    double keyval2 = Head.KeyVal("EXPTIME");
    string keyval3 = Head.KeyVal("TCSTGT");
    if (keyval.length()!=0)
      {
	if ((strstr(keyval.c_str(),"flat")) || (strstr(keyval.c_str(),"Blank"))) return FitsKey("TOADTYPE","flat");
	if (strstr(keyval.c_str(),"bias")||strstr(keyval.c_str(),"BIAS")||keyval2==0) return FitsKey("TOADTYPE","bias");
	return FitsKey("TOADTYPE","object");
      }
    if (keyval3.length()!=0)
      {
	if ((strstr(keyval3.c_str(),"flat")) || (strstr(keyval3.c_str(),"Blank"))) return FitsKey("TOADTYPE","flat");
	if (strstr(keyval3.c_str(),"bias")||strstr(keyval3.c_str(),"BIAS")||keyval2==0) return FitsKey("TOADTYPE","bias");
	return FitsKey("TOADTYPE","object");
      }
    if (keyval2==0) return FitsKey("TOADTYPE","bias"); 
    return FitsKey("TOADTYPE","undefined");
  }

  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  {
    int nx = Head.KeyVal("NAXIS1");
    int ny = Head.KeyVal("NAXIS2");
    string oscan;
    if ((nx==2102) && (ny==2052))
      {
	oscan = "[2052:2102,1:2052]";    
      }
    else if ((nx==1051) && (ny==1026))
      {
	oscan = "[1026:1051,1:1026]";    
      }
    else
      {
	cerr <<" NAXIS1 and/or NAXIS2 NOT USUAL !!!!!!!"<< endl;
	oscan = "[0:0;0:0]";
      }
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
    if ((nx==2102) && (ny==2052))
      {
	illu ="[251:1750,251:1750]";
      }
    else if((nx==1500) && (ny==1500))
      {
	illu ="[1:1500,1:1500]";
      }
    else if ((nx==1051) && (ny==1026))
      {
	illu ="[126:875,126:875]";
      }
    else if ((nx==750) && (ny==750))
      {
	illu ="[1:750,1:750]";
      }
    else
      {
	cerr <<" NAXIS1 and/or NAXIS2 NOT USUAL !!!!!!!"<< endl;
	illu = "[0:0;0:0]";
      }
    Frame result;
    fits_imregion_to_frame(illu, result);
    return result;
  }
};
#endif
