class KeckIIesi : public VirtualInstrument { /* TYPE_SELECTOR */

public :
  
  string TelInstName() const {return "KeckIIesi";}
  string InstName() const { return "esi";}
  string TelName() const {return "KECKII";}
  
  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { 
    const char *p =StringToUpper(Head.KeyVal("INSTRUME")).c_str(); 
    if  (strstr(p,"ESI")==p) return new KeckIIesi;
    else return NULL;
  }
  
  RETURN_A_VALUE(TOADPIXS,0.15);
  RETURN_A_VALUE(TOADNAMP,1);
  TRANSLATOR_DEC(TOADGAIN)
  {
    if (Head.HasKey("GAIN")) return FitsKey("TOADGAIN",string(Head.KeyVal("GAIN"))); 
    //    else return FitsKey("TOADGAIN",(1.97+2.1)/2.0);
    else return FitsKey("TOADGAIN",(1.97+2.1)/2.0);
  }
  RETURN_A_VALUE(TOADRDON,2.5); 
  SIMPLE_TRANSLATOR(TOADEXPO,"EXPOSURE");
  SIMPLE_TRANSLATOR(TOADTYPE,"OBJECT");
  SIMPLE_TRANSLATOR(TOADFILT,"DWFILNAM");
  SIMPLE_TRANSLATOR(TOADBAND,"DWFILNAM");
  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  {
    string keyval = "[3:12,1:1550]";
    Frame result;
    fits_imregion_to_frame(keyval, result);
    return result;
  }

  Frame TotalIlluRegion(const FitsHeader &Head) const
  {
    string illu ="[201:800,101:1400]";
    // This illu region has an impact on the GuessLinWCS routine

    Frame result;
    fits_imregion_to_frame(illu, result);
    return result;
  }

  bool GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const
  {
    //double angle =double(Head.KeyVal("ROTDEST"))*M_PI/180.0 + (M_PI/2.0);
    double angle =double(Head.KeyVal("ROTDEST"));
    double xc,yc;

    if (angle>45.0) 
      {
	angle = M_PI;
	// The values depend on the illuregion selected.
	xc = 310;
	yc = 975;
      }
    else
      {
	angle = M_PI/2.0;
	xc = 300;
	yc = 650;
      }
	
    Point center(xc,yc);
    double c = cos(angle);
    double s = sin(angle);
    GtransfoLin rotate(0,0,c,-s,s,c);
    
#ifdef DEBUG
    GtransfoLin rotate(0,0,-1,0,0,-1);
#endif
      
    return ComputeLinWCS (Head, center, rotate, Guess);
  }

  Frame SkyRegion(const FitsHeader &Head) const
  {
    int nx,ny;
    Head.ImageSizes(nx,ny);
    TanPix2RaDec Pix2RaDec;
    if (!GuessLinWCS(Head, Pix2RaDec)) return Frame();
    cout << " Lin WCS Guess " << Pix2RaDec << endl;
    Point p00 = Pix2RaDec.apply(Point(0,0));
    Point p11 = Pix2RaDec.apply(Point(nx,ny));
    Frame skyFrame(p00,p11);
    double factor = 1.0;
    if (getenv("DILATE_FRAME")) factor = atof(getenv("DILATE_FRAME"));
    return skyFrame.Rescale(factor);
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
