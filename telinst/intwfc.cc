
/* 3 leaf classes : 
   Int1Wfc (1 chip)
   Int4Wfc (4chips)
   IntWfcNewDaq (4 chips after daq upgrade in 99)
*/

/* CHECK : 1<=(TOADCHIP,IMAGEID)<=4 */
/* Illu overscan and GuessLinWcs are implemented for the 4-chip 
setup as it was in 2000.
*/

class IntWfc : public VirtualInstrument {    /* TYPE_SELECTOR */
public :
  string TelName() const {return "INT";};
  static VirtualInstrument *Acceptor(const FitsHeader &Head); // at the end


  SIMPLE_TRANSLATOR(TOADPIXS,"SECPPIX")
  TRANSLATOR_DEC(TOADFILT)
    {return FitsKey("TOADFILT", string(Head.KeyVal("WFFBAND"))+string(Head.KeyVal("WFFPSYS")));}
  virtual SIMPLE_TRANSLATOR(TOADRDON,"READNOIS") // overloaded for IntWfcNewDaq
  SIMPLE_TRANSLATOR(TOADRASC,"CAT-RA")
  SIMPLE_TRANSLATOR(TOADDECL,"CAT-DEC")
  TRANSLATOR_DEC(TOADEQUI)
    {
      string keyval = Head.KeyVal("CAT-EQUI");
      if (strstr(keyval.c_str(), "B1950")) return FitsKey("TOADEQUI",1950.0); 
      return FitsKey("TOADEQUI",2000.0); 
    }
  SIMPLE_TRANSLATOR(TOADUTIM,"UTSTART");
  TRANSLATOR_DEC(TOADBAND)
    { 
      string filter = Head.KeyVal("WFFBAND");
      return FitsKey("TOADBAND",ToadBand(StringToUpper(filter)));	  
    }

  bool GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const
  {  cerr << " GuessWCS still to be impemented for IntWfc " << endl;return false; }
};




class Int4Wfc : public IntWfc  /* TYPE_TESTER */
{
  string TelInstName() const { return "Int4Wfc";}
  string InstName() const { return "Wfc(4chips)";}


  virtual FitsKey TOADCHIP(const FitsHeader &Head, const bool Warn) const
  {
    int chip = Head.KeyVal("DASCHAN", Warn);
    if (chip == 5) chip = 2;
    return FitsKey("TOADCHIP",chip);
  };

  bool GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const;
  // Frame  SkyRegion(const FitsHeader &Head) const;

};

bool Int4Wfc::GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const
{
  int id = Head.KeyVal("TOADCHIP");
  switch (id)
    {
    case 1: ComputeLinWCS(Head,Point(-1120 ,1980.), RotationFlip(Left,Down), Guess); break;
    case 2: ComputeLinWCS(Head,Point( 4200 , 980.), RotationFlip(Down,Right),Guess); break;
    case 3: ComputeLinWCS(Head,Point( 3054.,1978.), RotationFlip(Left,Down), Guess); break;
    case 4: ComputeLinWCS(Head,Point(  993, 1995.), RotationFlip(Left,Down), Guess); break;
    default: VirtualInstrument::GuessLinWCS(Head,Guess); break;
    }
  return true;
}

#ifdef TENTATIVE
Frame IntWfcNewDaq::SkyRegion(const FitsHeader &Head) const
{
  double ra,dec;
  RaDec2000(Head, ra, dec);  
  double cosdec = cos(dec*M_PI/180.);
  return Frame(ra-0.38/cosdec, dec-0.30, ra+0.22/cosdec, dec+0.30);
}
#endif



class Int1Wfc : public IntWfc
{
  string TelInstName() const {return "Int1Wfc";}
  string InstName () const { return "Wfc1chip";}
  bool GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const
  {  cerr << " GuessLinWCS still to be impemented for Int1Wfc" << endl;return false; }
};  


class IntWfcNewDaq: public Int4Wfc {      /*  TYPE_TESTER */
public :
  string TelInstName() const {return "IntWfcNewDaq";}
  string InstName () const { return "WfcNewDaq";}
  SIMPLE_TRANSLATOR(TOADCHIP,"IMAGEID");

  FitsKey TOADRDON(const FitsHeader &Head, const bool Warn) const
  {
    double noise;
    if (Head.HasKey("RDNOISE") && ( (noise = Head.KeyVal("RDNOISE")) > 0. ))
	{
	    return FitsKey("TOADRDON",noise);
	}
    else 
      { 
	static double ronoise[4] = {11,13,9,11}; 
	return FitsKey("TOADRDON",ronoise[int(Head.KeyVal("IMAGEID"))-1]);
      }
    }

  Frame TotalIlluRegion(const FitsHeader &Head) const
  {Frame result; 
  fits_imregion_to_frame(string(Head.KeyVal("TRIMSEC")), result);
  return result;
  };

  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  {
    string oscan;
    oscan = "[2110:2150,1:2048]";    
    Frame result;
    fits_imregion_to_frame(oscan, result);
    return result;
  }
  

};



VirtualInstrument* IntWfc::Acceptor(const FitsHeader &Head)
{
  string KEYINST = StringToUpper(string(Head.KeyVal("INSTRUME")));
  if (KEYINST != "WFC") return NULL;
  double jd = Head.KeyVal("JD");
  if (Head.HasKey("JD") &&  jd > 2451410.5822505) /* roughly Aug 20th 1999 when the Daq was changed */
    { return new(IntWfcNewDaq);}
  if (Head.HasKey("CCDTYPE"))
    {
      string keyccd  = Head.KeyVal("CCDTYPE");
      keyccd = StringToLower(keyccd);
      if (keyccd=="loral") {return new Int1Wfc;}
      if (keyccd=="eev42") {return new Int4Wfc;}
    }
  // WARNING ! THIS IS THE ONLY WAY WE FOUND TO SEPARATE 
  // Int1Wfc from Int4Wfc using main header keywords..
  else if (Head.HasKey("JD"))
    {
      if (jd > 2451000.0) return new Int4Wfc; 
      else return new Int1Wfc; 
    }

return NULL;
}

#ifdef OLD_FITSTOAD

//*********************************************************************
FitsKey FitsHeader::IntWfcFormat(const string &KeyName) const
{
  if (KeyName == "TOADPIXS") return KeyVal("SECPPIX");
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADFILT") return FitsKey(KeyName, string(KeyVal("WFFBAND"))+string(KeyVal("WFFPSYS")));
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADRDON")
    {
    if (TelInst != IntWfcNewDaq) return KeyVal("READNOIS");
    // le RDNOISE est rempli a zero (31/08/99)
    if (HasKey("RDNOISE") && ( (double) ( KeyVal("RDNOISE")) > 0. ))
	{
	    return KeyVal("RDNOISE");
	}
    else 
      { 
	static double ronoise[4] = {11,13,9,11}; 
	return FitsKey(KeyName,ronoise[int(KeyVal("IMAGEID"))+1]);
      }
    }
  if (KeyName == "TOADRASC") return KeyVal("CAT-RA");
  if (KeyName == "TOADDECL") return KeyVal("CAT-DEC");
  if (KeyName == "TOADEQUI") 
    {
      string keyval = KeyVal("CAT-EQUI");
      if (strstr(keyval.c_str(), "B1950")) return FitsKey(KeyName,1950.0); 
      return FitsKey(KeyName,2000.0); 
    }
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UTSTART");
  if (KeyName == "TOADDATE")
    {
      int dd,mm,yy;
      string keyval = KeyVal("DATE-OBS");
      if (sscanf(keyval.c_str(),"%d-%d-%d",&yy,&mm,&dd) == 3)
 	{ 
	  char date_charstar[64];
	  sprintf(date_charstar,"%d/%d/%d",dd,mm,yy);
	  return FitsKey(KeyName,string(date_charstar));
	}
    } 
  if (KeyName == "TOADSCAN") 
    {
      if (HasKey("BIASSEC"))
	{
/*	  string keyval = KeyVal("BIASSEC");
	  int x0, x1, y0, y1;
	  if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
	    {
	      int nx = x1-x0+1;
	      int ny = y1-y0+1;
	      char sec_charstar[64];
	      sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	      return FitsKey(KeyName,string(sec_charstar));
	    }*/

	  return FitsKey(KeyName,"[10,2150;4105,4190]");

	}
      else /* BIASSEC absent on 31/Aug/99 */
      return FitsKey(KeyName,"[5,48;1,4200]");
 
      if (HasKey("DATASEC"))

	{

	  return FitsKey(KeyName, "[-52,2101;1,4200]");

	}

    }
  if (KeyName == "TOADILLU") 
    {
      if (HasKey("TRIMSEC"))
	{
	  string keyval = KeyVal("TRIMSEC");
	  int x0,x1,y0,y1,nx,ny;
	  if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
	    {
	      nx = x1-x0+1;
	      ny = y1-y0+1;
	      char sec_charstar[64];
	      sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	      return FitsKey(KeyName,string(sec_charstar)) ;
	    }
	}
      else /* TRIMSEC absent on 31/Aug/99 */
        return FitsKey(KeyName,"[54,2048;1,4100]");
    }
  if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP")
    {
    if (TelInst == IntWfcNewDaq)
      {
	return FitsKey("IMAGEID",fptr); 
        // KeyVal("IMAGEID") has the same effect, unless the key is absent and one prints it out
      }
    else
      {
	int ccd = KeyVal("DASCHAN");
	if (ccd==5) ccd = 2;
	return FitsKey(KeyName,ccd);
      }
    } 
  if (KeyName == "TOADBAND")
    {
      string filter = KeyVal("WFFBAND");
      return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
    }
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);

  return FitsKey(KeyName,NOVAL);      
}

#endif /*STORAGE */
