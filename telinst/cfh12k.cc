

/*
  TO BE DONE : implement the necessary swaps due to a bug in the BIASSEC : done in the default decoding
*/

#include "wcsutils.h"

class Cfht12K: public VirtualInstrument {  /* TYPE_SELECTOR TYPE_TESTER */

public :
  string TelInstName() const {return "Cfht12K";}
  string TelName () const {return "CFHT";}
  string InstName () const { return "12K";}

  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { if  (StringToUpper(Head.KeyVal("DETECTOR")) == "CFH12K") return new Cfht12K; else return NULL;}
  
  bool GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const;
  Frame SkyRegion(const FitsHeader &Head) const;


  // translators
  SIMPLE_TRANSLATOR(TOADPIXS,"PIXSCAL1");
  virtual RETURN_A_VALUE(TOADRDON,5);
  SIMPLE_TRANSLATOR(TOADUTIM,"UTC-OBS");
  SIMPLE_TRANSLATOR(TOADTYPE,"OBSTYPE");
  SIMPLE_TRANSLATOR(TOADCHIP,"IMAGEID");

  // in principle Illu and overscan defaults are correct

  // Add new key used to compute magnitude of the object

  TRANSLATOR_DEC(TOADPZPT)
  {
    double ccdqer[12] = {86.,77.,72.,95.,85.,90.,71.,76.,76.,79.,71.,84.};
    double ccdqei[12] = {57.,55.,52.,80.,70.,70.,52.,54.,51.,53.,54.,51.};
    double zpccd9r=26.52;
    double zpccd9i=26.22;
    string band = Head.KeyVal("TOADBAND");
    double zerop=0;
    int chip = Head.KeyVal("IMAGEID");
    double AirMass = Head.KeyVal("AIRMASS");
    if (band=="I") 
      {
	double a=0.05;
	zerop=zpccd9i+2.5*log10(ccdqei[chip]/ccdqei[9])-a*AirMass;
      }
    if (band=="R") 
      {
	double a=0.06;
	zerop=zpccd9r+2.5*log10(ccdqer[chip]/ccdqer[9])-a*AirMass;
      }
    double offset= 0;
    if (Head.HasKey("INIZEROP"))
      offset = Head.KeyVal("INIZEROP");
    else
      offset = 2.5*log10(double(Head.KeyVal("TOADEXPO")));
    return FitsKey("TOADPZPT",zerop+offset);
  };

};


bool Cfht12K::GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const
{  
  if (HasLinWCS(Head)) return TanLinWCSFromHeader(Head,Guess);
  else
    {// untested code : not clear that the offsets are correct
      // not clear either that this code is of any use.
      double ra,dec;
      RaDec2000(Head,ra,dec);
      double cosdec = cos(dec*M_PI/180.);
      int chip = Head.KeyVal("TOADCHIP");
      double pixs = Head.KeyVal("TOADPIXS");
      double crpix1 = ((chip%6) - 3) * 0.116 * 3600 /pixs;
      double crpix2 = - (chip/6) * 0.225 * 3600/ pixs;
      ComputeLinWCS(Head, Point(crpix1,crpix2), 
		    GtransfoIdentity(), Guess);
      return true;
    }
}

Frame Cfht12K::SkyRegion(const FitsHeader &Head) const
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

#ifdef OLD_FITSTOAD
//*********************************************************************
FitsKey FitsHeader::Cfht12KFormat(const string &KeyName) const
{
  if (KeyName == "TOADINST") return KeyVal("DETECTOR");
  if (KeyName == "TOADPIXS") return KeyVal("PIXSCAL1");
  if (KeyName == "TOADFILT") return KeyVal("FILTER");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADRDON") return FitsKey(KeyName,5);
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADEQUI") return KeyVal("EQUINOX");
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UTC-OBS");
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
      string keyval = KeyVal("BIASSEC");
      int x0,x1,y0,y1,nx,ny;
      if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
	{
	  nx = x1-x0+1;
	  ny = y1-y0+1;
	  char sec_charstar[64];
	  sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	  return FitsKey(KeyName,string(sec_charstar));
	}
#ifdef STORAGE // there are bad areas in the provided BIASSEC 
      return FitsKey(KeyName,"[2063,18;5,4120]");
#endif
    }
  if (KeyName == "TOADILLU") 
    {
      string keyval = KeyVal("DATASEC");
      int x0,x1,y0,y1,nx,ny ,xmin,ymin;
      // i betcha there's a bug in their DATASEC keyword.
      if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
	{
	    if (x0<x1) 
	       {
	         xmin = x0; 
	         nx = x1-x0+1;
	       }
	    else 
	      {
	       xmin = x1;
	       nx = x0-x1+1;
	      }
	    if (y0<y1) 
	       {
	         ymin = y0; 
	         ny = y1-y0+1;
	       }
	    else 
	      {
	       ymin = y1;
	       ny = y0-y1+1;
	      }
	
	  char sec_charstar[64];
	  sprintf(sec_charstar,"[%d,%d;%d,%d]",xmin,nx,ymin,ny);
	  return FitsKey(KeyName,string(sec_charstar)) ;
	}
    }
  if (KeyName == "TOADTYPE") return KeyVal("OBSTYPE");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return KeyVal("IMAGEID");
  if (KeyName == "TOADBAND")
    {
      string filter = KeyVal("FILTER");
      return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
    }
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);
  return FitsKey(KeyName,NOVAL);      
}
#endif

