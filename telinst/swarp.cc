

#include "wcsutils.h"
#include "astroutils.h" /* for Ra and Dec string <-> deg convertions */

class Swarp: public VirtualInstrument {  /* TYPE_SELECTOR */

public :
  string TelInstName() const {return "Swarp";}
  string TelName () const {return "no telescope";}
  string InstName () const { return "no instrument";}

  static VirtualInstrument *Acceptor(const FitsHeader &Head)
  { if  (string(Head.KeyVal("SOFTNAME")) == "SWarp") return new Swarp; 
  else return NULL;}
  
  bool GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const;


  // translators
  TRANSLATOR_DEC(TOADPIXS)
  {
    if (Head.HasKey("CDELT2", true))
      {
	double cdelt2 = Head.KeyVal("CDELT2");
	return FitsKey("TOADPIXS",fabs(cdelt2*3600.));
      }
    else return FitsKey("TOADPIXS", 0.);
  }

  // Add new key used to compute magnitude of the object

  SIMPLE_TRANSLATOR(TOADPZPT,"ZEROPT");

  TRANSLATOR_DEC(TOADRASC)
  {
    if (Head.HasKey("CRVAL1"))
      {
	double radeg = Head.KeyVal("CRVAL1");
	return FitsKey("TOADRASC",RaDegToString(radeg));
      }
    else return VirtualInstrument::TOADRASC(Head,Warn);
  }


  TRANSLATOR_DEC(TOADDECL)
  {
    if (Head.HasKey("CRVAL2"))
      {
	double radeg = Head.KeyVal("CRVAL2");
	return FitsKey("TOADDECL",DecDegToString(radeg));
      }
    else return VirtualInstrument::TOADDECL(Head,Warn);
  }

};


bool Swarp::GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const
{  
  if (HasLinWCS(Head)) return TanLinWCSFromHeader(Head,Guess);
  else return false;
}

