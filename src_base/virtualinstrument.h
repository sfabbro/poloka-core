// go through these simple "routines" in case we want to change

#define TRANSLATOR_DEC(RoutineName) \
      FitsKey RoutineName( const FitsHeader &Head, const bool Warn) const

#define SIMPLE_TRANSLATOR(RoutineName,KeyTag) \
       TRANSLATOR_DEC(RoutineName) { return Head.KeyVal(KeyTag,Warn);}

#define SIMPLE_TRANSLATOR_WITHDEF(RoutineName,KeyTag,DefaultValue) \
       TRANSLATOR_DEC(RoutineName) \
      { if (Head.HasKey(KeyTag)) return Head.KeyVal(KeyTag,Warn);\
        else return FitsKey(KeyTag,DefaultValue); }

#define RETURN_A_VALUE(RoutineName,Value) \
        TRANSLATOR_DEC(RoutineName) { return FitsKey(#RoutineName,Value);}



string toads_date(const string &FitsDate);
string ToadBand(const string &keyval);


// This is the class of which to derive the actual instruments.
class VirtualInstrument
{
public :
  // the routine that says if a given header belongs to this class
  // not in the base class
  //  static VirtualInstrument *Acceptor(const FitsHeader &);

  //! -
  virtual string TelInstName() const = 0;
  virtual string TelName() const = 0;
  virtual string InstName() const = 0;

  
  //! default uses BIASSEC
  virtual Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const;

  //! default : uses DATASEC. To be overwritten if not applicable
  virtual Frame TotalIlluRegion(const FitsHeader &Head) const;


  //! region corresponding to amplifier Iamp (only one amp -> iamp = 1)
  virtual Frame IlluRegion(const FitsHeader &Head, const int Iamp) const;


  //! region corresponding to amplifier Iamp once image trimed (only one amp -> iamp = 1)
  virtual Frame AmpRegion(const FitsHeader &Head, const int Iamp) const;

  
  //! gain corresponding to amplifier Iamp (only ona amp -> iamp = 1)
  virtual double AmpGain(const FitsHeader &Head, const int Iamp) const;

  //  virtual Frame GuessSkyArea(const FitsHeader &, const int IAmp) const = 0;

  // the default translators
  // if you change here, change also the documented default in ToadsKeys array...
  // I would advise NOT to change these defaults: derived classes rely on it
  virtual SIMPLE_TRANSLATOR_WITHDEF(TOADPIXS,"PIXSCALE",1.);
  virtual TRANSLATOR_DEC(TOADINST) { return FitsKey("TOADINST",TelInstName());}
  virtual TRANSLATOR_DEC(TOADTELE) { return FitsKey("TOADTELE",TelName());}
  virtual SIMPLE_TRANSLATOR(TOADFILT,"FILTER");
  virtual SIMPLE_TRANSLATOR(TOADEXPO,"EXPTIME");
  virtual SIMPLE_TRANSLATOR_WITHDEF(TOADGAIN,"GAIN",1.);
  virtual SIMPLE_TRANSLATOR_WITHDEF(TOADRDON,"RDNOISE",0.);
  virtual SIMPLE_TRANSLATOR(TOADRASC,"RA");
  virtual SIMPLE_TRANSLATOR(TOADDECL,"DEC");
  virtual SIMPLE_TRANSLATOR(TOADEQUI,"EQUINOX");
  virtual SIMPLE_TRANSLATOR(TOADAIRM,"AIRMASS");
  virtual SIMPLE_TRANSLATOR(TOADUTIM,"UT")
  virtual FitsKey   TOADDATE(const FitsHeader &Head, const bool Warn) const
  { string result = toads_date(string(Head.KeyVal("DATE-OBS")));
    if (result.length() == 0) result = toads_date(string(Head.KeyVal("DATE")));
    if (result.length() == 0) result = toads_date(string(Head.KeyVal("DATEOBS")));
    if (Warn && result.length() == 0) cout << " no DATE-OBS nor DATE nor DATEOBS " << endl;
    return FitsKey("TOADDATE", result);
  }

  virtual  FitsKey   TOADMJD(const FitsHeader &Head, const bool Warn) const
  {
    double val = ModJulianDay(Head);
    return FitsKey("MJD-OBS", val);
  }

  virtual FitsKey   TOADMMJD(const FitsHeader &Head, const bool Warn) const
  { 
    double val = ModifiedModifiedJulianDay(Head);
    return FitsKey("TOADMMJD", val);
  }
  
  virtual FitsKey   TOADBAND(const FitsHeader &Head, const bool Warn) const 
     { 
	 string filt_str = Head.KeyVal("TOADFILT");
	 RemovePattern(filt_str," ");	 
	 return FitsKey("TOADBAND",ToadBand(filt_str));
     }
  virtual SIMPLE_TRANSLATOR(TOADTYPE,"IMAGETYP");
  virtual SIMPLE_TRANSLATOR(TOADOBJE,"OBJECT"); 
  virtual RETURN_A_VALUE(TOADNAMP,1);
  virtual RETURN_A_VALUE(TOADCHIP,1);
  virtual SIMPLE_TRANSLATOR(TOADPZPT,"ZEROUSNO");
  virtual ~VirtualInstrument() {};

  /* if you add things in the functionnalities, update the TestTelInst function
  accordingly */

};
