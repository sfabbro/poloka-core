/*
  Attempting to add VLT ISAAC as one instrument.

*/
#ifdef VIRTUAL_INSTRUMENTS
class VltHawk : public VirtualInstrument { /* TYPE_SELECTOR */
public:
  string TelInstName() const {return "VltHawk";};
  string TelName() const {return "Vlt";}
  string InstName() const { return "Hawk";}
  
  static VirtualInstrument *Acceptor(const FitsHeader &Head) {
    if (CheckKeyToUpper(Head,"INSTRUME","HAWKI")) 
      return new VltHawk;
    return NULL;
  }

  SIMPLE_TRANSLATOR(TOADPIXS,"HIERARCH ESO INS PIXSCALE");
  SIMPLE_TRANSLATOR(TOADUTIM,"UTC");

  FitsKey TOADBAND(const FitsHeader &Head, const bool Warn) const {
    string filt1 = Head.KeyVal("FILTER1", true);
    if (filt1.find("OPEN")==string::npos) return FitsKey("TOADBAND",filt1);
    return FitsKey("TOADBAND",string(Head.KeyVal("FILTER2",true)));
  }


  FitsKey TOADFILT(const FitsHeader &Head, const bool Warn) const {
    return TOADBAND(Head,Warn);
  }

  TRANSLATOR_DEC(TOADGAIN) {

    // gains of the four HAWK chips (A courtesy of Chris Lidman):
    /*The gain is the product of two factors

    gain * DET.NDIT

    where DET.DIT is the of reads that are averaged into a single image.

    The gains are different for each chip. They are:

    CHIP1  1.705 
    CHIP2  1.870 
    CHIP3  1.735 
    CHIP4  2.110

    */
    if (Head.HasKey("GAIN")) return FitsKey("GAIN",double(Head.KeyVal("GAIN")));

    static double gains[4] = {1.705,1.87,1.73,2.11};
    int chip = Head.KeyVal("ESO DET CHIP NO");
    if (chip<1 || chip >4)
      {
	cout << "Invalid chip number in image " << Head.FileName() << endl;
	return FitsKey("TOADGAIN",-1);
      }
    double nreads = Head.KeyVal("ESO DET NDIT");
    return FitsKey("TOADGAIN",gains[chip-1]*nreads);
  }

  TRANSLATOR_DEC(TOADRDON) {
    // found on 
    // http://www.eso.org/observing/etc/doc/ut4/hawki/helphawki.html
    return FitsKey("TOADRDON" , 5); 
  }

  /*
    Convert the RA and DEC from degrees, to hours/degrees, minute and
    seconds.
   */
  /*
    There are two values for the air mass. One from the beginning of
    the observation and one for the end of it. The mean is returned.
   */
  /*
    Observation type.
   */
  SIMPLE_TRANSLATOR(TOADTYPE,"HIERARCH ESO DPR TYPE");


  SIMPLE_TRANSLATOR(TOADCHIP,"ESO DET CHIP NO");


  


};
#endif
