// vltfors.cc
//
// Last change: Time-stamp: <13-May-2002 15:48:19 seb>
//
//

class VltFors : public VirtualInstrument {/* TYPE_SELECTOR */
public:
  static VirtualInstrument *Acceptor(const FitsHeader &Head); // at the end for prototyping
  string TelName() const {return "Vlt";}

  SIMPLE_TRANSLATOR(TOADNAMP,"ESO DET OUTPUTS");
  //SIMPLE_TRANSLATOR(TOADOBJE,"ESO OBS TARG NAME");
  SIMPLE_TRANSLATOR(TOADUTIM,"UTC");
  //SIMPLE_TRANSLATOR(TOADFILT,"FILTER1");

  //SIMPLE_TRANSLATOR(TOADAIRM,"AIRMST");

  TRANSLATOR_DEC(TOADTYPE)
  {
    string type = Head.KeyVal("IMAGETYP");
    RemovePattern(type, " ");    
    if (type != "") return FitsKey("TOADTYPE",type);

    type = Head.KeyVal("DPRTYPE");
    RemovePattern(type, " ");    
    if (type != "") return FitsKey("TOADTYPE",type);
    
    type = Head.KeyVal("HIERARCH ESO DPR TYPE");
    RemovePattern(type, " ");    
    if (type != "") return FitsKey("TOADTYPE",type);
    return FitsKey("TOADTYPE","");
  }

  TRANSLATOR_DEC(TOADAIRM)
  {
    
    string  airmass = Head.KeyVal("AIRMASS");
    RemovePattern(airmass, " ");
    
    if (airmass != "") return FitsKey("TOADAIRM",airmass);

    airmass =  Head.KeyVal("AIRMST");
    
    if (airmass != "") return FitsKey("TOADAIRM",airmass);


    airmass =  Head.KeyVal("HIERARCH ESO TEL AIRM START");
    
    if (airmass != "") return FitsKey("TOADAIRM",airmass);

    return FitsKey("TOADAIRM","");
  }  
  TRANSLATOR_DEC(TOADDECL)
  {
    string decl = DecDegToString(Head.KeyVal("DEC"));
    return FitsKey("TOADDECL",decl);
    
  }

  TRANSLATOR_DEC(TOADRASC)
  {
    string rasc = RaDegToString(Head.KeyVal("RA"));
    return FitsKey("TOADRASC",rasc);
    
  }

  TRANSLATOR_DEC(TOADGAIN)
  {
    int namp = Head.KeyVal("ESO DET OUTPUTS");
    double totgain = 0;
    for (int i=1; i<namp+1; i++)
      {
	char keyname[20];
        sprintf(keyname,"ESO DET OUT%d CONAD",i);
	double gain = Head.KeyVal(keyname);
	totgain += gain;
      }
    if(totgain>0)
       return FitsKey("TOADGAIN",totgain/namp);      

    namp = Head.KeyVal("TOADNAMP");
    totgain = 0;
    for (int i=1; i<namp+1; i++)
      {
	char keyname[20];
        sprintf(keyname,"OUT%dGAIN",i);
	double gain = Head.KeyVal(keyname);
	totgain += gain;
      }
    
    return FitsKey("TOADGAIN",totgain/namp);      
  }

  TRANSLATOR_DEC(TOADRDON)
  {
    int namp = Head.KeyVal("ESO DET OUTPUTS");
    double totron = 0;
    for (int i=1; i<namp; i++)
      {
	char keyname[20];
        sprintf(keyname,"ESO DET OUT%d RON",i);
	double ron = Head.KeyVal(keyname);
	totron += ron;
      }
    if (totron>0)
      return FitsKey("TOADRDON",totron/namp);          
    
    namp = Head.KeyVal("TOADNAMP");
    
    for (int i=1; i<namp+1; i++)
      {
	char keyname[20];
        sprintf(keyname,"OUT%dRON",i);
	double ron = Head.KeyVal(keyname);
	totron += ron;
      }
    return FitsKey("TOADRDON",totron/namp);          
  }

};

class VltFors1 : public VltFors {
public :
  string InstName() const {return "Fors1";}
  string TelInstName() const { return "VltFors1";}

  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const
  {
    char sect[8];
    sprintf(sect,"BSEC%c",'A' + Iamp); 
    string secstr = Head.KeyVal(sect);
    int x0,x1,y0,y1;
    if (sscanf(secstr.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4) return Frame(x0-1, y0-1, x1-1, y1-1);
    cerr << "Overscan region for amp " << char('A'+Iamp) 
	 << " not found for " << Head.FileName() << endl;
    return Frame();
  }
  
  Frame IlluRegion(const FitsHeader &Head, const int Iamp) const
  {      
    switch (Iamp)
      {
      case 1 : return Frame(  16,     0, 1039, 1023);
      case 2 : return Frame(1040,     0, 2063, 1023);
      case 3 : return Frame(  16,  1024, 1039, 2047);
      case 4 : return Frame(1040,  1024, 2063, 2047);
      default: cerr<< "No such amplifier : "<< Iamp << endl; 
      }
    return Frame(0,2047,0,2047);
  }

  Frame TotalIlluRegion(const FitsHeader &Head) const
  {
    string illu = "[17:2063;1:2048]";
    Frame result;
    fits_imregion_to_frame(illu, result);
    return result;
  }


  //SIMPLE_TRANSLATOR(TOADPIXS,"PIXSCALE");
  TRANSLATOR_DEC(TOADPIXS)
  {
    
    string pixs = Head.KeyVal("PIXSCALE");
    RemovePattern(pixs, " ");

    if (pixs != "") return FitsKey("TOADPIXS",pixs);

    pixs =  Head.KeyVal("INSPIXSC");
    
    if (pixs != "") return FitsKey("TOADPIXS",pixs);

    pixs =  Head.KeyVal("HIERARCH ESO INS PIXSCALE");
    
    if (pixs != "") return FitsKey("TOADPIXS",pixs);

    return FitsKey("TOADFILT","");

  }
  TRANSLATOR_DEC(TOADFILT)
  {
    
    string band = Head.KeyVal("FILTER");
    RemovePattern(band, " ");

    if (band != "") return FitsKey("TOADFILT",band);
    band =  Head.KeyVal("FILT1NAM");
	
    if (band != "") return FitsKey("TOADFILT",band);

    band =  Head.KeyVal("FILTER1");
	
    if (band != "") return FitsKey("TOADFILT",band);

    band =  Head.KeyVal("HIERARCH ESO INS FILT1 NAME");

	
    if (band != "") return FitsKey("TOADFILT",band);

    
    return FitsKey("TOADFILT","");



  }

  TRANSLATOR_DEC(TOADOBJE)
  {
    string obje =  Head.KeyVal("TARGNAME");
    RemovePattern(obje, " ");
    if (obje != "") return FitsKey("TOADOBJE",obje); 

    obje =  Head.KeyVal("OBJECT");
    RemovePattern(obje, " ");
    if (obje != "") return FitsKey("TOADOBJE",obje); 

    obje = Head.KeyVal("ESO OBS TARG NAME");
    if (obje != "") return(FitsKey("TOADOBJE",obje));

    obje = Head.KeyVal("OBSNAME");
    if (obje != "") return(FitsKey("TOADOBJE",obje));

    return FitsKey("TOADOBJE","");
    
  }  

  //SIMPLE_TRANSLATOR(TOADOBJE,"OBSNAME");
  
  TRANSLATOR_DEC(TOADNAMP)
  {
    string namp = Head.KeyVal("OUTPUTS");
    RemovePattern(namp, " ");
    if (namp != "") return FitsKey("TOADNAMP",namp);
    namp = Head.KeyVal("DETOUTPU");
    if (namp != "") return FitsKey("TOADNAMP",namp);

    namp = Head.KeyVal("HIERARCH ESO DET OUTPUTS");
    if (namp != "") return FitsKey("TOADNAMP",namp);


    return FitsKey("TOADNAMP","");
  }    
  
};

VirtualInstrument* VltFors::Acceptor(const FitsHeader &Head)
{
  string inst = StringToUpper(string(Head.KeyVal("INSTRUME")));
  RemovePattern(inst, " ");
  if (inst == "FORS1") return new VltFors1;
  return NULL;
}
