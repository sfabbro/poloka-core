
string toads_date(const string &FitsDate)
{
  int yy,mm,dd;
  string result="";
  if ((sscanf(FitsDate.c_str(),"%d/%d/%d",&yy,&mm,&dd) != 3) &&
      (sscanf(FitsDate.c_str(),"%d-%d-%d",&yy,&mm,&dd) != 3)
      )
    return result;
  bool ok = false;
  if (dd > 31) {swap (yy,dd); ok = true;}// done
  else if (yy > 31) ok = true;
  if (!ok) return result;
  // convert 2 digits dates
  if (yy < 1949) yy += (yy>50)? 1900 : 2000;
  char date_string[16];
  sprintf(date_string,"%02d/%02d/%04d",dd,mm,yy);
  return string(date_string);
}
typedef map<string,char> StringCharMap;


string ToadBand(const string &keyval) 
{
  // This routine creates (only once) a map of filters to simple bands.
  static bool called = false;
  static StringCharMap ToadBandMap;
  if (!called)
    {
      called = true;
      string X = "UBVRIJHKLGYZ";
      char F;
      string sF;
      for (unsigned int i=0;i<X.length();i++)
	{
	  F = X[i];
	  sF = string(&F,1);
	  ToadBandMap["BESS_"+sF] = F;
	  ToadBandMap["BESS"+sF] = F;
	  ToadBandMap[sF+"_BESS"] = F;
	  ToadBandMap[sF+"BESS"] = F;
          ToadBandMap["guessed_"+sF] = F;
	  ToadBandMap[sF+"_guessed"] = F;
	  ToadBandMap[sF+"guessed"] = F;
	  ToadBandMap["HARRIS_"+sF] = F;
	  ToadBandMap["HARRIS"+sF] = F;
	  ToadBandMap[sF+"_HARRIS"] = F;
	  ToadBandMap[sF+"HARRIS"] = F;
	  ToadBandMap[sF+"_SLOAN"] = F;
	  ToadBandMap[sF+"SLOAN"] = F;
	  ToadBandMap["SLOAN"+sF] = F;
	  ToadBandMap["SLOAN_"+sF] = F;
	  ToadBandMap[sF+"_STROMGREN"] = F;
	  ToadBandMap[sF+"STROMGREN"] = F;
	  ToadBandMap["STROMGREN"+sF] = F;
	  ToadBandMap["STROMGREN_"+sF] = F;
	  ToadBandMap[sF+"CTIO"] = F;
	  ToadBandMap[sF+"_GUNN"] = F;
	  ToadBandMap[sF+"GUNN"] = F;
	  ToadBandMap["GUNN_"+sF] = F;
	  ToadBandMap["GUNN"+sF] = F;
	  ToadBandMap[sF+"_WIDE"] = F;
	  ToadBandMap[sF+"WIDE"] = F;
	  ToadBandMap["WIDE_"+sF] = F;
	  ToadBandMap["WIDE"+sF] = F;
	  ToadBandMap[sF] = F;
	  ToadBandMap["F850LP"] = 'Z';
	  ToadBandMap["F814W"] = 'I';
	  ToadBandMap["F625W"] = 'R';
	  ToadBandMap["F675W"] = 'R';
	  ToadBandMap["F555W"] = 'V';
	  ToadBandMap["F110W"] = 'J';
	  ToadBandMap[sF+"#810"] = F;
	  ToadBandMap[sF+"#811"] = F;
	  ToadBandMap[sF+"#812"] = F;
	  ToadBandMap[sF+"#813"] = F;
	  ToadBandMap[sF+"#814"] = F;
	  ToadBandMap[sF+"#815"] = F;
	  ToadBandMap["SLOGUN"+sF] = F;
	  ToadBandMap["6 "+sF] = F;
	  ToadBandMap["4 "+sF] = F;
	  ToadBandMap[sF+"s"] = F;
	  ToadBandMap["1"+sF] = F;
	  ToadBandMap[sF+"1"] = F;
	  ToadBandMap[sF+"3"] = F;
	  ToadBandMap["1"+sF+"#7"] = F;
	  ToadBandMap["2"+sF+"#74"] = F;
	  ToadBandMap["3"+sF+"#75"] = F;
	  ToadBandMap["4"+sF+"#76"] = F;
	  ToadBandMap["5"+sF+"#12"] = F;
	  ToadBandMap["Sloang`"] = 'G';
	  ToadBandMap["Sloanr`"] = 'R';
	  ToadBandMap["Sloani`"] = 'I';
	  ToadBandMap["Sloanu`"] = 'U';
	  ToadBandMap["Sloanz`"] = 'Z';
	  ToadBandMap[sF+"Johnson"] = F;
	  ToadBandMap["g.MP9401"] = 'G';
	  ToadBandMap["r.MP9601"] = 'R';
	  ToadBandMap["i.MP9701"] = 'I';
	  ToadBandMap["z.MP9801"] = 'Z';
	  ToadBandMap["u.MP9301"] = 'U';
	  ToadBandMap["macho_v"] = 'V';
	  ToadBandMap["macho_r"] = 'R';
	}  
    }
  StringCharMap::iterator match = ToadBandMap.find(keyval);
  if (match == ToadBandMap.end()) match = ToadBandMap.find(StringToUpper(keyval));
  if (match == ToadBandMap.end())
    {
      cerr << " Did not find the associate band with filter " << keyval << endl;
      return keyval;
    }    
  char result = (*match).second; /*ToadBandMap[keyval];*/
  return string(&result,1);
}

Frame VirtualInstrument::AmpRegion(const FitsHeader &Head, const int Iamp) const
{
  Frame illuRegion = IlluRegion(Head, Iamp);
  Frame totIlluRegion = TotalIlluRegion(Head);
  // do not use gtransfo here
  illuRegion.xMin -= totIlluRegion.xMin;
  illuRegion.xMax -= totIlluRegion.xMin;
  illuRegion.yMin -= totIlluRegion.yMin;
  illuRegion.yMax -= totIlluRegion.yMin;
  return illuRegion;
  //return illuRegion.ApplyTransfo(GtransfoLinShift(-totIlluRegion.xMin, - totIlluRegion.yMin));
}

Frame VirtualInstrument::IlluRegion(const FitsHeader &Head, const int Iamp) const
{
  if (Iamp != 1)
    {
      cerr << "using default IlluRegion for iamp != 1. expect serious troubles " << endl;
      cerr << Head.FileName() << ' ' << "tel/inst" << TelInstName() << endl;
    }
  return TotalIlluRegion(Head);
}

Frame VirtualInstrument::TotalIlluRegion(const FitsHeader &Head) const
{
  Frame result;
  if (Head.HasKey("DATASEC"))
    fits_imregion_to_frame(string(Head.KeyVal("DATASEC")), result);
  else result = Frame(Head,WholeSizeFrame);
  return result;
}

Frame VirtualInstrument::OverscanRegion(const FitsHeader &Head, const int Iamp) const
{
  if (Head.HasKey("BIASSEC"))
    {
      string keyval = Head.KeyVal("BIASSEC");
      Frame result;
      fits_imregion_to_frame(keyval, result);
      return result;
    }
  else return Frame();
}

double VirtualInstrument::AmpGain(const FitsHeader &Head, const int Iamp) const
{
  double gain = Head.KeyVal("TOADGAIN");
  return gain;
}

void VirtualInstrumentDestructor(VirtualInstrument *p)
{
  delete p;
}

/******************** implementation of Unknown intrument *************/

class Unknown : public VirtualInstrument {

public :
  static VirtualInstrument *Acceptor(const FitsHeader &Head) {return new Unknown();};
  string TelInstName() const {return "Unknown";};
  string InstName() const {return "Unknown";};
  string TelName() const {return "Unknown";};

};

