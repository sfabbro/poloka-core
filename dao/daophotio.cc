#include "daophotio.h"

int DaoFileNumber(const DbImageCatalogKind FileType)
{
  switch (FileType)
    {      
    case DaophotAls: return 1;
    case DaophotAp:  return 2;
    case DaophotLst: 
    case DaophotNei: return 3;
    default: cerr << " DaoFileNumber(" << FileType << ") Error:  file type not implemented \n";
    }
  return 0;
}

string DaoFileExtension(const DbImageCatalogKind FileType)
{
  switch (FileType)
    {
    case DaophotAls: return "als";
    case DaophotNst: return "nst";
    case DaophotPk : return  "pk";
    case DaophotNei: return "nei";
    case DaophotLst: return "lst";
    case DaophotAp : return  "ap";
    default: cerr << " DaoExtension(" << FileType << ") : Error : file type not implemented\n";
    }
  return "";
}

// partial specialization for daophot aperture file
template<> void read_dao_star<DaophotAp>(istream &daostream, SEStar &star)
{
  double skyerror, skyskew, dmag;
  daostream >> star.Fond() >> skyerror>> skyskew >> dmag;
  star.EFlux() = 0.921034 * dmag * star.flux;
}

// partial specialization for daophot pick file
template<> void read_dao_star<DaophotLst>(istream &daostream, SEStar &star)
{
  daostream >> star.Fond();
}

// partial specialization for daophot neighbor file
template<> void read_dao_star<DaophotNei>(istream &daostream, SEStar &star)
{
  daostream >> star.Fond() >> star.Chi();
  char flag;
  daostream >> flag;
  if (isdigit(flag)) daostream.unget();
  else star.Flag() |= 256; 
}

// partial specialization for daophot  allstar file
template<> void read_dao_star<DaophotAls>(istream &daostream, SEStar &star)
{
  double dmag;      
  daostream >> dmag;
  star.EFlux() = 0.921034 * dmag * fabs(star.flux);
  daostream >> star.Fond();
  double iter;
  daostream >> iter;
  star.Iter() = int(iter);
  daostream >> star.Chi();
  daostream >> star.Sharp();
}

// partial specialization for daophot peak file
template<> void read_dao_star<DaophotPk>(istream &daostream, SEStar &star)
{
  read_dao_star<DaophotAls>(daostream, star);
}

// partial specialization for daophot  nstar file
template<> void read_dao_star<DaophotNst>(istream &daostream, SEStar &star)
{
  read_dao_star<DaophotAls>(daostream, star);
}

// partial specialization for daophot aperture file
template<> void write_dao_star<DaophotAp>(ostream &daostream, const SEStar &star)
{
  // only one aperture is supported! Yet to change.
  double eflux = (star.EFlux() > 0) ? star.EFlux() : 10000;
  double magerr = min(1.0857 * eflux / fabs(star.flux), 9.9999);
  daostream << endl << "    " << setw(9) << star.Fond() 
	    << setprecision(2) << setw(6) << 0.5 << setw(6) << 0.5  // bogus errors and skewness on sky
	    << setprecision(4) << setw(9) << magerr << endl;
}

// partial specialization for daophot peak file
template<> void write_dao_star<DaophotLst>(ostream &daostream, const SEStar &star)
{
  int ndigits=3;
  if (star.Fond() > 9999.999) ndigits=2;
  daostream << setw(9) << setprecision(ndigits) << star.Fond();
}

// partial specialization for daophot neighbor file
template<> void write_dao_star<DaophotNei>(ostream &daostream, const SEStar &star)
{
  char flag = ' ';
  if (star.Flag() != 0) flag = '*';
  daostream << setw(9) << star.Fond()
	    << setw(9) << setprecision(3) << star.Chi() 
	    << "   " << flag;
}

// partial specialization for daophot allstar file
template<> void write_dao_star<DaophotAls>(ostream &daostream, const SEStar &star)
{
  int ndigits;
  if (star.Fond() > 0) ndigits = min(3, max(0, 6-int(log10(star.Fond()))));
  else ndigits = min(3, max(0, 5-int(log10(-star.Fond()+0.001))));
  daostream << setw(9) << setprecision(4) << 1.0857 * star.EFlux() / fabs(star.flux)
	    << setw(9) << setprecision(ndigits) << star.Fond() 
	    << setw(9) << setprecision(0) << star.Iter()
	    << setw(9) << setprecision(3) << star.Chi()
	    << setw(9) << star.Sharp();
}

void read_dao_header(istream &daostream, int& Ncol, int& Nrow, float& LowBad, float& Threshold, 
		     float& Ap1, float& Gain, float& ReadNoise, float& FitRad)
{
  string dummy;
  getline(daostream, dummy);
  // have not seen an utility of nl yet, so do not pass it as argument
  int nl;
  daostream >> nl >> Ncol >> Nrow >> LowBad >> Threshold >> Ap1 >> Gain >> ReadNoise >> FitRad;
}

void read_dao_header(istream &daostream, string& header)
{
  getline(daostream, header);
  string dummy;
  getline(daostream, dummy);
  header += "\n" + dummy;
}
