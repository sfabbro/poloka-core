#include "daophotio.h"

int DaoFileNumber(const DbImageCatalogKind FileType)
{
  switch (FileType)
    {      
    case DaophotAls: return 1;
    case DaophotAp:  return 2;
    case DaophotLst: 
    case DaophotNei: return 3;
    default: { cerr << " DaoFileNumber : file type " << FileType << " not found or not yet implemented\n";}
    }
  return 0;
}

template<> void read_dao_star<DaophotAp>(istream &daostream, SEStar &star)
{
  double skyerror, skyskew, dmag;
  daostream >> star.Fond() >> skyerror>> skyskew >> dmag;
  star.EFlux() = 0.921034 * dmag * star.flux;
}

template<> void read_dao_star<DaophotLst>(istream &daostream, SEStar &star)
{
  daostream >> star.Fond();
}

template<> void read_dao_star<DaophotNei>(istream &daostream, SEStar &star)
{
  daostream >> star.Fond() >> star.Chi();
  char flag;
  daostream >> flag;
  if (isdigit(flag)) daostream.unget();
  else star.Flag() |= 256; 
}

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

template<> void read_dao_star<DaophotPk>(istream &daostream, SEStar &star)
{
  read_dao_star<DaophotAls>(daostream, star);
}

template<> void read_dao_star<DaophotNst>(istream &daostream, SEStar &star)
{
  read_dao_star<DaophotAls>(daostream, star);
}

template<> void write_dao_star<DaophotAp>(ostream &daostream, const SEStar &star)
{
  // only one aperture is supported! Yet to change.
  double eflux = (star.EFlux() > 0) ? star.EFlux() : 10000;
  double magerr = min(1.0857 * eflux / fabs(star.flux), 9.9999);
  daostream << endl << "    " << setw(9) << star.Fond() 
	    << setprecision(2) << setw(6) << 0.5 << setw(6) << 0.5  // bogus errors and skewness on sky
	    << setprecision(4) << setw(9) << magerr << endl;
}

template<> void write_dao_star<DaophotLst>(ostream &daostream, const SEStar &star)
{
  int ndigits=3;
  if (star.Fond() > 9999.999) ndigits=2;
  daostream << setw(9) << setprecision(ndigits) << star.Fond();
}

template<> void write_dao_star<DaophotNei>(ostream &daostream, const SEStar &star)
{
  char flag = ' ';
  if (star.Flag() != 0) flag = '*';
  daostream << setw(9) << star.Fond()
	    << "   " << flag;
}

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
