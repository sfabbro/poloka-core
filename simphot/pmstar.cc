#include "pmstar.h"
#include "starlistexception.h"
#include "fastifstream.h"

std::string PmStar::WriteHeader_(std::ostream &s, const char* i) const
{
  if (i== NULL) i= "";
  string baseStarFormat =  BaseStar::WriteHeader_(s, i);
  s << "# pmx"<< i <<" : proper motion along x " << endl;
  s << "# pmy"<< i <<" : proper motion along y " << endl;
  return baseStarFormat+" PmStar 1";
}

void PmStar::writen(std::ostream &s) const
{
    BaseStar::writen(s);
    s << pmx << ' ' << pmy << ' ';
}

PmStar* PmStar::read(fastifstream & Rd, const char *Format)
{
    PmStar *s= new PmStar();
    s->read_it(Rd,Format);
    return s;
}

void PmStar::read_it(fastifstream & Rd, const char *Format)
{
  int formatValue = 0;
  if (Format) 
    formatValue = DecodeFormat(Format,"PmStar");
  if (formatValue == 1)
    {
      BaseStar::read_it(Rd , Format);
      Rd >> pmx >> pmy;
    }
 else throw(StarListException(" Unknown format value for PmStar "));
}

#include "starlist.cc"
template class StarList<PmStar>;


PmStarList::PmStarList(std::string FileName)
{
  read(FileName);
}

void PmStarList::read(const std::string &FileName)
{
  StarList<PmStar>::read(FileName);
  const GlobalVal &glob = GlobVal();
  if (!glob.HasKey("REFDATE"))
    throw(PolokaException("PmStarList::PmStarList : "+FileName+ " misses a REFDATE global key "));
  refDate = glob.getDoubleValue("REFDATE");
}

void PmStarList::SetRefDate(const double &RefDate)
{
  refDate = RefDate;
  GlobVal().setDoubleValue("REFDATE",refDate);
}

#ifdef STORAGE
void PmStarList::ChangeRefDate(const double &NewRefDate)
{
  double delta = NewRefDate-refDate;
  for(PmStarList::iterator i= begin(); i != end(); ++i)
    {
      PmStar &s = **i;
      s.x += s.pmx*delta;
      s.y += s.pmy*delta;
    }
  refDate = NewRefDate;
  GlobVal().setDoubleValue("REFDATE",NewRefDate);
}
#endif
