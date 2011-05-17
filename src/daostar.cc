#include <iterator>
#include "daostar.h"
#include "fastifstream.h"

static const int DAO_SATUR_FLAG = 1;
static const int DAO_NEIGHBOUR_FLAG = 2;

DaoStar::DaoStar() : 
  BaseStar(0., 0., 0.), num(0), sky(0.), esky(0.), skyskew(0.), sharp(0.5), round(0.), chi(0.), flag(0), iter(0) {}

//! The object has neighbours, bright and close enough to bias aperture photometry
bool DaoStar::HasNeighbours() const {
  return (flag & DAO_NEIGHBOUR_FLAG);
}

bool DaoStar::IsSaturated() const {
  return (flag & DAO_SATUR_FLAG);
}

bool DaoStar::IsCosmic() const {
  return (sharp < 0.1 || sharp > 1.2);
}

void DaoStar::FlagAsBlended() {
  flag |= DAO_NEIGHBOUR_FLAG;
}

void DaoStar::FlagAsSaturated() {
  flag |= DAO_SATUR_FLAG;
}

void DaoStar::writen(ostream& s) const {
  BaseStar::writen(s);
  s << num << ' '
    << sky << ' '
    << esky << ' '
    << skyskew << ' '
    << round << ' '
    << iter << ' '
    << chi << ' '
    << sharp << ' '
    << flag << ' '
    << apers.size() << ' ';
  copy(apers.begin(), apers.end(), ostream_iterator<double>(s, " "));
  copy(eapers.begin(), eapers.end(), ostream_iterator<double>(s, " "));
}

void DaoStar::read_it(fastifstream& r, const char* Format) {
  BaseStar::read_it(r, Format);
  int format = DecodeFormat(Format, "DaoStar");
  r >> num
    >> sky
    >> esky
    >> skyskew
    >> round
    >> iter
    >> chi
    >> sharp
    >> flag;
  size_t naper;
  r >> naper;
  apers.reserve(naper);
  eapers.reserve(naper);
  for (size_t i=0; i<naper; ++i) { r >> apers[i]; }
  for (size_t i=0; i<naper; ++i) { r >> eapers[i]; }
}

string DaoStar::WriteHeader_(ostream& s, const char* c) const {
  string baseStarFormat = BaseStar::WriteHeader_(s, c);
  if (!c) c = "";
  s << "# num" << c << " : \n"
    << "# sky" << c << " : \n"
    << "# esky" << c << " : \n"
    << "# skyskew" << c << " : \n"
    << "# round" << c << " : \n"
    << "# iter" << c << " : \n"
    << "# chi" << c << " : \n"
    << "# sharp" << c << " : \n"
    << "# flag" << c << " : \n"
    << "# naper" << c << " : \n";
  size_t naper = apers.size();
  for (size_t i=0; i<naper; ++i) {
    s << "# aper"  << i << c << " : aperture flux\n"
      << "# eaper" << i << c << " : aperture flux error\n";
  }
  
  return baseStarFormat + " DaoStar 1 ";
}

#include "starlist.cc"
template class StarList<DaoStar>;

BaseStarList* Dao2Base(DaoStarList* List) { 
  return (BaseStarList *) List;
}

const BaseStarList* Dao2Base(const DaoStarList* List) { 
  return (const BaseStarList *) List;
}

BaseStarList& Dao2Base(DaoStarList& List) { 
  return (BaseStarList &) List;
}

const BaseStarList& Dao2Base(const DaoStarList& List) { 
  return (const BaseStarList &) List;
}
