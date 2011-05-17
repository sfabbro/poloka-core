// This may look like C code, but it is really -*- C++ -*-
#ifndef DAOSTAR__H
#define DAOSTAR__H

#include <vector>
#include "sestar.h"

using namespace std;

class DaoStar : public BaseStar {

 public:

  DaoStar();

  size_t num;
  double eflux, sky, esky, skyskew, sharp, round, chi;
  int flag, iter;
  vector<double> apers, eapers;

  bool IsSaturated() const;
  bool HasNeighbours() const;
  bool IsCosmic() const;
  void FlagAsBlended();
  void FlagAsSaturated();

  // all the i/o crap
  static const char *TypeName() { return "DaoStar"; }
  void read_it(fastifstream & rd, const char *format);
  void writen(ostream &s) const ;
  string WriteHeader_(ostream& s, const char* c=0) const;
};

#include "starlist.h"

class DaoStarList : public StarList<DaoStar> {
 public:
  DaoStarList() {} 
  DaoStarList(const string &FileName) { read(FileName); }
  int write(const string& FileName) const;
};

typedef DaoStarList::const_iterator DaoStarCIterator;
typedef DaoStarList::iterator DaoStarIterator;
typedef CountedRef<DaoStar> DaoStarRef;

const BaseStarList* Dao2Base(const DaoStarList* This);
const BaseStarList& Dao2Base(const DaoStarList& This);

BaseStarList* Dao2Base(DaoStarList* This);
BaseStarList& Dao2Base(DaoStarList& This);

SEStarList* Dao2SE(DaoStarList* This);
SEStarList& Dao2SE(DaoStarList& This);

const SEStarList* Dao2SE(const DaoStarList* This);
const SEStarList& Dao2SE(const DaoStarList& This);

#endif // DAOSTAR__H
