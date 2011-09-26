#ifndef PMSTAR__H
#define PMSTAR__H

#include "basestar.h"
#include "iostream"

class fastifstream;

class PmStar : public BaseStar
{
 public :
  double pmx, pmy; // proper motions, in pixels/day;

 public:
 PmStar() : pmx(0), pmy(0) {};

 PmStar(const BaseStar &B, const double &Pmx=0 , const double &Pmy = 0) : BaseStar(B) , pmx(Pmx) , pmy(Pmy) {};

    //I/O's

  std::string WriteHeader_(std::ostream &s, const char* i) const;
  void writen(std::ostream &s) const;

  static PmStar* read(fastifstream & Rd, const char *Format); 
  void read_it(fastifstream & Rd, const char *Format); 

};
  
#include "starlist.h"

class PmStarList : public StarList<PmStar> 
{

 public :

  double refDate; // the idea is that if refdate is not set, there is no refdate in the output file

  PmStarList() : refDate(-1e30) {};

  PmStarList(const double RefDate) {SetRefDate(RefDate);}

  PmStarList(std::string FileName);

  void read(const std::string &FileName);

  void SetRefDate(const double &RefDate);

};



#endif /* PMSTAR__H */
