#ifndef ADDLCFAKES__H
#define ADDLCFAKES__H

#include <string>
#include <list>

/* tools to add fake SN to already aligned images, 
   as the differential photometry  processing produces them.
*/



#include "basestar.h"

using namespace std;

struct Observation
{
  string inputName;
  double photFactor;
  double mmjd;  
  Observation(const string &S, const double &pf);
};


struct  SNObservations : public list<Observation>
{
  BaseStar modelStar;
  int xShift, yShift;

  
  SNObservations(const BaseStar &B, int Dx, int Dy) :
       modelStar(B), xShift(Dx), yShift(Dy) {};
  
  void AddObservation(const string &ImageName, const double PhotFactor); 


  void GetMinMaxDates(double &mindate, double &maxdate) const;
       

};


typedef SNObservations::const_iterator ObsCIterator;

class SNList : public list<SNObservations>
{
};

typedef SNList::const_iterator SNCIterator;


class SNAdder
{
 private:
  SNList snList;

 public:
  void AddObject(const SNObservations &LightCurve);
  void GenerateImages() const;
  bool GenerateLightCurveFile(const std::string& FileName, 
			      const string &PhotRef) const;

};



#endif
