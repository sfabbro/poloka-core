// This may look like C code, but it is really -*- C++ -*-


#ifndef  CANDSTAR__H 
#define  CANDSTAR__H

#include "datdetec.h"
#include "candidatestar.h"
#include "sestar.h"

#include "persistence.h"

class CandStar {
  CLASS_VERSION(CandStar,1);
  #define CandStar__is__persistent
public:
int numero;

CandidateStar star;
CandidateStar star1;
CandidateStar star2;

SEStar star_ref;

double dass_ref;
double dass_1;
double dass_2;
double dass_12;



public :
  
  CandStar(){numero=0; dass_1 = -1 ; dass_2 = -1 ; dass_12 = -1 ; 
  dass_ref = -1;};
  
  bool  IsGood(DatDetec & dat);
  
  void  ComputeFaperRef(Image & imgref, double rayon,double  fd_ref);
  
  int write(ostream & pr) const;
  int write_nice(ostream & pr) const;
  void write_scan(ostream & pr) const;
  void WriteHeader(ostream &pr);
  static void WriteHeader_Scan(ostream &pr);
  
  void read_it(istream& r, const char *Format);
  
private:
};


class CandStarList : public list<CandStar> {
  CLASS_VERSION(CandStarList,1);
  #define CandStarList__is__persistent
public :
  CandStarList() {};
  CandStarList(const string &FileName) {read(FileName);};
  int write(const string &FileName) const ;
  int write_nice(const string &FileName) const ;
  int write_scan(const string &FileName) const ;
  int read (const string &FileName);
  void ComputeFAperRef(Image & imgref, double rayonref, double fd_ref);
  void Construct_List(CandidateStarList & stl,  CandidateStarList & stl1,  
		      CandidateStarList & stl2, SEStarList & stlref, double dist_min);
  void Construct_Simple(CandidateStarList & stl,  CandidateStarList & stl1,  
			CandidateStarList & stl2, SEStarList & stlref, Image & imgref,
			double dass_min, double rayon_ref, double fd_ref);  
  void Cut(CandStarList & stlcut, DatDetec & dat);
};


typedef list<CandStar>::const_iterator CandStarCIterator;
typedef list<CandStar>::iterator CandStarIterator;


//void
//Construct_List(CandStarList & stlcand, CandidateStarList & stl,  CandidateStarList & stl1,  
//	       CandidateStarList & stl2, SEStarList & stlref, double dist_min);

//void
//Construct(CandStarList & stlcand, CandidateStarList & stl,  CandidateStarList & stl1,  
//	  CandidateStarList & stl2, SEStarList & stlref, DatDetec & dat,
//	  Image & imgref, Image & imgnew1, Image & imgnew2, Image & imgnew);

//void
//Construct_Simple(CandStarList & stlcand, CandidateStarList & stl,  CandidateStarList & stl1,  
//		 CandidateStarList & stl2, SEStarList & stlref, Image & imgref,
//		 double dass_min, double rayon_ref, double fd_ref);


//void
//Cut(CandStarList & stlcand, CandStarList & stlcut, DatDetec & dat);



//void
//ComputeFAperRef(CandStarList & stlcand,  
//		Image & imgref, double rayonref, double fd_ref);

//CandStar*
//FindClosest(CandStarList & stl, double x, double y, double & dass);


#endif


