#ifndef  YQUEMSTAR__H 
#define  YQUEMSTAR__H

#include "persistence.h"

#include "candidatestar.h"
#include "datdetec.h"


//********************   DEFINITION  YquemStar   *********************

// CandidateStar + objet associe sur la reference


class YquemStar : public CandidateStar {
  CLASS_VERSION(YquemStar,1);
  #define YquemStar__is__persistent

public:
SEStar StarRef;
double dass_ref ;

public:
  YquemStar(){dass_ref=0;};
  YquemStar(CandidateStar const & sestar);
  


  int write_nice(ostream & pr) const;
  
  
  /* DOCF for dump with NO end-of-line*/
  virtual void    dumpn(ostream& s = cout) const;
  /* DOCF for dump  */
  virtual void    dump(ostream& s = cout) const ;
  
  /* DOCF for write with NO end-of-line*/
  virtual void    writen(ostream& s = cout) const ;
  
  /* DOCF for write with NO end-of-line*/
  virtual void    writen_scan(ostream& s = cout) const ;
  
  /* DOCF to read once the object is created */
  virtual void    read_it(istream& r, const char *Format); 
  
  /* DOCF to read and create the object  */
  
  static YquemStar* read(istream& r, const char* Format); 
  
  /* DOCF  to write the StarList header with the string ${}^{\star}i$
     appended to every ntuple variable (with no end)  */
  static const char * WriteHeader_(const char *i, ostream & pr = cout);
  
  
  /* DOCF  to write the StarList header (with no end)  */
  
  std::string WriteHeader_(ostream & pr = cout, const char *i = NULL) const;
  
  const char *WriteHeader_scan(ostream & pr = cout , const char*i = NULL ) const;

  /* DOCF to write the list header in a usable form for the scanning program */
    static void WriteHeader_Scan(ostream & pr) ;

  /* DOCF to write the list in a usable form for the scanning program */
  void write_scan(ostream & pr) const;  

  static const char *TypeName() { return "YquemStar";}


};

/* what concerns the YquemStarList's : */
#include "starlist.h"


class YquemStarList : public StarList<YquemStar> {
  public :
    
  YquemStarList(){}
  YquemStarList(const string FileName); 
  //  ~YquemStarList(){cout << " ~YquemStarList : " << this << endl;};
  
  // to write the list for the scaning from the complete list  
  void write_scan(const string &FileName);
 
 void write_nice(const string &FileName); 
 
 
 void ComputeFAperRef(Image & imgref, double rayonref, double fd_ref);
 
 // Construction of the list: add the ref object to the candidatestar
 void Construct_List(CandidateStarList & stl, SEStarList & stlref);
 
 // Apply the cut on the YquemStarList : sigtonoise
 void Cut(YquemStarList & stlcut, DatDetec & dat);
 
 void Construct(CandidateStarList & stl, SEStarList & stlref, 
		Image & imgref, double rayon_ref, double fd_ref);
};


//typedef StarList<YquemStar> YquemStarList;


typedef YquemStarList::const_iterator YquemStarCIterator;
typedef YquemStarList::iterator YquemStarIterator;

BaseStarList* Yquem2Base(YquemStarList * This);
const BaseStarList* Yquem2Base(const YquemStarList * This);



#endif

