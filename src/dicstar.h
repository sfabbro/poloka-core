#ifndef DICSTAR_SEEN
#define DICSTAR_SEEN

#include <iostream>
#include <fstream>
#include <vector>

#include "basestar.h"

/*! \file */
//! Dictionnary Star  
/*! with arbitrary number of (key,value), value=double 
 * 
*/

class DicStar : public BaseStar {
  
  public :
  
  
    DicStar();
  DicStar(double xx, double yy, double ff); 
  DicStar(const std::vector<string>& firstKeys, const std::vector<string>& newkeys); 

  void AddKey(const string &KeyName, const double &Val);

  double getval(const string& key) const;
  int setval(const string& key,double val);

  //! to read once the object is created 
  virtual void    Read(istream& r, const char *Format); 
  
  //! to read and create the object  
  static DicStar* read(const std::vector<string> &firstKeys, const std::vector<string>& newkeys, istream& r, const char *Format); 

  static DicStar* read( istream& r, const char *Format); 
  

  //!  to write the StarList header with the string appended to every ntuple variable (with no end)
  string WriteHeader_(ostream & pr = cout, const char *i=NULL) const;
 
  //! for dump with NO end-of-line
  virtual void    dumpn(ostream& s = cout) const;
  
  //! for dump  
  virtual void    dump(ostream& s = cout) const ;
  
  //! for write with NO end-of-line
  virtual void    writen(ostream& s = cout) const ;
  
  //! for write 
  virtual void    write(ostream& s= cout)  const ;
  
  private :
    
    void Set_to_Zero();

  protected :
   
     std::vector<string> firstkeys;
    std::vector<string> key;
   std::vector<double> val;


   friend class DicStarList;

#ifdef USE_ROOT
  ClassDef(DicStar,1) // Root's stuff
#endif 
    };
    
    
/********************   FIN DEFINITION DicStar   **********************/


/* what concerns the DicStarList's : */
#include "starlist.h"


class DicStarList : public  StarList<DicStar> {
 public: 

  DicStar *EmptyStar() const;
  DicStarList() {};
  ~DicStarList() {};
  DicStarList(const string &FileName);

 private:
  std::vector<string> key;
  std::vector<string> firstKey;
};

typedef DicStarList::const_iterator DicStarCIterator;
typedef DicStarList::iterator DicStarIterator;



#endif
