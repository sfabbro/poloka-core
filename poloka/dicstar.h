#ifndef DICSTAR_SEEN
#define DICSTAR_SEEN

#include <iostream>
#include <fstream>
#include <vector>

#include <poloka/basestar.h>
#include <poloka/starlist.h>

class fastifstream;

/*! \file */
//! Dictionnary Star  
/*! with arbitrary number of (key,value), value=double 
 * 
*/

class DicStar : public BaseStar {
 
  int rank; // position in the file.
  const char *format; // it could be a CountedRef on a string
 
  public :
  
  
    DicStar();
  DicStar(double xx, double yy, double ff); 
  DicStar(const std::vector<string>& firstKeys, const std::vector<string>& newkeys,const char *Format ); 


  bool HasKey(const string &Key) const;
  void AddKey(const string &KeyName, const double &Val);

  unsigned NKeys() const;

  double getval(const string& key) const;
  int setval(const string& key,double val);

  //! position in the file
  int Rank() const { return rank;}

  //! to read once the object is created 
  virtual void    Read(fastifstream& r, const char *Format); 
  
  //! to read and create the object  
  static DicStar* read(const std::vector<string> &firstKeys, const std::vector<string>& newkeys, fastifstream& r, const char *Format); 

  //! this routine should probably not be used
  static DicStar* read(fastifstream& r, const char *Format); // should no be used, in fact

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

class DicStarList : public  StarList<DicStar> {
 public: 

  DicStar *EmptyStar() const;
  DicStarList() {};
  ~DicStarList() {};
  DicStarList(const string &FileName);
  void init(const string &FileName);

  bool HasKey(const string &Key) const ;
 private:
  std::vector<string> key;
  std::vector<string> firstKey;
  //  unsigned baseStarLength;

  int read(fastifstream & r); // should not be used

};

typedef DicStarList::const_iterator DicStarCIterator;
typedef DicStarList::iterator DicStarIterator;

BaseStarList* Dic2Base(DicStarList*);
const BaseStarList* Dic2Base(const DicStarList*);
BaseStarList& Dic2Base(DicStarList&);
const BaseStarList& Dic2Base(const DicStarList&);

#endif
