// This may look like C code, but it is really -*- C++ -*-
#ifndef SUB__H
#define SUB__H

#include "reducedimage.h"

#include <list>
#include <vector>




class ImageSum;
class ReducedImage;

#include "imagesubtraction.h"


typedef enum StackType { RegularKind = 0, SwarpKind = 1};

class NewStack : public StringList {
public :
  string name;
  string Name() const { return name;}
  ReducedImageRef newStack; // will points in fact on an ImageSum.
  ReducedImageRef original_newStack; 
  ImageSubtractionRef sub; // inherits from the PSFMatch
  ImageSubtractionRef original_sub;
  StackType stackType;
  NewStack(const string &Name = "")
  { name = Name; stackType = RegularKind;}

  // this class will not work if you copy objects after having assigned the
  // pointers newStack and sub.

  //! Name of the subtraction from Name of the stack: if the stackname is new1, then the SubName is new1. If the stackname is blabla, then the SubName is blabla_sub.
  string SubName() const;

  ~NewStack();
};


/* \file */

//! Handling of an actual subtraction (shift-coadd-subtract-detect). See \ref subfile for the way to drive it.
class Sub {
  

  /* there are no pointers to ReducedImage's in this class but
     CountedRef's ro ReducedImage (called ReducedImageRef and co), so
     that destructors (eg ~ImageSum()) are automatically called when
     needed, in practise by ~Sub 
  */

  //private :
protected :
  StringList Ref;
  vector<NewStack> AllNew;
  StringList AllInputImages;
  string ImageNameToExtract;
  bool overwrite;
  bool onlyDet;
  bool onlyOneSub;
  bool detectOnAllSub;
  ReducedImageRef RefStack;
  ReducedImageRef Original_New;
  ImageSubtractionRef GlobalSub;
  ImageSubtractionRef Original_Sub;
  ReducedImageRef GlobalNew; // en fait une ImageSum en pratique.
  string globnewname ;
  string globsubname ;

  //string DatacardsName;
  bool FixRef;
  
  string GeomRefName;
  ReducedImageRef GeometricReference;

public :
  StringList AllNewNames()const ;
  string GlobalNewName() const { return globnewname;}
  string GlobalSubName() const { return globsubname;}
  string& GlobalNewName() { return globnewname;}
  string& GlobalSubName() { return globsubname;}
  //! the constructor. see \ref subfile for the syntax of the "subfile"
  Sub(const string &SubFileName, const bool Overwrite = false, const bool OnlyDet = false);


  ~Sub();
  ReducedImage * ExtractSubimage(const string &SubImageName);
  void RemoveImage (const string & ToRemove);
  void SubstituteName(const string &Original, const string &Substitution);
  int  DoIt();
  int  ExpectedMagLim();
  ReducedImage* DoOneStack(const StringList &Components, 
			   const string &StackName, 
			   const int ToDo,
			   const StackType ST= RegularKind);

  void FindGeometricReference();
  int  CheckNewNames();

  void RunDetection();

 
};



#endif /* SUB__H */
