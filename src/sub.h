// This may look like C code, but it is really -*- C++ -*-
#ifndef SUB__H
#define SUB__H

#include "reducedimage.h"
#include "candidatestar.h"
#include "candstar.h"
#include "yquemstar.h"
#include <list>
#include <vector>




class SENearStarList;
class ImageSum;
class ReducedImage;

#include "imagesubtraction.h"


typedef enum StackType { RegularKind = 0, SwarpKind = 1};

class NewStack : public StringList {
public :
  string name;
  string Name() const { return name;}
  ReducedImageRef newStack;
  ImageSubtractionRef sub;
  StackType stackType;

  NewStack(const string &Name = "")
  { name = Name; stackType = RegularKind;}

  // this class will not work if you copy objects after having assigned the
  // pointers newStack and sub.

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

private :
  StringList Ref;
  vector<NewStack> AllNew;
  StringList AllNewNames() const;
  StringList AllInputImages;
  string ImageNameToExtract;
  bool overwrite;
  bool onlyDet;
  bool onlyOneSub;
  bool detectOnAllSub;
  ReducedImageRef RefStack;
  ImageSubtractionRef GlobalSub;
  ReducedImageRef GlobalNew;
  SENearStarList *FakeList;

  //string DatacardsName;
  bool AddFakes;
  bool AssociateGal;
  bool FixRef;
  
  string GeomRefName;
  ReducedImageRef GeometricReference;

public :
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

  /*
  int DoOneSub(const ReducedImage *RefStack, const ImageSum *NewStack, const string &SubName, ImageSubtraction *&Sub);
  

  
  void MatchDetectionsWithFakes(SEStarList *Detections, const string &MatchListName);
  void ApplyCutsAndWrite(ImageSubtraction &ASubtraction);
  void Construction(ImageSubtraction &ASub, YquemStarList & stlcand);
  void ConstructMatchAndCut();
  void Cut_Write(YquemStarList & stlcand, string & cutscanname);
  void Cut_Write(CandStarList & stlcand);
  */
};



#endif /* SUB__H */
