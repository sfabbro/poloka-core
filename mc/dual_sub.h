// This may look like C code, but it is really -*- C++ -*-
#ifndef DUALSUB__H
#define DUALSUB__H

#include "sub.h"
//#include "mcimage.h"
#include "transformedimage.h"

class SimSNStarList ;


class DualSub : public Sub 
{

  string original_globnewname ;
  string original_globsubname ;
  string globmcresultname ;

  StringList Original_AllNewNames;

  // they will be generated in the coordinates system and photometric system of GeometricReference

  SimSNStarList *SNList;
  string Dual_Image_Path;
public :


  string MCResultName() const ;
  string GlobalMCResultName() const { return globmcresultname; }
  string& GlobalMCResultName()  { return globmcresultname; }
  string Original_GlobNewName() const { return original_globnewname;}
  string & Original_GlobNewName() { return original_globnewname;}
  string Original_GlobSubName() const { return original_globsubname;}
  string & Original_GlobSubName() { return original_globsubname;}
  string MCDir() const { string dir = "./mc/" ; return(dir) ;}

  //! redo the Sub that was done in MC mode.
  //! imc is in case we are looping to execute the MCMode several times, to increase the statistics. Then the directories and result list will contains
  //! the number imc in their names, and a global result file will be produced, //! by appending during the loop all the partial result files.
  DualSub(Sub const & ASub, string image_path, string list_name);
  void DoIt();
  void MakeFakeList(string list_name);
  void MatchDetectionsWithFakes();
  void MakeDetectionsWithFakes(string &ListName);

};

void ProcessDualSub(Sub const & ASub, string image_path, string list_name);




#endif
