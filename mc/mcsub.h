// This may look like C code, but it is really -*- C++ -*-
#ifndef MCSUB__H
#define MCSUB__H

#include "sub.h"
#include "mcimage.h"
#include "transformedimage.h"

class SimSNStarList ;


class MCSub : public Sub 
{
  AddMethod Addition_Method;
  string original_globnewname ;
  string original_globsubname ;
  string globmcresultname ;

  StringList Original_AllNewNames;

  // they will be generated in the coordinates system and photometric system of GeometricReference
  SimSNStarList *SNList;
public :
  int I_MC ; // if the MC is run in a loop to increase the statistics.
  string FakeListName() const ;
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
  MCSub(Sub const & ASub, const int imc=-1);
  void DoIt();
  void MakeFakeList();
  void MatchDetectionsWithFakes();
private :
  ImageGtransfoRef  GetFirstTransfo(); 
  bool CheckTransfoHomogeneity();

};

void MCProcess(Sub const & ASub);




#endif
