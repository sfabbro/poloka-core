// This may look like C code, but it is really -*- C++ -*-
//
// $Id: datacards.h,v 1.1 2004/02/20 10:48:36 nrl Exp $
//
// Datacards, acquisition EROS II
//
//
// Eric Aubourg, Decembre 95
// Reza Ansari, Aout 96
//
// DAPNIA/SPP (Saclay) / CEA    LAL - IN2P3/CNRS  (Orsay)

#ifndef DATACARDS_SEEN
#define DATACARDS_SEEN

#include <string>
#include <functional>
#include <list>
#include <vector>
#include <string>

using namespace std;


typedef int (*ProcCard)(string const& key, string const& toks);
 
class DataCards  {
public:
   DataCards();                  
   DataCards(string const& fn);

           // nom dans variable d'environnement PEIDA_DATACARDS
           // par defaut, peida.datacards                        
           // Si pas chemin complet, on tente dans repertoire 
           // en cours, puis dans PEIDA_WORK

   virtual ~DataCards() {}

   void    AddProcF(ProcCard f, string const& mtch="*");

   void    Clear();
   void    ReadFile(string const& fn);
   void    AppendCard(string const& line);

   int     NbCards();
   bool    HasKey(string const& key);
   int     NbParam(string const& key);
   string  SParam(string const& key, int numero = 0, string def="");
   long    IParam(string const& key, int numero = 0, long def = 0);
   double  DParam(string const& key, int numero = 0, double def = 0);
   
   friend ostream& operator << (ostream& s, DataCards c);

public:
   struct Card {
     string kw;
     vector<string> tokens;
     //	 STRUCTCOMPF(Card,kw)
   };
   typedef list<Card> CardList;
   struct CrdPF {
     ProcCard pf;
     string  patt;
     // STRUCTCOMPF(CrdPF,pf)
   };
   typedef list<CrdPF> CrdPFList;
protected:
   CardList cards;
   CrdPFList cpfs;

   void  DoReadFile(string const& fn);

   int   ApplyPF(CrdPF & cpf, string const& key, string const& toks);
   int   ApplyPFL(string const& key, string const& toks);

   void  RemoveCard(string const& key);
   
   Card* FindKey(string const& key);

#ifndef SWIG  
  struct KeyEq : binary_function<Card, string, bool> {
    bool operator()(const Card& x, const string& y) const 
    { return x.kw == y; }
   };
#endif

};
#endif
