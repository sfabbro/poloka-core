// $Id: datacards.cc,v 1.3 2004/10/20 10:43:43 astier Exp $
//
// Datacards, acquisition EROS II
//
//
// Eric Aubourg, Decembre 95
//
// DAPNIA/SPP (Saclay) / CEA

#include "datacards.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstdio>

//++
// Class	DataCards
// Lib		Outils++ 
// include	datacards.h
//
//   	Cette classe permet la gestion des parametres d'un programme a partir 
//   	de mot-cle (lecture d'un fichier par exemple)
//--

//++
// Titre	Constructeurs
//--
//++
//
// DataCards()
// DataCards(string const& fn)
//    Createur avec lecture des parametres ds le fichier de nom "fn"
//--

DataCards::DataCards()
{
}

DataCards::DataCards(string const& fn)
{
    ReadFile(fn);
}

//++
// Titre	Methodes
//--
//++
// AddProcF(ProcCard f, string const& mtch="*")
//	Ajoute une fonction de traitement a la liste pour les mots cle
//      compatibles avec la chaine "mtch" ("mtch" peut contenir "*" en debut 
//      fin de mot) 
//
// Clear()
//	Supprime les cartes existantes
//
// ReadFile(string const& fn) 
//      Lit le contenu du fichiers "fn" et ajoute les cartes a la liste
//
// AppendCard(string const& line)
//	Rajoute la carte "line" a la liste
//--

void
DataCards::AddProcF(ProcCard f, string const& mtch)
{
CrdPF mpf;
if (f == NULL)  return;
mpf.pf = f;
if (mtch.length() <= 0)  mpf.patt = "*";
else  mpf.patt = mtch;
cpfs.push_back(mpf);

// On applique cette nouvelle fonction aux cartes existantes
CardList::iterator  ic;
for(ic = cards.begin(); ic != cards.end(); ic ++)
  {
  vector<string>::iterator ik;
  string tks;
  for(ik = (*ic).tokens.begin(); ik != (*ic).tokens.end(); ik++)
    tks = tks + " " + (*ik);
  ApplyPF(mpf, (*ic).kw, tks);
  }
}

void 
DataCards::Clear()
{
cards.erase(cards.begin(), cards.end());
}

void
DataCards::ReadFile(string const& fn)
{
  /* was used in Peida to read default datacards files. 
     senseless in our case */
  DoReadFile(fn);
}

void
DataCards::AppendCard(string const& crd)
{
Card c;
size_t p = 1;
size_t q = crd.find_first_of(" \t");
size_t l = crd.length();

string toks;
if (l < 2)  return;
if (crd[0] != '@')  return;

if (q < l)
  {  c.kw = crd.substr(p,q-p);  toks = crd.substr(q, l-q); }
else { c.kw = crd.substr(p,l-p);  toks = ""; }
//  On applique les ProcFunc's
ApplyPFL(c.kw, toks);
while (q < l) 
  {
  p = crd.find_first_not_of(" \t",q+1); // au debut d'un token
  if (p>=l) break;
  q = crd.find_first_of(" \t",p); // la fin du token;
  string token = crd.substr(p,q-p);
  c.tokens.push_back(token);
  }
// On supprime la carte de la liste, si elle existe deja ...
RemoveCard(c.kw);
cards.push_back(c);
}

void
DataCards::DoReadFile(string const& fn)
{
char line_buff[512];
FILE *fip;

if ( (fip = fopen(fn.c_str(),"r")) == NULL )
  {
    cerr << " DataCards::DoReadFile cannot open file " << fn << endl;
    return;
  }
while (fgets(line_buff,511,fip) != NULL)
  {
    char *last = line_buff+strlen(line_buff) -1;
    if (*last == '\n') *last = '\0'; /* CR at end of line : remove it */
    string line(line_buff);
    AppendCard(line);
  }
 fclose(fip);
}

int
DataCards::ApplyPF(CrdPF & cpf, string const& key, string const& toks)
{
size_t l,lk;
int rc = 0;
// On verifie si le "pattern" correspond
bool mtch = false;
l = cpf.patt.length(); 
if (cpf.patt == "*")  mtch = true;
else if (cpf.patt[0] == '*')   
  {
  lk = key.length();
  if (cpf.patt[l-1] != '*')
    {
    if (strcmp(key.c_str()+(lk-l+1), cpf.patt.c_str()+1) == 0)   mtch = true;
    }
  else if  (key.find(cpf.patt.substr(1,l-2)) < lk)  mtch = true;
  }
else if (cpf.patt[l-1] == '*')
  {
  if ( strncmp(key.c_str(), cpf.patt.c_str(),l-1) == 0)  mtch = true;
  }
else if (key == cpf.patt)  mtch = true;

// Si oui, on appelle la fonction correspondante
if (mtch)  rc = cpf.pf(key, toks); 

return(rc);
}


int
DataCards::ApplyPFL(string const& key, string const& toks)
{
int rc = 0;
CrdPFList::iterator icf;
for(icf = cpfs.begin(); icf != cpfs.end(); icf++)
  rc += ApplyPF((*icf), key, toks);
return(rc);
}

void  
DataCards::RemoveCard(string const& key)
{
CardList::iterator i;
for(i=cards.begin(); i != cards.end(); i++)
  if ((*i).kw == key) { cards.erase(i);  break;  }
}

DataCards::Card *
DataCards::FindKey(string const& key)
{
/*
  CardList::iterator i = find_if(cards.begin(), cards.end(), bind2nd(KeyEq(),key));
  if (i == cards.end() ) return NULL;
*/
  CardList::iterator i;
  for(i=cards.begin(); i != cards.end(); i++)
    if ((*i).kw == key) return &*i;

  return NULL;
}

//++
// Titre	Acces aux parametres
//--
//++
// int   NbCards()
//	Renvoie le nombre de cartes data
// bool	 HasKey(string const& key) 
//	Indique l'existence d'une carte avec la cle "key"
// int   NbParam(string const& key)
//	Indique le nombre de parametre (separes par des espaces) pour la cle "key"
// string  SParam(string const& key, int num = 0, string def="")
//	Renvoie la valeur du parametre numero "num" ( 0..(NbParam()-1) ) sous forme de
//	chaine de caracteres ("string")
// long    IParam(string const& key, int numero = 0, long def = 0)
//	Renvoie la valeur du parametre numero "num" ( 0..(NbParam()-1) ) convertie 
//	en entier ("long")
// double  DParam(string const& key, int numero = 0, double def = 0)
//	Renvoie la valeur du parametre numero "num" ( 0..(NbParam()-1) ) convertie 
//	en flottant ("double")
//--


bool
DataCards::HasKey(string const& key)
{
  return FindKey(key) != NULL;
}

int 
DataCards::NbCards()
{
return(cards.size());
}

int 
DataCards::NbParam(string const& key)
{
  DataCards::Card * p = FindKey(key);
  if (!p) return(-1);
  else return(p->tokens.size());
}

string
DataCards::SParam(string const& key, int numero, string def)
{
  DataCards::Card * p = FindKey(key);
  if (!p) return def;
  if ( (numero < 0) || (numero >= int(p->tokens.size())) )  return def;
  return p->tokens[numero];
}

long
DataCards::IParam(string const& key, int numero, long def)
{
  string p = SParam(key, numero, "");
  if (p == "") return def;
  long i;
  //istrstream(p.c_str(), p.length()) >> i;
  sscanf(p.c_str(),"%ld",&i);
  return i;
}

double
DataCards::DParam(string const& key, int numero, double def)
{
  string p = SParam(key, numero, "");
  if (p == "") return def;
  double i;
  //istrstream(p.c_str(), p.length()) >> i;
  sscanf(p.c_str(),"%lg",&i);
  return i;
}

   
ostream& operator << (ostream& s, DataCards c)
{
  for (DataCards::CardList::iterator i = c.cards.begin(); i != c.cards.end(); i++) {
    s << setw(10) << (*i).kw << " ";
    for (vector<string>::iterator j = (*i).tokens.begin(); j != (*i).tokens.end(); j++)
      s << (*j) << " ";
    s << endl;
  }
  return s;
}



