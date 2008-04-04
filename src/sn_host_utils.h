#ifndef SN_HOST_UTILS__CC
#define SN_HOST_UTILS__CC
#include <string>
#include <iostream>
#include <fstream>

#include <iomanip> 

#include "dictfile.h"


using namespace std;


// dist ell normalisee
#define DELLN_LIM 1.3
// dist d'assoc entre obj notre cat - cat cfhtls aligne
#define DLIM_CFHTLS 1.5

using namespace std;

void SetSNHostTags(DictFile & lhost);

class ReducedImage ;
void Vignettes(ReducedImage & rim, string dir, DictFile & lhost, int halfsize);
void UpDateFromZPegOut(DictFile & lhost, DictFile & zpeg_out, string suffixe="");


// ============= pour analyse ============= 

int KeepFirstHost(DictFile & file);
int KeepOneSN(DictFile & file, string name);
void EraseStringKeys(DictFile & file);
void WriteShort(DictFile & file, string nom);
void WriteShort(DictFile & file, ofstream & pr , bool oui_append);

void FlagSN(DictFile & file,DictFile &file_flags, DictFile & file_tout);
// trouver la SN dans un eliste avec le nom
bool FindSN(DictFile & file,string name_sn, DictFileIterator & sn_found);



int OK_FLAG(DictFileIterator line);
bool IsWithHost(DictFileIterator line);
// 0 1 ou -1, -100etc.
int Host_Ell_Spir_Type(DictFileIterator line);
// pour les nearby, pour que ce soit coherent
int Host_Ell_Spir_Type(int type) ;
// separe entre type <= Sb et type>Sb, type_cut = 4.
int Host_Earlier_Later_Type(int type, int type_cut);
int Host_red_blue_Type(DictFileIterator line, double BV_cut);
int Host_passive_active_Type(DictFileIterator line, double ssfr_cut);


// 
int PegaseId_3(DictFileIterator line);
int PegaseId_7(DictFileIterator line);
string PegaseId_3String(int id);
string PegaseId_7String(int id);
int PegaseType(DictFileIterator line);
string PegaseString(DictFileIterator line);
string PegaseString(int id);


// ============= cfhtls host============= 

bool Host_is_in_CFHTLS(DictFileIterator line);
bool IsCFHTLSEll(DictFileIterator line);
bool IsCFHTLSSbc(DictFileIterator line);
bool IsCFHTLSScd(DictFileIterator line);
bool IsCFHTLSIrrSB(DictFileIterator line);
int CFHTLS_Id(DictFileIterator line);
string  CFHTLS_String(DictFileIterator line);
int CFHTLSHost_Ell_Spir_Type(DictFileIterator line);

// ==== warnings
bool Warning_Uncertain_Attribution(DictFileIterator line);
bool Warning_Disagreement_w_CFHTLS(DictFileIterator line);
#endif /* SN_HOST_UTILS__CC */
