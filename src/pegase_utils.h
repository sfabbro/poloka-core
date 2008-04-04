#ifndef PEGASE_UTILS__CC
#define PEGASE_UTILS__CC
#include <string>
#include <iostream>
#include <fstream>



#include "dictfile.h"

using namespace std;



string * PegaseKeys(int & N, bool with_u, bool with_BVUR, bool first_sol=false);
void AddPegaseKeys(DictFile & peg_file_out,  string suffixe, bool with_u, bool with_BVUR);
void CheckPegaseKeys(DictFile & peg_file_out,  string suffixe, bool with_u, bool with_BVUR);

// print pegase file header with calibtype  VEGA or AB
void PegaseHeader_(ofstream & prp,string peg_bands, string peg_filters, int nband, bool use_AB);

// in case filters are (u) g r i z, and file names are effective_filter_u/g/r/i/z.dat
void PegaseHeader(ofstream & prp, bool with_u, bool with_BVUR, bool use_AB);
// print the beginning of the line
void PrintPegase(ofstream & prp, string sn_name, double in_x, double in_y, double RA, double Dec);
// print the mag part
void PrintPegaseMag(ofstream & prp, double mag , double emag);


// to read the file printed as above
void ReadPegaseInFile(string file_name, DictFile & peg_infile, bool & with_u , bool & with_BVUR);
void WritePegaseInFile(string file_name, DictFile & peg_infile,  bool with_u , bool with_BVUR, bool is_host_list, bool use_AB=false);
// to read the out file in a DictFile
void ReadPegaseOutFile(string file_name, DictFile & peg_outfile, bool & with_u ,  bool & with_BVUR, bool oui_append = false);

// to run zpeg, templates and parameter file are in dir.
void RunZpeg(string file_in, string file_out, string dir, bool is_age_cstr, bool is_z_fixed, bool with_u, bool with_BVUR);

//=================== pour classification =================== 

bool IsPegaseEll(DictFileIterator line);
bool IsPegaseEll(double tmpf);
bool IsPegaseSa(DictFileIterator line);
bool IsPegaseSb(DictFileIterator line);
bool IsPegaseSbc(DictFileIterator line);
bool IsPegaseSc(DictFileIterator line);
bool IsPegaseSd(DictFileIterator line);
bool IsPegaseSpir(DictFileIterator line);
bool IsPegaseEarlySpir(DictFileIterator line);
bool IsPegaseLateSpir(DictFileIterator line);
bool IsPegaseIrr(DictFileIterator line);







#endif /* PEGASE_UTILS__CC */
