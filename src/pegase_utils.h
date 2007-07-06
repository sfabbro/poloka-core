#ifndef PEGASE_UTILS__CC
#define PEGASE_UTILS__CC
#include <string>
#include <iostream>
#include <fstream>



#include "dictfile.h"

using namespace std;


// print pegase file header with calibtype  VEGA
void PegaseHeader_(ofstream & prp,string peg_bands, string peg_filters, int nband);
// in case filters are (u) g r i z, and file names are effective_filter_u/g/r/i/z.dat
void PegaseHeader(ofstream & prp, bool with_u, bool with_BV);
// print the beginning of the line
void PrintPegase(ofstream & prp, string sn_name, double in_x, double in_y, double RA, double Dec);
// print the mag part
void PrintPegaseMag(ofstream & prp, double mag , double emag);


// to read the file printed as above
void ReadPegaseInFile(string file_name, DictFile & peg_infile, bool & with_u , bool & with_BV);
void WritePegaseInFile(string file_name, DictFile & peg_infile,  bool with_u , bool with_BV);
bool LookInFile(DictFile & in_file, string sn_name, DictFileIterator &vu);
// to read the out file in a DictFile
void ReadPegaseOutFile(string file_name, DictFile & peg_outfile, bool & with_u ,  bool & with_BV);

// to run zpeg, templates and parameter file are in dir.
void RunZpeg(string file_in, string file_out, string dir, bool is_age_cstr, bool is_z_fixed, bool with_u, bool with_BV);

void Analysis_DZ_over_1plusZ(DictFile & peg_file_out, double & dz_mean, double &  dz_sig, double zmin=0, double zmax=1.5,bool quiet=false);
void Analysis_DZ(DictFile & peg_file_out, double & dz_mean, double &  dz_sig, double zmin=0, double zmax=1.5, bool quiet=false);
// calcule les shifts en dm pour que tmp mag = obs mag, et corrige (corr) ou non
// a appeler apres run zpeg avec z fixe !
void Analysis_DM(DictFile & peg_file_out, double *dm_mean, double *dm_sig, int nband, double zmin=0, double zmax=1.5,bool corr=false);

#endif /* PEGASE_UTILS__CC */
