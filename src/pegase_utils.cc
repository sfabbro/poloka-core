#include <iomanip> 

#include "pegase_utils.h"
#include "fileutils.h"
#include "vutils.h"

#define DELL_LIM 1.5

static bool Isu_InFile(string file_name)
{
  ifstream rd(file_name.c_str());
  char buff[4096];
  string clef ;
  while ( (rd >> clef ) && (clef != "##WITH_u") )
    {
      // lire la fine de la ligne
      rd.getline(buff,4096);
    }
  if (clef == "##WITH_u")
    {
      
      cerr << "filter u in file " << file_name << endl ;
      return(true);
    }
  return(false);
}

static bool IsBV_InFile(string file_name)
{
  ifstream rd(file_name.c_str());
  char buff[4096];
  string clef ;
  while ( (rd >> clef ) && (clef != "##WITH_BV") )
    {
      // lire la fine de la ligne
      rd.getline(buff,4096);
    }
  if (clef == "##WITH_BV")
    {
      
      cerr << "filters B&V in file " << file_name << endl ;
      return(true);
    }
  return(false);
}


// (BV) + (u) + griz
void PegaseHeader(ofstream & prp, bool with_u, bool with_BV)
{

  string pegase_bands = "g eg r er i ei z ez" ;
  string pegase_filters = "effective_filter_g.dat,effective_filter_r.dat,effective_filter_i.dat,effective_filter_z.dat" ; 
  // sert au format et ajout des VEGA
  int nband = 4 ;
  if (with_u)
    {
      pegase_bands = "u eu "+pegase_bands; 
      pegase_filters = "effective_filter_u.dat,"+pegase_filters;
      nband++;
    }
  if (with_BV)
    {
      pegase_bands = "b eb v ev " +  pegase_bands ;
      pegase_filters = "effective_filter_B.dat,effective_filter_V.dat," + pegase_filters ;
      nband += 2 ;
    }
  PegaseHeader(prp, pegase_bands, pegase_filters, nband);
}



void PegaseHeader(ofstream & prp,string peg_bands, string peg_filters, int nband)
{
  prp << "##name z d ra dec " << peg_bands << endl ;
  prp << "##NBANDS " << nband << endl ;
  prp << "# FORMAT (a12,1x,f12.4,1x,f12.2,1x,f12.5,1x,f12.5,1x,"<< nband<<"(1x,f12.3,1x,f12.3))" << endl ;
  prp << "# FILTERS " << peg_filters << endl ;
  prp << "# CALIBTYPES " ;
  for(int k=0;k<nband;k++)
    {
      prp << "VEGA";
      if (k<nband-1) prp <<",";
    }
  prp << endl ;
  prp << setfill(' ');
  prp  << setiosflags(ios::fixed) ;
}
void PrintPegaseMag(ofstream & prp, double mag , double emag)
{
  prp << " " << setw(12) << mag << " "  << setw(12) << emag ;
}

void PrintPegase(ofstream & prp, string sn_name, double in_x, double in_y, double RA, double Dec)
{  
  prp << setw(12) << sn_name << " " ;
  prp<< setprecision(4);
  prp << setw(12) << in_x << " " ;
  prp<< setprecision(2);
  prp << setw(12)<< in_y << " " ;
  prp << setprecision(8) ;
  prp << setw(12) << RA << " " << setw(12) << Dec ;
  prp<< setprecision(3);
}


bool LookInFile(DictFile & in_file, string sn_name, DictFileIterator &vu)
{
  //cerr << "Looking for " << sn_name << endl ;
  for(DictFileIterator sn =  in_file.begin(); sn != in_file.end(); sn++)
    {
      DictFileEntry entry = *sn ;
      string name = entry.Value("name");
      if (name == sn_name)
	{
	  //cerr << sn_name << " found " << endl ;
	  vu = sn;
	  return(true);
	}
    }
  //cerr << sn_name << " NOT found " << endl ;
  return(false);
}


void ReadPegaseInFile(string file_name, DictFile & peg_file, bool & with_u, bool & with_BV)
{
  with_u = Isu_InFile(file_name);
  with_BV = IsBV_InFile(file_name);
  ifstream rd(file_name.c_str());
  char c ;
  char buff[4096];
  peg_file.AddKey("name");
  peg_file.AddKey("z");  //x
  peg_file.AddKey("delln"); //y
  peg_file.AddKey("ra");
  peg_file.AddKey("dec");
  if (with_BV)
    {
      peg_file.AddKey("mb");
      peg_file.AddKey("emb");  
      peg_file.AddKey("mv");
      peg_file.AddKey("emv");
    }
  if (with_u)
    {
      peg_file.AddKey("mu");
      peg_file.AddKey("emu");
    }
  peg_file.AddKey("mg");
  peg_file.AddKey("emg");
  peg_file.AddKey("mr");
  peg_file.AddKey("emr");
  peg_file.AddKey("mi");
  peg_file.AddKey("emi");
  peg_file.AddKey("mz");
  peg_file.AddKey("emz");

  
  while( rd >> c ) // to test eof
    {
      rd.unget() ; 
      rd.getline(buff,4096);
      if ( !(c == '#') )
	{
	   
	  peg_file.push_back(DictFileEntry(buff,peg_file));
	}
    }

}


void WritePegaseInFile(string file_name, DictFile & peg_file, bool with_u, bool with_BV)
{
  ofstream pr(file_name.c_str());
  PegaseHeader(pr,with_u, with_BV); 
  for(DictFileIterator sn =  peg_file.begin(); sn != peg_file.end(); sn++)
    {
      string name = sn->Value("name");
      double z = sn->Value("z");
      double delln = sn->Value("delln");
      double ra = sn->Value("ra");
      double dec = sn->Value("dec");
      PrintPegase(pr, name, z, delln, ra,dec);
      double mag, emag ;
      // B et V
      if (with_BV)
	{
	  PrintPegaseMag(pr, 99. , 99.);
	  PrintPegaseMag(pr, 99. , 99.);
	}
      if (with_u)
	{
	  mag = sn->Value("mu") ; emag = sn->Value("emu") ; 
	  PrintPegaseMag(pr, mag , emag);
	}
      mag = sn->Value("mg") ; emag = sn->Value("emg") ; 
      PrintPegaseMag(pr, mag , emag);
      mag = sn->Value("mr") ; emag = sn->Value("emr") ; 
      PrintPegaseMag(pr, mag , emag);
      mag = sn->Value("mi") ; emag = sn->Value("emi") ; 
      PrintPegaseMag(pr, mag , emag);
      mag = sn->Value("mz") ; emag = sn->Value("emz") ; 
      PrintPegaseMag(pr, mag , emag);
      pr << endl ;
    }
  pr.close();
}
  

  
void RunZpeg(string file_in, string file_out, string dir, bool is_age_cstr, 
	     bool is_z_fixed, bool with_u, bool with_BV)
{
  string template_file = dir + "/templates";
  string monpar_file = dir + "/monpar";
  if (is_age_cstr)
    {
      template_file += "_age_cstr" ;
      monpar_file += "_age_cstr" ;
    }
  if (is_z_fixed)
    {
      template_file += "_zfix" ;
      monpar_file += "_zfix" ;
    }
  if (with_u)
    template_file += "_u" ;
  if (with_BV)
    template_file += "_BV" ;
  template_file += ".tmp" ;
  monpar_file += ".par" ;

  if ( !FileExists( template_file ) )
    cerr << "missing " <<  template_file << " for running zpeg " << endl ;
  if ( !FileExists(monpar_file ) )
    cerr << "missing " <<  monpar_file << " for running zpeg " << endl ;


  string command = "cp " + dir + "/effective_filter*.dat . ; rm -f " + file_out + " ; zpeg " + file_in + " -t " + template_file + " -p " +  monpar_file + " -o " + file_out + " >&! zpeg.log" ;

  cerr << "Running : " << command ;
  system(command.c_str());

  string command_to_plot_1 = "cp " + template_file + " " + monpar_file + " . ";
  system(command_to_plot_1.c_str()); 
  string idl_file = "essai.pro" ;
  {ofstream pr(idl_file.c_str());
  pr << "plot_zpegfits_silence,'"+file_out+"',/zpeg_scen,/reread" << endl ;
  pr << "retall" << endl << "exit" << endl ;
  pr.close();
  }
  string command_to_plot_2 = "setenv IDL_STARTUP " + idl_file + " ; idl " ;
  cerr << command_to_plot_2 << endl ;
  //system(command_to_plot_2.c_str()); 
}


static void AddPegaseKeys(DictFile & peg_file_out,  string suffixe, bool with_u, bool with_BV)
{

  peg_file_out.AddKey("zp"+suffixe);
  peg_file_out.AddKey("zp_min"+suffixe);
  peg_file_out.AddKey("zp_max"+suffixe);
  peg_file_out.AddKey("sm"+suffixe);//stellar mass
  peg_file_out.AddKey("sm_min"+suffixe);
  peg_file_out.AddKey("sm_max"+suffixe);
  peg_file_out.AddKey("tmp"+suffixe);
  peg_file_out.AddKey("age"+suffixe);
  peg_file_out.AddKey("a_dist"+suffixe); //alpha parameter  
  peg_file_out.AddKey("ebv"+suffixe); 
  peg_file_out.AddKey("ssfr"+suffixe);
  peg_file_out.AddKey("ssfr_min"+suffixe);
  peg_file_out.AddKey("ssfr_max"+suffixe);
  peg_file_out.AddKey("sfr"+suffixe);
  peg_file_out.AddKey("sfr_min"+suffixe);
  peg_file_out.AddKey("sfr_max"+suffixe);
  peg_file_out.AddKey("scm"+suffixe);//stellar contain mass M star + M WD
  peg_file_out.AddKey("scm_min"+suffixe);
  peg_file_out.AddKey("scm_max"+suffixe);
  peg_file_out.AddKey("stm"+suffixe);// total mass formed stars
  peg_file_out.AddKey("stm_min"+suffixe);
  peg_file_out.AddKey("stm_max"+suffixe);
  if (with_u)
    {
      peg_file_out.AddKey("mabsu"+suffixe); // mag abs
    }
  peg_file_out.AddKey("mabsg"+suffixe);
  peg_file_out.AddKey("mabsr"+suffixe);
  peg_file_out.AddKey("mabsi"+suffixe);
  peg_file_out.AddKey("mabsz"+suffixe);
  if (with_BV)
    {
      peg_file_out.AddKey("tmp_mb"+suffixe); // template mag
      peg_file_out.AddKey("tmp_mb_min"+suffixe);
      peg_file_out.AddKey("tmp_mb_max"+suffixe);
      peg_file_out.AddKey("tmp_mv"+suffixe); // template mag
      peg_file_out.AddKey("tmp_mv_min"+suffixe);
      peg_file_out.AddKey("tmp_mv_max"+suffixe);
    }
  if (with_u)
    {
      peg_file_out.AddKey("tmp_mu"+suffixe); // template mag
      peg_file_out.AddKey("tmp_mu_min"+suffixe);
      peg_file_out.AddKey("tmp_mu_max"+suffixe);
    }
  peg_file_out.AddKey("tmp_mg"+suffixe);
  peg_file_out.AddKey("tmp_mg_min"+suffixe);
  peg_file_out.AddKey("tmp_mg_max"+suffixe);
  peg_file_out.AddKey("tmp_mr"+suffixe);
  peg_file_out.AddKey("tmp_mr_min"+suffixe);
  peg_file_out.AddKey("tmp_mr_max"+suffixe);
  peg_file_out.AddKey("tmp_mi"+suffixe);
  peg_file_out.AddKey("tmp_mi_min"+suffixe);
  peg_file_out.AddKey("tmp_mi_max"+suffixe);
  peg_file_out.AddKey("tmp_mz"+suffixe);
  peg_file_out.AddKey("tmp_mz_min"+suffixe);
  peg_file_out.AddKey("tmp_mz_max"+suffixe);
  return ;
}

static void is_ubv_in_outfile(string file_name, bool & with_u, bool & with_BV)
{
  ifstream rd(file_name.c_str());
  char c ;
  char buff[4096];
  with_u = false ;
  bool with_B = false,  with_V = false ;
  while( (rd >> c)  && ( c!= '#'))
    {
      rd.unget() ; 
      rd.getline(buff,4096);
      if (StringMatchPattern(buff,"effective_filter_u.dat"))
	{
	  with_u = true ;
	  cerr << "u in " << file_name << endl ;
	}
      if (StringMatchPattern(buff,"effective_filter_B.dat"))
	{
	with_B = true ;
	  cerr << "B in " << file_name << endl ;
	}
      if (StringMatchPattern(buff,"effective_filter_V.dat"))
	{
	with_V = true ;
	  cerr << "V in " << file_name << endl ;
	}
	
    }
  if ( with_B != with_V)
    cerr << "Trouble, B&V are't both absent or present" << endl ;
  with_BV = with_B ;
  return ;
}
void ReadPegaseOutFile(string file_name, DictFile & peg_file_out, bool & with_u, bool & with_BV)
{
  is_ubv_in_outfile(file_name, with_u,with_BV);
  ifstream rd(file_name.c_str());
  char c ;
  char buff[4096];
  peg_file_out.AddKey("name");
  peg_file_out.AddKey("z");  //x
  peg_file_out.AddKey("delln"); //y
  peg_file_out.AddKey("ra");
  peg_file_out.AddKey("dec");
  peg_file_out.AddKey("bestxi2");
  peg_file_out.AddKey("nband");
  peg_file_out.AddKey("nsol");
  peg_file_out.AddKey("zmax_vu");
  peg_file_out.AddKey("vmax");
  if (with_BV)
     {
       peg_file_out.AddKey("mb");
       peg_file_out.AddKey("emb");
       peg_file_out.AddKey("mv");
       peg_file_out.AddKey("emv");
     }
  // magnitudes en input
  if (with_u)
    {
      peg_file_out.AddKey("mu");
      peg_file_out.AddKey("emu");
    }
  peg_file_out.AddKey("mg");
  peg_file_out.AddKey("emg");
  peg_file_out.AddKey("mr");
  peg_file_out.AddKey("emr");
  peg_file_out.AddKey("mi");
  peg_file_out.AddKey("emi");
  peg_file_out.AddKey("mz");
  peg_file_out.AddKey("emz");
  //resultats best fit
  AddPegaseKeys(peg_file_out, "", with_u, with_BV);
  //resultats 4 autres best solutions
  AddPegaseKeys(peg_file_out, "1", with_u, with_BV);
  AddPegaseKeys(peg_file_out,"2", with_u, with_BV);
  AddPegaseKeys(peg_file_out, "3", with_u, with_BV);
  AddPegaseKeys(peg_file_out, "4", with_u, with_BV);

  
  while( rd >> c ) // to test eof
    {
      rd.unget() ; 
      rd.getline(buff,4096);
      if ( !(c == '#') )
	{
	  // on enleve les separateurs
	  string buffer = buff ;
	  RemovePattern(buffer,"|");
	  peg_file_out.push_back(DictFileEntry(buffer.c_str(),peg_file_out));
	}
    }

}


void Analysis_DZ(DictFile & peg_file_out, double & dz_mean, double &  dz_sig, double zmin, double zmax, bool quiet)
{
  double dz_max = 0.1 ;
  int Noff = 0 ;
  int n = 0 ;
  double *DZtab = new double[peg_file_out.size()] ;
  for(DictFileIterator sn =  peg_file_out.begin(); sn != peg_file_out.end(); sn++)
    {
      string name = sn->Value("name");
      double z = sn->Value("z");
      double delln = sn->Value("delln");
      if ((z > zmin) && (z < zmax) && (delln < DELL_LIM ))
	{
	  double zp  = sn->Value("zp");
	  double dz  = zp - z ;
	  if (fabs(dz) > dz_max)
	    {
	      if (!quiet) cerr << "OFFDZ : " << name << " " << z << " " << zp << " " << dz << endl ;
	      Noff++;
	    }
	  DZtab[n] = dz ;
	  n++;
	}
    }
  cerr <<  "OFFDZ : " << peg_file_out.size() <<  " " << Noff << endl ;
  double sigdz= 0 , mdz =0;
  double k = 3 ;
  int niter = 3 ;
  mdz = clipmean(DZtab, n , sigdz,k, niter);
  cerr << "dz moyenne, rms, nval : " << mdz << " " << sigdz << " " << n << endl ;
  delete[] DZtab ;
  dz_mean = mdz ;
  dz_sig = sigdz ;
  return ;
}
void Analysis_DZ_over_1plusZ(DictFile & peg_file_out, double & dz_mean, double &  dz_sig, double zmin, double zmax, bool quiet)
{
  double dzz_max = 0.2 ;
  int Noff = 0 ;
  int n = 0 ;
  double *DZtab = new double[peg_file_out.size()] ;
  string pegmag5[5] = {"mu","mg","mr","mi","mz"} ;
  int nok=0;
  for(DictFileIterator sn =  peg_file_out.begin(); sn != peg_file_out.end(); sn++)
    {
      string name = sn->Value("name");
      double z = sn->Value("z");
      double delln = sn->Value("delln");
      bool is_ok  = true ;
      if ( z < zmin )
	is_ok  = false ;
      if ( z > zmax )
	is_ok  = false ;
      if (delln > DELL_LIM )
	is_ok  = false ;
      for(int ib = 0 ; ib <5; ib++)
	{
	  if (sn->HasKey(pegmag5[ib]))
	    {
	      double em =  sn->Value("e"+pegmag5[ib]);
	      if (em > 0.1)
		is_ok  = false ;		  
	    }
	    }
      double mg =  sn->Value("mg");
      if (mg > 25.5 )
	is_ok  = false ;

      if (is_ok)
	{
	  nok++;
	  double zp  = sn->Value("zp");
	  double dz  = (zp - z)/(1+z) ;
	  if (fabs(dz) > dzz_max)
	    {
	      if (!quiet) cerr << "OFFDZ/(1+Z) : " << name << " " << z << " " << zp << " " << dz << endl ;
	      Noff++;
	    }
	  
	  DZtab[n] = dz ;
	  n++;
	    
	}
    }

  cerr <<  "DZ/(1+Z) : " << peg_file_out.size() <<  " nok : " << nok << "  noff : " << Noff << endl ;
  double sigdz= 0 , mdz =0;
  double k = 3 ;
  int niter = 3 ;
  mdz = clipmean(DZtab, n , sigdz,k, niter);
  cerr << "dz/(1+z) moyenne, rms, nval : " << mdz << " " << sigdz << " " << n << endl ;
  delete[] DZtab ;
  dz_mean = mdz ;
  dz_sig = sigdz ;
  return ;
}

void Analysis_DM(DictFile & peg_file_out, double *dm_mean, double *dm_sig, int nband, double zmin, double zmax, bool corr)
{
  string pegmag4[4] = {"mg","mr","mi","mz"} ;
  string pegmag5[5] = {"mu","mg","mr","mi","mz"} ;
  string *pegmag = pegmag4 ;
  if (nband==5) pegmag = pegmag5 ;

  for (int ib = 0 ; ib < nband; ib++)
    {
      
      int n = 0 ;
      double *DMtab = new double[peg_file_out.size()] ;
      for(DictFileIterator sn =  peg_file_out.begin(); sn != peg_file_out.end(); sn++)
	{
	  
	  double m_mes = sn->Value(pegmag[ib]);
	  double m_t = sn->Value("tmp_"+pegmag[ib]);
	  double z = sn->Value("z");
	  //cerr << m_mes << " " << m_t << endl ;
	  if ((z > zmin) && (z < zmax) && (m_t > 0 ) && (m_t < 90))
	    {
	      DMtab[n] = m_t - m_mes ;
	      n++;
	    }
	}
    
      double sigdm= 0 , mdm =0;
      double k = 3 ;
      int niter = 3 ;
      mdm = clipmean(DMtab, n , sigdm,k, niter);
      cerr << "dm " << pegmag[ib] << " moyenne, rms, nval : " << mdm << " " << sigdm << " " << n << endl ;
      delete[] DMtab ;
      dm_mean[ib] = mdm ;
      dm_sig[ib] = sigdm ;
    }

  if (corr)
    {
      for(DictFileIterator sn =  peg_file_out.begin(); sn != peg_file_out.end(); sn++)
	{
	  for (int ib = 0 ; ib < nband; ib++)
	    {	      
	      double m_mes = sn->Value(pegmag[ib]);
	      sn->ModKey(pegmag[ib], m_mes+dm_mean[ib]) ;
	    }
	}
    }
  return ;
}
