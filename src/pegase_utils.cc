#include <iomanip> 

#include "pegase_utils.h"
#include "fileutils.h"
#include "vutils.h"
#include <cmath>

#define DELL_LIM 1.5
#define ZPEG_COMMANDE "/afs/in2p3.fr/throng/snovae/softsnls/pegase/zpeg_5.1/src/zpeg "
//#define ZPEG_COMMANDE "/afs/in2p3.fr/throng/snovae/softsnls/pegase/zpeg/project/bin/zpeg "

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

static bool IsBVUR_InFile(string file_name)
{
  ifstream rd(file_name.c_str());
  char buff[4096];
  string clef ;
  while ( (rd >> clef ) && (clef != "##WITH_BVUR") )
    {
      // lire la fine de la ligne
      rd.getline(buff,4096);
    }
  if (clef == "##WITH_BVUR")
    {
      
      cerr << "filters B&V&U&R in file " << file_name << endl ;
      return(true);
    }
  return(false);
}


// (BVUR) + (u) + griz
void PegaseHeader(ofstream & prp, bool with_u, bool with_BVUR, bool use_AB)
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
  if (with_BVUR)
    {
      pegase_bands = "b eb v ev uj uje rc rce " +  pegase_bands ;
      pegase_filters = "effective_filter_B.dat,effective_filter_V.dat,effective_filter_U.dat,effective_filter_R.dat," + pegase_filters ;
      nband += 4 ;
    }
  if (with_u)
    prp << "##WITH_u" << endl ;
  if (with_BVUR)
    prp << "##WITH_BVUR" << endl ;
  PegaseHeader_(prp, pegase_bands, pegase_filters, nband, use_AB);
}



void PegaseHeader_(ofstream & prp,string peg_bands, string peg_filters, int nband, bool use_AB)
{
  prp << "##name z d ra dec " << peg_bands << endl ;
  prp << "# FORMAT (a12,1x,f12.4,1x,f12.2,1x,f12.5,1x,f12.5,1x,"<< nband<<"(1x,f12.3,1x,f12.3))" << endl ;
  prp << "# FILTERS " << peg_filters << endl ;
  prp << "# CALIBTYPES " ;
  if (! use_AB )
    {
      for(int k=0;k<nband;k++)
	{
	  prp << "VEGA";
	  if (k<nband-1) prp <<",";
	}
    }
  else
    {
      for(int k=0;k<nband;k++)
	{
	  prp << "AB";
	  if (k<nband-1) prp <<",";
	}
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





void ReadPegaseInFile(string file_name, DictFile & peg_file, bool & with_u, bool & with_BVUR)
{
  with_u = Isu_InFile(file_name);
  with_BVUR = IsBVUR_InFile(file_name);
  ifstream rd(file_name.c_str());
  char c ;
  char buff[4096];
  peg_file.AddKey("name");
  peg_file.AddKey("z");  //x
  peg_file.AddKey("delln"); //y
  peg_file.AddKey("ra");
  peg_file.AddKey("dec");
  if (with_BVUR)
    {
      peg_file.AddKey("mb");
      peg_file.AddKey("emb");  
      peg_file.AddKey("mv");
      peg_file.AddKey("emv");
      peg_file.AddKey("muj");
      peg_file.AddKey("emuj");  
      peg_file.AddKey("mrc");
      peg_file.AddKey("emrc");
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


void WritePegaseInFile(string file_name, DictFile & peg_file, bool with_u, bool with_BVUR, bool is_host_list, bool use_AB)
{
  ofstream pr(file_name.c_str());
  PegaseHeader(pr,with_u, with_BVUR, use_AB); 
  for(DictFileIterator sn =  peg_file.begin(); sn != peg_file.end(); sn++)
    {
      string name ;
      double z, delln, ra, dec ;
      if (is_host_list)
	{
	  if (sn->HasKey("hostname"))
	    {
	       string s = sn->Value("hostname");
	       name = s ;
	    }
	  else
	    {
	      if (sn->HasKey("name"))
		{
		  string s = sn->Value("name");
		  name = s ;
		}
	      else
		cerr << "no name in dictfile : pegase_utils " << endl ;
	    }
	  z = sn->Value("z");
	  delln = sn->Value("d_elln");
	  ra = sn->Value("rah");
	  dec = sn->Value("dech");
	}
      else
	{
	  string s  = sn->Value("name");
	  name = s ;
	  z = sn->Value("z");
	  delln = sn->Value("delln");
	  if ( !(sn->HasKey("ra")) && !(sn->HasKey("RA")) )
	    cerr << "list has no ra or RA key " << endl ;
	  if ( !(sn->HasKey("dec")) && !(sn->HasKey("Dec")) )
	    cerr << "list has no dec or Dec key " << endl ;
	  if (sn->HasKey("ra"))
	      ra = sn->Value("ra");
	  if (sn->HasKey("RA"))
	      ra = sn->Value("RA");
	  if (sn->HasKey("dec"))
	      dec = sn->Value("dec");
	  if (sn->HasKey("Dec"))
	      dec = sn->Value("Dec");
	}
      PrintPegase(pr, name, z, delln, ra,dec);
      double mag, emag ;
      // B et V
      if (with_BVUR)
	{
	  PrintPegaseMag(pr, 99. , 99.);
	  PrintPegaseMag(pr, 99. , 99.);
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
	     bool is_z_fixed, bool with_u, bool with_BVUR)
{ 
  string file_in1 = file_in ;
  string file_out1 = file_out ;

  string workdir = "." ;
  char *tdir = getenv("ZPEG_WORK_DIR");
  if (tdir) 
    {
      workdir = tdir ;
      file_in1 = workdir + "/temp_pegase.in" ;
      file_out1 = workdir + "/temp_pegase.out" ; 
      string command1 = "cp " + file_in + " " + file_in1;
      cerr << "Running : " << command1 << endl;
      system(command1.c_str());
    }


  string template_file = "templates";
  string monpar_file = "monpar"; 
  string idl_file =  workdir+"/essai" ;

 
  if (is_age_cstr)
    {
      template_file += "_age_cstr" ;
      monpar_file += "_age_cstr" ;
    }
  if (is_z_fixed)
    {
      template_file += "_zfix" ;
      monpar_file += "_zfix" ; 
      idl_file += "_zfix" ;
    }
  if (with_u)
    template_file += "_u" ;
  if (with_BVUR)
    template_file += "_BVUR" ;

  template_file += ".tmp" ;
  monpar_file += ".par" ; 
  idl_file += ".pro" ;
  bool is_template_file_here = true ;
  if ( !FileExists( dir+"/"+template_file ) )
    {
      cerr << "missing " <<  dir+"/"+template_file << " for running zpeg " << endl ;
      cerr << workdir+"/"+template_file << " will be created"  << endl ;
      is_template_file_here = false ;
    }
  if ( !FileExists(dir+"/"+monpar_file ) )
    cerr << "missing " <<  dir+"/"+monpar_file << " for running zpeg " << endl ;


  string command = "cp " + dir + "/effective_filter*.dat  " + workdir;
  cerr << "Running : " << command << endl;
  system(command.c_str());

  command = "cd " + workdir + " ; rm -f " + file_out1 ;
  cerr << "Running : " << command << endl;
  system(command.c_str());

  if (is_template_file_here)
    {
      command = "cd " + workdir + " ; ln -fs " + dir + "/" + template_file + " " + template_file ;
      cerr << "Running : " << command<< endl ;
      system(command.c_str());
    }
  else
    {
      // comme le template n'est pas ds le repertoire des file de pegase
      // on le creera sur place - il faut s'assurer qu'il n'y en a pas 
      // deja un qui pourtrait ne pas correspondre
      command = "cd " + workdir + " ; rm -f " + template_file ;
      cerr << "Running : " << command<< endl ;
      system(command.c_str());
    }

  command = "cd " + workdir + " ;ln -fs " + dir + "/" + monpar_file + " " + monpar_file ;
  cerr << "Running : " << command << endl ;
  system(command.c_str());

  
  command = "cd " + workdir + " ; " + ZPEG_COMMANDE + file_in1 + " -t " + template_file + " -p " +  monpar_file + " -o " + file_out1 + " > zpeg.log 2>&1" ;
  cerr << "Running : " << command << endl;
  system(command.c_str());

  command = "cd " + workdir + " ; cat zpeg.log " ;
  cerr << "Running : " << command << endl;
  system(command.c_str());

  {ofstream pr(idl_file.c_str());
  pr << "plot_zpegfits_silence,'"+file_out1+"',/zpeg_scen,/reread" << endl ;
  pr << "retall" << endl << "exit" << endl ;
  pr.close();
  }
  string command_to_plot_2 = "setenv IDL_STARTUP " + idl_file + " ; idl " ;
  cerr << "to obtain the plots with idl : " << command_to_plot_2 << endl ;
  //system(command_to_plot_2.c_str()); 

  command = "cd " + workdir + " ; rm -f  " + monpar_file + " " + monpar_file ;
  cerr << "Running : " << command << endl ;
  system(command.c_str());
 command = "cd " + workdir + " ; rm -f  effective_filter*.dat" ;
  cerr << "Running : " << command << endl ;
  system(command.c_str());

  if (tdir) 
    {     
      string command1 = "cp " + file_out1 + " " + file_out;
      cerr << "Running : " << command1 << endl;
      system(command1.c_str());
    }



}


string *PegaseKeys(int & N, bool with_u, bool with_BVUR, bool first_sol)
{
  N = 42 ;
  if (with_BVUR) 
    N += 16 ;
  if (with_u) 
    N += 4 ;

  string * tab_keys = new string[500] ;

  int n = 0 ;

  /*if (first_sol)
    {
      
      tab_keys[n] = "xi2" ; n++ ;
      tab_keys[n] = "nsol" ; n++ ;
      N +=2 ;
      }*/

  tab_keys[n] = "zp" ; n++ ;
 tab_keys[n] = "zp_min" ; n++ ;
 tab_keys[n] = "zp_max" ; n++ ;
 tab_keys[n] = "sm" ; n++ ;//stellar mass
 tab_keys[n] = "sm_min" ; n++ ;
 tab_keys[n] = "sm_max" ; n++ ;
 tab_keys[n] = "tmp" ; n++ ;
 tab_keys[n] = "age" ; n++ ;
 tab_keys[n] = "a_d" ; n++ ; //alpha parameter  
 tab_keys[n] = "ebv" ; n++ ;
 
 tab_keys[n] = "xi2" ; n++ ; 

 tab_keys[n] = "ssfr" ; n++ ; // sfr - stellar mass = sfr - sm
 tab_keys[n] = "ssfr_min" ; n++ ;
 tab_keys[n] = "ssfr_max" ; n++ ;
 tab_keys[n] = "sfr" ; n++ ;
 tab_keys[n] = "sfr_min" ; n++ ;
 tab_keys[n] = "sfr_max" ; n++ ;
 tab_keys[n] = "scm" ; n++ ;//stellar contain mass M star + M WD
 tab_keys[n] = "scm_min" ; n++ ;
 tab_keys[n] = "scm_max" ; n++ ;
 tab_keys[n] = "stm" ; n++ ;// total mass formed stars
 tab_keys[n] = "stm_min" ; n++ ;
 tab_keys[n] = "stm_max" ; n++ ;
 tab_keys[n] = "sta" ; n++ ;// ages of  formed stars
 tab_keys[n] = "sta_min" ; n++ ;
 tab_keys[n] = "sta_max" ; n++ ;
  if (with_BVUR)
    {
     tab_keys[n] = "mab" ; n++ ; // mag abs B
     tab_keys[n] = "mav" ; n++ ; // mag abs V
     tab_keys[n] = "mauj" ; n++ ; // mag abs U
     tab_keys[n] = "marc" ; n++ ; // mag abs R
    }
  if (with_u)
    {
     tab_keys[n] = "mau" ; n++ ; // mag abs
    }
 tab_keys[n] = "mag" ; n++ ;
 tab_keys[n] = "mar" ; n++ ;
 tab_keys[n] = "mai" ; n++ ;
 tab_keys[n] = "maz" ; n++ ;
 // mag abs min
  if (with_BVUR)
    {
     tab_keys[n] = "mab_min" ; n++ ;
     tab_keys[n] = "mav_min" ; n++ ;
     tab_keys[n] = "mauj_min" ; n++ ;
     tab_keys[n] = "marc_min" ; n++ ;
    }
  if (with_u)
    {
     tab_keys[n] = "mau_min" ; n++ ;
    }
 tab_keys[n] = "mag_min" ; n++ ;
 tab_keys[n] = "mar_min" ; n++ ;
 tab_keys[n] = "mai_min" ; n++ ;
 tab_keys[n] = "maz_min" ; n++ ;
 // mag abs max
  if (with_BVUR)
    {
     tab_keys[n] = "mab_max" ; n++ ; 
     tab_keys[n] = "mav_max" ; n++ ; 
     tab_keys[n] = "mauj_max" ; n++ ; 
     tab_keys[n] = "marc_max" ; n++ ; 
    }
  if (with_u)
    {
     tab_keys[n] = "mau_max" ; n++ ;
    }
 tab_keys[n] = "mag_max" ; n++ ;
 tab_keys[n] = "mar_max" ; n++ ;
 tab_keys[n] = "mai_max" ; n++ ;
 tab_keys[n] = "maz_max" ; n++ ;



  if (with_BVUR)
    {
     tab_keys[n] = "t_mb" ; n++ ; // template mag
     //tab_keys[n] = "t_mb_min" ; n++ ;
     //tab_keys[n] = "t_mb_max" ; n++ ;
     tab_keys[n] = "t_mv" ; n++ ; // template mag
     //tab_keys[n] = "t_mv_min" ; n++ ;
     //tab_keys[n] = "t_mv_max" ; n++ ;
     tab_keys[n] = "t_muj" ; n++ ; // template mag
     tab_keys[n] = "t_mrc" ; n++ ; // template mag
    }
  if (with_u)
    {
     tab_keys[n] = "t_mu" ; n++ ; // template mag
     //tab_keys[n] = "t_mu_min" ; n++ ;
     //tab_keys[n] = "t_mu_max" ; n++ ;
    }
 tab_keys[n] = "t_mg" ; n++ ;
 //tab_keys[n] = "t_mg_min" ; n++ ;
 //tab_keys[n] = "t_mg_max" ; n++ ;
 tab_keys[n] = "t_mr" ; n++ ;
 //tab_keys[n] = "t_mr_min" ; n++ ;
 //tab_keys[n] = "t_mr_max" ; n++ ;
 tab_keys[n] = "t_mi" ; n++ ;
 //tab_keys[n] = "t_mi_min" ; n++ ;
 //tab_keys[n] = "t_mi_max" ; n++ ;
 tab_keys[n] = "t_mz" ; n++ ;
 //tab_keys[n] = "t_mz_min" ; n++ ;
 //tab_keys[n] = "t_mz_max" ; n++ ;
 cerr << "Pegase keys : " << n << " " << N << endl ;
 //for (int ii = 0 ; ii < n ; ii++)
 //  cerr << tab_keys[ii] << endl ;
 if (n != N )
   cerr << "Problem in pegase keys : " << n << " " << N << endl ;
 
 return (tab_keys);
}

void AddPegaseKeys(DictFile & peg_file_out,  string suffixe, bool with_u, bool with_BVUR)
{
  
  int Nkeys = 0 ;
  cerr << "Creating pegase keys " << endl ;
  string *tab_keys  = PegaseKeys(Nkeys, with_u, with_BVUR);
  cerr << Nkeys << " keys created " << endl ;
  for(int ii = 0 ; ii < Nkeys ; ii++)
    {
      string key = tab_keys[ii]+suffixe ;
      //cerr << key << endl ;
      peg_file_out.AddKey(key);
    }
  delete [] tab_keys ;
  return ;
}





void CheckPegaseKeys(DictFile & peg_file_out,  string suffixe, bool with_u, bool with_BVUR)
{
  
  int Nkeys = 0 ;
  cerr << "Check pegase keys " << endl ;
  string *tab_keys  = PegaseKeys(Nkeys, with_u, with_BVUR);
  cerr << Nkeys << " keys created " << endl ;
  for(int ii = 0 ; ii < Nkeys ; ii++)
    {
      string key = tab_keys[ii]+suffixe ;

      for(DictFileIterator it =  peg_file_out.begin() ; it != peg_file_out.end() ; it++)
	{
	 
	  string res = it->Value(key); 
	  //cerr << res << endl ;
	  if ( res == "*******")
	    {
	      cerr << "### " << key << " " << res << endl ;
	      double u = -99 ;
	      it->ModKey(key,u);
	    }
	  if (StringMatchPattern(res.c_str(),"*+*")   && !StringMatchPattern(res.c_str(),"*E+*") && !StringMatchPattern(res.c_str(),"+*") )
	    {
	      cerr << "### " << key << " " << res << endl ;
	      double u = -99 ;
	      it->ModKey(key,u);
	    }
	}
    }
  delete [] tab_keys ;
  return ;
}




static void is_ubv_in_outfile(string file_name, bool & with_u, bool & with_BVUR)
{
  ifstream rd(file_name.c_str());
  char c ;
  char buff[4096];
  with_u = false ;
  bool with_B = false,  with_V = false, with_U = false, with_R = false ;
  while( (rd >> c)  && ( c== '#'))
    {
      rd.unget() ; 
      rd.getline(buff,4096);
      string sbuff = buff ;
      //cerr << buff << endl ;
      if (sbuff.find("effective_filter_u.dat") < sbuff.length())
	{
	  with_u = true ;
	  cerr << sbuff << endl ;
	  cerr << "u in " << file_name << endl ;
	}
      if (sbuff.find("effective_filter_B.dat") < sbuff.length())
	{
	  with_B = true ;
	  cerr << sbuff << endl ;
	  cerr << "B in " << file_name << endl ;
	}
      if (sbuff.find("effective_filter_V.dat") < sbuff.length())
	{
	  with_V = true ;
	  cerr << sbuff << endl ;
	  cerr << "V in " << file_name << endl ;
	}
      if (sbuff.find("effective_filter_U.dat") < sbuff.length())
	{
	  with_U = true ;
	  cerr << sbuff << endl ;
	  cerr << "U in " << file_name << endl ;
	}
      if (sbuff.find("effective_filter_R.dat") < sbuff.length())
	{
	  with_R = true ;
	  cerr << sbuff << endl ;
	  cerr << "R in " << file_name << endl ;
	}
	
    }
  if ( ! (  (with_B==false && with_V==false && with_U==false && with_R==false)
	    || (with_B==true && with_V==true && with_U==true && with_R==true) ) )
    cerr << "Trouble, B&V&U&R aren't all absent or present" << endl ;
  with_BVUR = with_B ;
  return ;
}
void ReadPegaseOutFile(string file_name, DictFile & peg_file_out, bool & with_u, bool & with_BVUR, bool oui_append)
{
  is_ubv_in_outfile(file_name, with_u,with_BVUR);
  ifstream rd(file_name.c_str());
  char c ;
  char buff[4096];
  if ( ! oui_append )
    {
  peg_file_out.AddKey("name");
  peg_file_out.AddKey("z");  //x
  peg_file_out.AddKey("delln"); //y
  peg_file_out.AddKey("ra");
  peg_file_out.AddKey("dec");
  peg_file_out.AddKey("bxi2"); // best xi2
  peg_file_out.AddKey("nband");
  peg_file_out.AddKey("nsol");
  peg_file_out.AddKey("zmax_vu");
  peg_file_out.AddKey("vmax");
  if (with_BVUR)
     {
       peg_file_out.AddKey("mb");
       peg_file_out.AddKey("emb");
       peg_file_out.AddKey("mv");
       peg_file_out.AddKey("emv");
       peg_file_out.AddKey("muj");
       peg_file_out.AddKey("emuj");
       peg_file_out.AddKey("mrc");
       peg_file_out.AddKey("emrc");
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
  AddPegaseKeys(peg_file_out, "", with_u, with_BVUR);
  //resultats 4 autres best solutions
  AddPegaseKeys(peg_file_out, "1", with_u, with_BVUR);
  AddPegaseKeys(peg_file_out,"2", with_u, with_BVUR);
  AddPegaseKeys(peg_file_out, "3", with_u, with_BVUR);
  AddPegaseKeys(peg_file_out, "4", with_u, with_BVUR);
    }
  
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
  CheckPegaseKeys(peg_file_out, "", with_u, with_BVUR);
  CheckPegaseKeys(peg_file_out, "1", with_u, with_BVUR);
  CheckPegaseKeys(peg_file_out, "2", with_u, with_BVUR);
  CheckPegaseKeys(peg_file_out, "3", with_u, with_BVUR);
  CheckPegaseKeys(peg_file_out, "4", with_u, with_BVUR);
}



//============== decodages des types ===============


//E-S0 =0
bool IsPegaseEll(DictFileIterator line)
{  
  double tmpf = line->Value("tmpf") ;
  return(IsPegaseEll(tmpf));
}

bool IsPegaseEll(double tmpf)
{  
  
  if ( (tmpf == 1 ) || (tmpf == 2 ))
    return true ;
  else
    return false ;
}

bool IsPegaseSa(DictFileIterator line)
{  
  double tmpf = line->Value("tmpf") ;
  if ( fabs(tmpf-3)<1.e-2 ) //==3
    return true ;
  else
    return false ;
}
bool IsPegaseSb(DictFileIterator line)
{  
  double tmpf = line->Value("tmpf") ;
  if ( fabs(tmpf-4)<1.e-2 ) //==4
    return true ;
  else
    return false ;
}
bool IsPegaseSbc(DictFileIterator line)
{  
  double tmpf = line->Value("tmpf") ;
  if ( fabs(tmpf-5)<1.e-2 ) //==5
    return true ;
  else
    return false ;
}

bool IsPegaseSc(DictFileIterator line)
{  
  double tmpf = line->Value("tmpf") ;
  if ( fabs(tmpf-6)<1.e-2 ) //==6
    return true ;
  else
    return false ;
}


bool IsPegaseSd(DictFileIterator line)
{  
  double tmpf = line->Value("tmpf") ;
  if (  fabs(tmpf-7)<1.e-2 ) //==7
    return true ;
  else
    return false ;
}

bool IsPegaseSpir(DictFileIterator line)
{  
  double tmpf = line->Value("tmpf") ;
  if ( (tmpf >= 3 ) && (tmpf <= 7 ))
    return true ;
  else
    return false ;
}


//Early Spir : Sa Sb Sbc = 1
bool IsPegaseEarlySpir(DictFileIterator line)
{  
  double tmpf = line->Value("tmpf") ;
  if ( (tmpf >= 3 ) && (tmpf <= 5 ))
    return true ;
  else
    return false ;
}

//Late Spir : Sc Sd = 2
bool IsPegaseLateSpir(DictFileIterator line)
{  
  double tmpf = line->Value("tmpf") ;
  if ( (tmpf >= 6 ) && (tmpf <= 7 ))
    return true ;
  else
    return false ;
}

//Irr SB = 3
bool IsPegaseIrr(DictFileIterator line)
{  
  double tmpf = line->Value("tmpf") ;
  if ( (tmpf >= 8 ) && (tmpf <= 9 ))
    return true ;
  else
    return false ;
}



