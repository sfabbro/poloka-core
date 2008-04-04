#include <iomanip> 

#include "sn_host_utils.h"
#include "pegase_utils.h"
#include "fileutils.h"
#include "reducedimage.h"
#include "fitsimage.h"

void SetSNHostTags(DictFile & lhost)
{
  // SN
  lhost.AddKey("name"); 
  lhost.AddKey("type");  
  lhost.AddKey("field");
  lhost.AddKey("season");
  lhost.AddKey("z");
  lhost.AddKey("ra_sn");
  lhost.AddKey("dec_sn");
  lhost.AddKey("xsn");
  lhost.AddKey("ysn");
  lhost.AddKey("ok");
  lhost.AddKey("nhost");
  lhost.AddKey("d_elln1");
  lhost.AddKey("d_elln2");
  lhost.AddKey("d_elln3");
  // Host (or neighbour)
  // in multimag cat
  lhost.AddKey("hostname"); 
  lhost.AddKey("host_n"); 
  lhost.AddKey("xh");
  lhost.AddKey("yh"); 
  lhost.AddKey("rah");
  lhost.AddKey("dech");
  lhost.AddKey("nh"); 
  lhost.AddKey("ah"); 
  lhost.AddKey("bh"); 
  lhost.AddKey("gxh");  
  lhost.AddKey("gyh"); 
  lhost.AddKey("gmxxh"); 
  lhost.AddKey("gmyyh"); 
  lhost.AddKey("gmxyh"); 
  lhost.AddKey("gmaah"); 
  lhost.AddKey("gmbbh");
  lhost.AddKey("gmxxlh"); 
  lhost.AddKey("gmyylh"); 
  lhost.AddKey("gmxylh"); 
  lhost.AddKey("gmaalh"); 
  lhost.AddKey("gmbblh");  
  lhost.AddKey("ell_a");
  lhost.AddKey("ell_b");
  lhost.AddKey("ell_angle");
  lhost.AddKey("circle_r");
  lhost.AddKey("g_ell_a");
  lhost.AddKey("g_ell_b");
  lhost.AddKey("g_ell_angle");
  lhost.AddKey("g_circle_r");
  lhost.AddKey("g_ok");
  lhost.AddKey("d");
  lhost.AddKey("d_ell");
  lhost.AddKey("d_elln");
  lhost.AddKey("mu");
  lhost.AddKey("emu");
  lhost.AddKey("mg");
  lhost.AddKey("emg");
  lhost.AddKey("mr");
  lhost.AddKey("emr");
  lhost.AddKey("mi");
  lhost.AddKey("emi");
  lhost.AddKey("mz");
  lhost.AddKey("emz");
  // info in cfhtls cat
  lhost.AddKey("d_cf");
  lhost.AddKey("ra_cf");
  lhost.AddKey("dec_cf");
  lhost.AddKey("z_cf");
  lhost.AddKey("zmin_cf");
  lhost.AddKey("zmax_cf");
  lhost.AddKey("tmp_cf");
  lhost.AddKey("ebv_cf");
  // zpeg result with zfix
  bool  with_u=true, with_BV=true ;
  cerr << "Adding pegase keys " << endl ;
  lhost.AddKey("nsolf"); // nsolutions trouvees zfix
  lhost.AddKey("nsol"); // nsolutions trouvees
  // zpeg result with zfix
  AddPegaseKeys(lhost,"f", with_u, with_BV);
  // 2eme solution 
  AddPegaseKeys(lhost,"1f", with_u, with_BV);
  // zpeg result with z not fix
  AddPegaseKeys(lhost,"", with_u, with_BV);
  // 2eme solution 
  AddPegaseKeys(lhost,"1", with_u, with_BV);
  lhost.AddKey("nsolna"); // nsolutions trouvees
  lhost.AddKey("nsolnaf"); // nsolutions trouvees
  // zpeg result with z fix no age cstr
  AddPegaseKeys(lhost,"naf", with_u, with_BV);
  // 2eme solution  
  AddPegaseKeys(lhost,"1naf", with_u, with_BV);
  // zpeg result with z not fix no age cstr
  AddPegaseKeys(lhost,"na", with_u, with_BV);
  // 2eme solution  
  AddPegaseKeys(lhost,"1na", with_u, with_BV);
  // zpeg result with z fix at z_cf ?
}


// a partir du dictfile, fabrications de toutes les images avec ds9
// attention, ttes les sn doivent correspondre a une meme image : meme champ, meme saison
// et elles doivent etre classees en no de hotes

void Vignettes(ReducedImage & rim, string dir, DictFile & lhost, int halfsize)
{
  for(DictFileIterator sn = lhost.begin(); sn != lhost.end(); )
    {
      string snname = sn->Value("name");
      int nhost =  sn->Value("nhost");
      double xsn =  sn->Value("xsn");
      double ysn =  sn->Value("ysn");
      string nom_img =  dir+"/"+snname+"_img.fits";
      string nom_img_seg =  dir+"/"+snname+"_seg.fits";
      double  d_numero_hote = sn->Value("nh"); 
      // c'est le numero de segmentation du 1st host
      // a inscrire dans l'entete fits de l'image de segmentation
      int numero_hote = (int) d_numero_hote ;
      if (FileExists(nom_img.c_str()))
	{
	  string com = "rm -f " + nom_img ;
	  system(com.c_str());
	}
      if (FileExists(nom_img_seg.c_str()))
	{
	  string com = "rm -f " + nom_img_seg ;
	  system(com.c_str());
	}
      // decoupage image ici pour imcopy
      int x1 = (int) (xsn+0.5) - halfsize + 1; 
      int y1 = (int) (ysn+0.5) - halfsize + 1;
      int x2 = x1 + 2 *  halfsize; 
      int y2 = y1 + 2 *  halfsize; 
      if ( x1 <= 0) x1 = 1;
      if ( y1 <= 0) y1 = 1;
      if ( x2 > rim.XSize()) x2 = rim.XSize();
      if ( y2 > rim.YSize()) y2 = rim.YSize();
      char ss[200];
      sprintf(ss,"[%d:%d,%d:%d]",x1,x2,y1,y2);
      string sss = ss ;
      string commande = "/afs/in2p3.fr/throng/snovae/softsnls/bin/imcopy '" + rim.FitsName() + sss +"' "+nom_img;  
      cerr << "Executing : " << commande << endl ;
      system(commande.c_str());
      string commande_seg = "/afs/in2p3.fr/throng/snovae/softsnls/bin/imcopy '" + rim.Dir()+ "/segmentation.fits" + sss +"' "+nom_img_seg;  
      cerr << "Executing : " << commande_seg << endl ;
      system(commande_seg.c_str());
      {
	FitsHeader head(nom_img_seg,RW);
	head.AddOrModKey("NSEG_HOST",numero_hote);
      }

      // decalage a apporter a x,y puisqu'a priori une image decoupee en (1,1) en bas a droite
      // est l'image de depart en (0,0)
      x1 -=1;
      y1 -=1;

      string nomreg =   dir+"/"+snname+"_img.reg";
      ofstream pr(nomreg.c_str());
      string color = " # color = red" ;
      double fact = 1. ;
      pr << "global color=black font=\"helvetica 12 bold\" width=2 select=1 edit=1 move=1 delete=1 include=1 fixed=0 source" << endl ;
      pr << "diamond point " << xsn -x1+1 << " " << ysn -y1+1 
     << " # color = red" << endl ;
      pr << "text " << x2-x1-20 << " " << y2-y1-20
	 << " {"<< snname << "} " << color << endl ;
      int n = 0 ;
      while ( n < nhost)
	{
	  n++;
	  string hostname = sn->Value("hostname");
	  cerr << "trace ellipses : " << snname << " " << hostname << " " << n << endl ;
	  double ell_a = sn->Value("ell_a");
	  double xh = sn->Value("xh");
	  double yh = sn->Value("yh");
	  double ell_b = sn->Value("ell_b");
	  double ell_angle = sn->Value("ell_angle");
	  double circle_r = sn->Value("circle_r");
	  if (circle_r < 0 ) // c'est une ellipse
	  pr << "ellipse " << xh - x1 + 1<< " " 
	     <<  yh  -y1 + 1<< " " << fact*ell_a << " " << fact*ell_b 
	     << " " <<  ell_angle << " # color = blue" << endl ;
	  else
	    pr << "circle " << xh -x1+ 1 << " " << yh  -y1 + 1<< " " << fact*circle_r <<  " # color = blue" << endl ;

	  // trace ellipse gaussienne
	  double g_ell_a = sn->Value("g_ell_a");
	  double gxh = sn->Value("gxh");
	  double gyh = sn->Value("gyh");
	  double g_ell_b = sn->Value("g_ell_b");
	  double g_ell_angle = sn->Value("g_ell_angle");
	  double g_circle_r = sn->Value("g_circle_r");
	  string colcol = " # color = green";
	  int isok = sn->Value("g_ok");
	  if (isok < 0 )
	    colcol = " # color = red";
	  if (g_circle_r < 0 ) // c'est une ellipse
	  pr << "ellipse " << gxh - x1 + 1<< " " 
	     <<  gyh  -y1 + 1<< " " << fact*g_ell_a << " " << fact*g_ell_b 
	     << " " <<  g_ell_angle << colcol << endl ;
	  else
	    pr << "circle " << gxh -x1+ 1 << " " << gyh  -y1 + 1<< " " << fact*g_circle_r <<  colcol << endl ;



	  // trace des numeros
	  pr << "text " << xh -x1 << " " << yh -y1
	     << " {"<< n << "} " << color << endl ;
	  
	  sn++;
	}
      pr.close() ;
      // fabrication du jpeg et du ps
      string nom_jpeg =  dir+"/"+snname+"_img.png";
      string nom_ps =  dir+"/"+snname+"_img.ps";
      string nom_jpeg3 =  dir+"/"+snname+"_img_zoom.png";
      string nom_ps3 =  dir+"/"+snname+"_img_zoom.ps";
      //ds9 -tile 03D1au_img.fits -histequ (ou -zscale)  -zoom 2 -cmap invert yes -region 03D1au_img.reg -saveas jpeg vv.jpeg -print filename "vv.ps" -print palette gray -print destination file -print -quit
      //string scale = " -zscale " ;
      //string scale = " -histequ " ;
      string scale = " -log -scalemode 99. -sqrt  " ;
      string commande_2 = "ds9 -view colorbar no -wcs wcs " + nom_img + " -zoom to fit -cmap invert yes " + scale + " -region " + nomreg + " -saveas png " + nom_jpeg + "  -print filename \"" + nom_ps +"\" -print destination file -print -quit" ;
      cerr << "Executing : " << commande_2 << endl ;
      system(commande_2.c_str());
      //scale = " -log -sqrt -mode 98. " ;
      scale = " -zscale " ;
      string commande_3 = "ds9 -view colorbar no -wcs wcs " + nom_img + " -zoom 6 -cmap invert yes " + scale + " -region " + nomreg + " -saveas png " + nom_jpeg3 + "  -print filename \"" + nom_ps3 +"\" -print destination file -print -quit" ;
      cerr << "Executing : " << commande_3 << endl ;
      system(commande_3.c_str());
      
    }
  return ;
}

void UpDateFromZPegOut(DictFile & lhost, DictFile & zpeg_out, string suffixe )
{
  if (lhost.size() != zpeg_out.size())
    {
      cerr << "Error : host list and zpeg out file must have same size : " 
	   << lhost.size() << " " << zpeg_out.size() << endl ;
      return ;
    }
  
  int Nkeys = 0 ; 
  bool with_u=true, with_BV=true ;
  string *tab_keys = PegaseKeys(Nkeys, with_u, with_BV); 
  DictFileIterator from_zpeg = zpeg_out.begin();
  DictFileIterator to_copy = lhost.begin();
  // copy 1st and 2nd sol
  string suff = suffixe;
  while(from_zpeg != zpeg_out.end() && to_copy != lhost.end())
    {
      for(int ii = 0 ; ii < Nkeys ; ii++) 
	    {
	      string val = (string) from_zpeg->Value(tab_keys[ii]+"") ;
	      to_copy->ModKey(tab_keys[ii]+suff,  val);
	      val = (string) from_zpeg->Value(tab_keys[ii]+"1") ;
	      to_copy->ModKey(tab_keys[ii]+"1"+suff, val);
	    }
	  string nsol = from_zpeg->Value("nsol");
	  to_copy->ModKey("nsol"+suff, nsol);
	  from_zpeg++;
	  to_copy++;
    }
  
  delete [] tab_keys ;
}

void EraseStringKeys(DictFile & file)
{
 for(DictFileIterator sn = file.begin(); sn != file.end(); sn++)
    {
      sn->ModKey("name",1) ;
      sn->ModKey("hostname", 1) ;
      sn->ModKey("field", 1) ;
      sn->ModKey("type", 1) ;
    }
 return ;
}


//--------- ecriture (pour paw par ex) de la version abregee : 
// sans les string, et avec juste la 1st solution

bool OKShort( string & key)
{
  if ( key == "hostname" )
    return false ;
  if ( key == "type" )
    return false ;
  if ( StringMatchPattern(key.c_str(),"*1") || StringMatchPattern(key.c_str(),"*1f") || StringMatchPattern(key.c_str(),"*1naf") || StringMatchPattern(key.c_str(),"*1na") )
    return(false); // enleve aussi d_elln1 !!
  return(true);
}
  
void WriteShort(DictFile & file,string nom)
{
  ofstream pr(nom.c_str());
  WriteShort(file, pr , false);
  pr.close();
  return ;
}

void WriteShort(DictFile & file, ofstream & pr , bool oui_append)
{
  map<int,string> tags;
  for (Dictionnary::const_iterator it =  file.Dict().begin(); it !=  file.Dict().end(); ++it)
    tags[it->second] = it->first;

  if (! oui_append ) // on ecrit le header
    {
      for (unsigned i = 0; i <file.Dict().size() ; ++i)
	{
	  if ( OKShort( tags[i] ) && (tags[i] != "name" ) )
	    pr << "# " << tags[i] << " : " << endl ;
	}
      pr << "#name : " << endl ;
      pr << "#end " << endl ;
    }

 for(DictFileIterator sn = file.begin(); sn != file.end(); sn++)
   {
     string field = sn->Value("field");
     if (field == "D1")
       sn->ModKey("field", 1) ;
     if (field == "D2")
       sn->ModKey("field", 2) ;
     if (field == "D3")
       sn->ModKey("field", 3) ;
     if (field == "D4")
       sn->ModKey("field", 4) ;
    }
 for(DictFileIterator sn = file.begin(); sn != file.end(); sn++)
   {
     for (unsigned i = 0; i <file.Dict().size() ; ++i)
       {
	 if ( OKShort(tags[i]) && (tags[i] != "name" ) )
	   {
	     string val =  sn->Value(tags[i]);
	     pr << val << " "  ;
	   }
       }
     pr <<  sn->Value("name") << endl ;
       
    }

 return ;
}


// on ne garde que les 1ers hotes
int KeepFirstHost(DictFile & file)
{
 for(DictFileIterator sn = file.begin(); sn != file.end(); )
    {
      double nn = sn->Value("host_n");
      if (nn >= 1)
	 sn = file.erase(sn);
      else
	sn++;
    }
 return(file.size());
}

// on ne garde que la sn de nom name
int KeepOneSN(DictFile & file, string name)
{
 for(DictFileIterator sn = file.begin(); sn != file.end(); )
    {
       string snname = sn->Value("name");
      if (snname != name )
	 sn = file.erase(sn);
      else
	sn++;
    }
 return(file.size());
}


// flagger a partir d'une liste
// file est la liste avec les first host seulement, file tout l aliste avec les hotes no i, au cas ou on doivent echanger.
void FlagSN(DictFile & file,DictFile &file_flags,DictFile &file_tout)
{
  for(DictFileIterator sn_flag = file_flags.begin(); sn_flag != file_flags.end(); sn_flag++)
    {
    string name_flag = sn_flag->Value("name");
    for(DictFileIterator sn = file.begin(); sn != file.end(); sn++)
      {
	string name = sn->Value("name");
	if (name == name_flag)
	  {
	    string val = sn_flag->Value("flag");
	    cerr << name << " flaggee " << val << endl ;
	    if (val == "NOT_OK")
	      sn->ModKey("ok",-100);
	    if (val == "NO_HOST")
	      sn->ModKey("ok",-1);
	    if (val == "SWAP_1_2")
	      {
		string name_h = name + "_1" ;
		// file et file_tout doivenet avoir le meme dictionnaire
		for(DictFileIterator sn2 = file_tout.begin(); sn2 != file_tout.end(); sn2++)
		  {
		    string name2 = sn2->Value("hostname");
		    if (name2 == name_h)
		      {
			string former = sn->Value("hostname") ;
			cerr << "Swaping " << name2 << " and " << former << endl ;
			map<int,string> tags;
			for (Dictionnary::const_iterator it =  file.Dict().begin(); it !=  file.Dict().end(); ++it)
			  tags[it->second] = it->first;
			for (unsigned i = 0; i <file.Dict().size() ; ++i)
			  if (tags[i] != "hostname")
			    {
			      string val = sn2->Value(tags[i]) ;
			      sn->ModKey(tags[i],val);
			    }
		      }
		  }
		//sn->ModKey("ok",-100);// pas encore de swap implemente
	      }
	  }
      }
    }
}


// trouver la SN dans un eliste avec le nom
bool FindSN(DictFile & file,string name_sn, DictFileIterator & sn_found)
{
  for(DictFileIterator sn = file.begin(); sn != file.end(); sn++)
    {
      string name = sn->Value("name");
      double nn = sn->Value("host_n"); 
      if ((name == name_sn) && (nn < 0.8))
	{
	  sn_found = sn ;
	  return true ;
	}
    }
  return false ;
}
   

// ============= pour analyse ============= 

// 1 si ok, 0 sinon
int OK_FLAG(DictFileIterator line)
{
  bool has_ok = line->HasKey("ok") ;
  if (has_ok)
    {
      double oktype=line->Value("ok") ;
      return( (int) oktype);
    }
  else
    return(1) ;
}

bool IsWithHost(DictFileIterator line)
{
  double zspec = line->Value("z") ;
  // double dell = line->Value("d_ell") ;
  double delln = line->Value("d_elln") ;
  //if ((zspec > 0) && (dell < 5 ));
 if ((zspec > 0) && (delln < DELLN_LIM ))
    return true ;
  else
    return false ;
}

// pour les nearby, pour que ce soit coherent
int Host_Ell_Spir_Type(int type)
{
  int id ;
  if ( type > -1  )
    {
      if (IsPegaseEll(type)) id = 0 ;
      else id = 1 ;
    }
  else
    id=-100 ; // pb
  return(id);
}

// pour les nearby, pour que ce soit coherent
int Host_Earlier_Later_Type(int type, int type_cut)
{
  int id ;
  if ( type > -1  )
    {
      if ( type <= type_cut ) id = 0 ;
      else id = 1 ;
    }
  else
    id=-100 ; // pb
  return(id);
}

int Host_Ell_Spir_Type(DictFileIterator line)
{
  int id ;
  if ( IsWithHost(line) && (OK_FLAG(line) == 1) )
    {
      if (IsPegaseEll(line)) id = 0 ;
      else id = 1 ;
    }
  if ( !IsWithHost(line) && (OK_FLAG(line) == 1) )
    {
      id = -1 ;
    }
  if ( OK_FLAG(line) != 1 )
    id = OK_FLAG(line);
  return(id);
}

int Host_red_blue_Type(DictFileIterator line, double BV_cut)
{
  int id ;
  if ( IsWithHost(line) && (OK_FLAG(line) == 1) )
    {
     double mabf = line->Value("mabf");
     double mavf  = line->Value("mavf");
     if ( mabf - mavf >= BV_cut ) // B grand dc flux B petit dc rouge
	    id = 0 ;
      else id = 1 ;
    }
  if ( !IsWithHost(line) && (OK_FLAG(line) == 1) )
    {
      id = -1 ;
    }
  if ( OK_FLAG(line) != 1 )
    id = OK_FLAG(line);
  return(id);
}


int Host_passive_active_Type(DictFileIterator line, double ssfr_cut)
{
  int id ;
  if ( IsWithHost(line) && (OK_FLAG(line) == 1) )
    {
      double sfrf = line->Value("ssfrf") ;
     if ( sfrf < ssfr_cut) //passive
	    id = 0 ;
      else id = 1 ;
    }
  if ( !IsWithHost(line) && (OK_FLAG(line) == 1) )
    {
      id = -1 ;
    }
  if ( OK_FLAG(line) != 1 )
    id = OK_FLAG(line);
  return(id);
}


int PegaseId_3(DictFileIterator line)
{
  int id = -1 ;
  if (IsWithHost(line))
    {
      if (IsPegaseEll(line))
	id = 0 ;
      if ( IsPegaseEarlySpir(line))
	id = 1 ;
      if ( IsPegaseLateSpir(line))
	id = 2 ;
      if ( IsPegaseIrr(line))
	id = 3 ;
    }
  return id ;
}
  
int PegaseId_7(DictFileIterator line)
{
  int id = -1 ;
  if (IsWithHost(line))
    {
      if (IsPegaseEll(line)) //type 1 2 -> 0
	id = 0 ;
      if ( IsPegaseSpir(line)) // type 3 4 5 6 7 -> 1 2 3 4 5  
	{
	  double tmpf = line->Value("tmpf") ;
	  id = (int) (tmpf+0.5)-2 ;
	}
      if ( IsPegaseIrr(line)) // type 8 9  -> 6
	id =6  ;
    }
  return id ;
}
string PegaseId_3String(int id)
{
  switch(id){
  case -1 :
    return("not_found");
    break ;
  case 0 :
    return("Ell");
    break ;
  case 1 :
    return("EarlySpir");
    break ;
  case 2 :
    return("LateSpir");
    break ;
  case 3 :
    return("Irr/SB");
    break ;
  }
  return("");
}
 
string PegaseId_7String(int id)
{
  switch(id){
  case -1 :
    return("not_found");
    break ;
  case 0 :
    return("E/S0");
    break ;
  case 1 :
    return("Sa");
    break ;
  case 2 :
    return("Sb");
    break ;
  case 3 :
    return("Sbc");
    break ;
  case 4 :
    return("Sc");
    break ;
  case 5 :
    return("Sd");
    break ;
  case 6 :
    return("Irr-SB");
    break ;
  }
  return("");
}

int PegaseType(DictFileIterator line)
{ 
  int id ;
  if ( IsWithHost(line) && (OK_FLAG(line) == 1) )
    {
      double tmpf = line->Value("tmpf") ;
      id = (int) tmpf ;
    }
  if ( !IsWithHost(line) && (OK_FLAG(line) == 1) )
    {
      id = -1 ;
    }
  if ( !(OK_FLAG(line) == 1) )
    id = OK_FLAG(line);
  return(id);
}
string PegaseString(DictFileIterator line)
{ 
  int id = PegaseType(line);
  return(PegaseString(id));
}

string PegaseString(int id)
{ 
 
  switch(id){
  case -1 :
    return("no id. host");
    break ;
  case -100 :
    return("prob. case");
    break ;
  case 1 :
    return("E");
    break ;
  case 2 :
    return("S0");
    break ;
  case 3 :
    return("Sa");
    break ;
  case 4 :
    return("Sb");
    break ;
  case 5 :
    return("Sbc");
    break ;
  case 6 :
    return("Sc");
    break ;
  case 7 :
    return("Sd");
    break ;
  case 8 :
    return("Irr-SB");
    break ;
  case 9 :
    return("Irr-SB");
    break ;
  }
  return("");
}


// ============= pour CFHTLS


bool Host_is_in_CFHTLS(DictFileIterator line)
{
  bool d_cf = line->Value("d_cf");
  if ((d_cf > 0 ) && (d_cf < DLIM_CFHTLS))
    return true ;
  else
    return false ;
}

//E-S0
bool IsCFHTLSEll(DictFileIterator line)
{  
  double tmp_cf = line->Value("tmp_cf") ;
  if ( tmp_cf <=16 )
    return true ;
  else
    return false ;
}

//Sbc
bool IsCFHTLSSbc(DictFileIterator line)
{
  double tmp_cf = line->Value("tmp_cf") ;
  if (( tmp_cf >=17 ) && ( tmp_cf <=29 ))
    return true ;
  else
    return false ;
}

//Scd
bool IsCFHTLSScd(DictFileIterator line)
{
  double tmp_cf = line->Value("tmp_cf") ;
  if (( tmp_cf >=30 ) && ( tmp_cf <=43 ))
    return true ;
  else
    return false ;
}

//Irr + SB 
bool IsCFHTLSIrrSB(DictFileIterator line)
{
  double tmp_cf = line->Value("tmp_cf") ;
  if (( tmp_cf >=44 ) && ( tmp_cf <=62 ))
    return true ;
  else
    return false ;
}





int CFHTLS_Id(DictFileIterator line)
{
 int id = -1 ;
  if (Host_is_in_CFHTLS(line))
    {
      if (IsCFHTLSEll(line))
	id = 0 ;  
      if (IsCFHTLSSbc(line))
	id = 1 ;
      if (IsCFHTLSScd(line))
	id = 2 ;
      if (IsCFHTLSIrrSB(line))
	id = 3 ;
    }
  return(id);
}
string  CFHTLS_String(DictFileIterator line)
{
  int id = CFHTLS_Id(line);
 if (id == -1)
    return("not_found");
  if (id == 0)
    return("Ell");
  if (id == 1)
    return("Sbc");
  if (id == 2)
    return("Scd");
  if (id == 3)
    return("Irr-SB");
    
  return("");

}



int CFHTLSHost_Ell_Spir_Type(DictFileIterator line)
{
  if(IsWithHost(line) && Host_is_in_CFHTLS(line))
    {if (IsCFHTLSEll(line))  return(0); else return(1);}
  else
    return(-1);
}
	

// ================== les warnings


// 3 hotes potentiels tout pres
bool Warning_Uncertain_Attribution(DictFileIterator line)
{
  double d_elln1 = line->Value("d_elln");
  double d_elln2 = line->Value("d_elln2");
  double d_elln3 = line->Value("d_elln3");
  string name =  line->Value("name");
  string field = line->Value("field") ;
  double season = line->Value("season") ;
  bool warning = false ;
  if ( (d_elln1 > 0 ) && (d_elln1 < DELLN_LIM ) && ( d_elln2 >0 ) && (d_elln2 <  DELLN_LIM ))
    {
      warning = true ;
      if ((d_elln3 > 0 ) &&  (d_elln3 <  DELLN_LIM ))
	cerr << name << " " << field << " S" << int(season) <<  " Uncertain (3 close hosts): " 
	     << d_elln1 << " " << d_elln2  << " " << d_elln3 << endl ; 
      else	
	cerr << name << " " << field << " S" << int(season) <<  " Uncertain (2 close hosts): " 
	     << d_elln1 << " " << d_elln2 << endl ; 
    }
  return(warning);
}

bool Warning_Disagreement_w_CFHTLS(DictFileIterator line)
{
  string name = line->Value("name") ;
  string field = line->Value("field") ;
  string season = line->Value("season") ;
  int host_type = Host_Ell_Spir_Type(line) ; // 0 ou 1, -1 si pas trouve
  int cfhtls_type = CFHTLSHost_Ell_Spir_Type(line) ; // 0 ou 1, -1 si pas
  if (host_type > 0 ) // de ttes facons, cfhtls host presque tjrs trouve car c'est la sn !
    if (host_type != cfhtls_type)
      {
	cerr << name << " " << field << " " << season << " DESACCORD " 
	     << PegaseString(line)  <<  " " << CFHTLS_String(line)  << endl ; 
	return(true);
      }
  return(false);
}
