#include <cmath>

#include "sestar.h"
#include "standardstar.h"
#include "fitsimage.h"
#include "fileutils.h"
#include "dbimage.h"
#include "usnoutils.h"
#include "wcsutils.h"

#ifndef M_PI
#define  M_PI           3.14159265358979323846  /* pi */
#endif


// To match an image with the Standard catalogue

/***************************************************************************/
/***************************************************************************/
/***********************************      **********************************/
/*********************************** MAIN **********************************/
/***********************************      **********************************/
/***************************************************************************/
/***************************************************************************/

void ImageProcess(const string &fitsFileName, const string &sexCatalogName, double *absorption, double *zeropoints, double *err_zeropoints, int &nzeros, StandardStarList * &standardList)
{

   FitsHeader header(fitsFileName);

   cout << endl << "Matching with Standard catalog..."<< endl;
   cout << fitsFileName << endl;
   cout << "INSTRUMENT = " << header.KeyVal("TOADINST") << endl;
   cout << "CCD no = " << header.KeyVal("TOADCHIP") << endl;
   cout << "FILTER = " << header.KeyVal("TOADFILT") << endl;
   cout << "OBJECT = " << header.KeyVal("TOADOBJE") << endl;
   cout << endl;
   
   if (!HasLinWCS(header))
     {
       cerr << "The image " << fitsFileName << " doesn't have WCS !! " << endl;
       return;
     }
   standardList = GetSelectedStandardStarList(header);
   if (!standardList || standardList->size() == 0)
     {
       cout << " No standard star selected. Stop here !" << endl;
       return;
     }
   cout << " Read " << standardList->size() << " objects from standard table in the selected window" << endl;

   SEStarList sestarList(sexCatalogName);
   cout << " Read " << sestarList.size() << " objects from SexCat " << endl;   

   //   StandardColor couleur = GetColor(header);
   //   double expo = header.KeyVal("TOADEXPO");

   string filt = header.KeyVal("TAODBAND");
   double zerotheo = ReadTheoZeroPoint(header);
   //   double zerotheo = ReadTheoZeroPoint(header) - 2.5*log10(header.KeyVal("TOADEXPO"));
   //   zerotheo += 2.5*log(header.KeyVal("OLDGAIN"));
   cout <<" Theorical Zero Point : "<< zerotheo <<endl;
   //   double zerotheo = header.KeyVal("ZEROUSNO");
   //   zerotheo -= 2.5*log(expo);
   //   cout << "zero theo : "<<zerotheo1<<" zero usno : "<<zerotheo2<<endl;

   GetStandardZeroPoint(standardList, sestarList, header, zeropoints, err_zeropoints, nzeros);

   for (int i=0 ; i<nzeros ; i++) absorption[i]=zerotheo - zeropoints[i];
}


static void usage(const char *pgname)
{
  cout << pgname << ' ' << " < standard dbimage name .... > (for matching with Standard list)" << endl;
  cout <<"      -f : create a list of absorption values + infos from used stars "<<endl;
  cout <<"      -ff : fit the airmass and color term and save them in a file "<<endl;
}


int main(int argc, char **argv)
{

  int fitparam = 0;
  int outfile = 0;
  int chip_ref = 0;

  list<string> dbImageList;
  if (argc < 2)
    {
      usage(argv[0]);  return 0;
    }
  
  for (int i=1;  i<argc; ++i)
    {
      char *arg = argv[i];
      if (arg[0] != '-') 
	{
	  dbImageList.push_back(string(arg));
	  continue;
	}
      if (arg[1] == 'f') 
	{
	  outfile=1;
	  if (arg[2] == 'f') fitparam=1;
	}
      if (arg[1] == 'h') {usage(argv[0]); exit(0);}
      if (arg[1] == 'c')
	{
	  i++;
	  chip_ref = atoi(argv[i]);
	  break;
	}
    }
  
  //  cout << " the all list of images will be added to calculate ONE zero point. Be sure you add the right images..."<<endl;
  
  double *totalabs = new double[dbImageList.size()*10];
  double *normtotalabs = new double[dbImageList.size()*10];
  // double total_weight = 0.0;
  int ntotal=0;
  string band;
  string chip;
  string instName;
  string date;
  string outname="";
  double zerotheo;

  StandardStarList totalstdstarlist;
  for (list<string>::iterator i=dbImageList.begin(); i!= dbImageList.end(); ++i)
    {
      string name = *i;
      DbImage dbimage(name);
      if (!dbimage.IsValid())
	{
	  cerr << " Be careful ! " << name << " must be an image name ! " << endl;
	}
      
      string sexCatName = dbimage.ImageCatalogName(SExtractor);
      if (!FileExists(sexCatName.c_str()))
	{
	  cerr << "The SExtractor Catalog associated to image " << name << " doesn't exist !! " << endl;
	  continue;
	}
      
      string fitsFileName = dbimage.FitsImageName(Calibrated);
      if (!FileExists(fitsFileName))
	{
	  cerr << "The  calibrated fits image " << name << " doesn't exist !! " << endl;
	  continue;
	}
      FitsHeader head(fitsFileName);
      band = head.KeyVal("TOADBAND");
      chip = head.KeyVal("TOADCHIP");
      instName = head.KeyVal("TOADINST");
      date = head.KeyVal("TOADDATE");

      zerotheo = ReadTheoZeroPoint(head);
      //      zerotheo = ReadTheoZeroPoint(head) - 2.5*log10(head.KeyVal("TOADEXPO"));
      //      zerotheo += 2.5*log(head.KeyVal("OLDGAIN"));

      double abs[100];
      double zeros[100];
      double err_zeros[100];
      int nzeros=0;
      StandardStarList * stdstarlist;
      ImageProcess(fitsFileName, sexCatName, abs, zeros, err_zeros, nzeros, stdstarlist);
      for (int i=0; i<nzeros; i++)
	{
	  totalabs[ntotal]=abs[i];
	  //	 normtotalzeros[ntotal]=zeros[i]/err_zeros[i];
	  normtotalabs[ntotal]=abs[i];
	  //	 total_err_zeros[ntotal]=err_zeros[i];
	  //	 total_weight += 1.0/err_zeros[i];
	  ntotal++;
	}
      if (nzeros > 0)
	{
	  for (StandardStarIterator si = stdstarlist->begin(); si != stdstarlist->end(); si++)
	    {
	      StandardStar * pstar = (StandardStar *) *si;
	      StandardStar * newstdstar = new StandardStar(*pstar);
	      totalstdstarlist.push_back(newstdstar);
	    }
       }
    }
  
  // zeropoint is the clipped-mean at 2sigmas with 1 iterations.
  double absorption,zeropoint,errzero;

  if (outfile)
	{
	  unsigned int l = date.length();
	  for (unsigned int i=0; i<l; i++)
	    if (date[i]=='/') date[i]='_';
	  outname += "Abs_"+instName+"_"+band+"_"+date+".list";
	  ofstream out(outname.c_str());
	  out << "#absorb : absorption" << endl;
	  out << "#fluxpsec : flux per second" << endl;
	  out << "#x : x " << endl;
	  out << "#y : y " << endl;
	  out << "#flag : flag " << endl;
	  out << "#end" << endl;
	  int i = 0;
	  for (StandardStarIterator si = totalstdstarlist.begin(); si != totalstdstarlist.end(); si++)
	    {
	      StandardStar * pstar = (StandardStar *) *si;
	      out << normtotalabs[i] <<" " <<pstar->fluxpersec <<" "<< pstar->x <<" "<< pstar->y <<" "<<pstar->Flag() << endl;
	      i++;
	    }
    }
  absorption = clipmean(normtotalabs,ntotal,errzero,2.0,2);
  // zeropoint = (ntotal/total_weight)*clipmean(normtotalzeros,ntotal,3,errzero,2);
  // errzero *= ntotal/total_weight;
  
  // errzero /=sqrt(ntotal);

  cout << "******************************************************"<<endl;
  cout << "******************************************************"<<endl;
  cout << "  ABS("<<band<<") FOUND USING "<< ntotal<<" STARS : "<< absorption << " +- " << errzero << endl;
  zeropoint = zerotheo - absorption;
  cout << "  ZP("<<band<<","<<chip_ref<<") : "<< zeropoint << " +- " << errzero << endl;
  cout << "******************************************************"<<endl;
  cout << "******************************************************"<<endl;
  delete [] totalabs;
 
  if (fitparam)
    {
      if (ntotal >= 3)
	{
	  double a0, a1, a2, err_a0, err_a1, err_a2;
	  if (outfile)
	    {
	      unsigned int l = date.length();
	      for (unsigned int i=0; i<l; i++)
		if (date[i]=='/') date[i]='_';
	      outname += "ZP_"+instName+"_"+chip+"_"+band+"_"+date+".txt";
	    }
	  a0 = FitZeroPoint(totalstdstarlist,a1,a2,zeropoint,errzero,err_a0,err_a1,err_a2,outname);
	  cout << endl;
	  cout << "Result of fit (with the following convention) : "<<endl;
	  cout << "======================================================"<<endl;
	  cout << "M = -2.5log10(flux/time) + zp + a1*Color + a2*Airmass)"<<endl;
	  cout << "======================================================"<<endl;
	  cout << "zp("<< band <<","<< chip <<") = "<< a0 << " +/- "<< err_a0 
	       << "  a1 = " << a1 << " +/- " << err_a1
	       << "  a2 = "<<a2<< " +/- " << err_a2
	       << endl;
	}
      else
	{
	  cout << "not enough data to fit the 3 parameters (only "<< ntotal <<" stars)."<< endl;
	}
    }
  return 1;
}

