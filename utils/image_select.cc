#include <sstream>
#include <iostream>
#include <vector>

#include "stringlist.h"
#include "frame.h"
#include "gtransfo.h"
#include "fileutils.h"
#include "astroutils.h"
#include "wcsutils.h"
#include "fitsimage.h"
#include "dbimage.h"
#include "imageutils.h"

static void usage(const char *progName)
{

  cerr << "Usage: " << progName << " [OPTIONS] <image(s)>\n"
       << "  list FITS images with criteria of keys in header \n"
       << "   OPTIONS are \n"
       << "     -db : flag for <image(s)> to be dbimages (calibrated, then raw) \n"
       << "     -band BAND : select by specified one letter filter BAND \n"
       << "     -chip CHIP : select by specified chip CHIP \n"
       << "     -ra HH:MM:SS  : select by right ascencion \n"
       << "     -dec DD:MM:SS : select by declination \n"
       << "     -rad VALUE : select an radius in arcmin around specified ra and dec (default is 1)\n"
       << "     -date DD/MM/YYYY-DD/MM/YYYY : select between 2 UTC dates \n"
       << "     -area VALUE : minimum overlap area in armin2 to use with -ref (default is 20) \n"
       << "     -ref IMAGE : match WCS among images and a reference IMAGE \n"
       << "     -s 'KEY OP VAL' : select with a formatted expression: \n"
       << "                        KEY -> FITS key to select \n"
       << "                        OP  -> operator (==,<,>,<=,>=,~=) ~= is search string\n"
       << "                        VAL -> value to search for\n"
       << "   Exemple:\n"
       << "     image_select -chip 00 -date '22/03/2003-24/03/2003' *.fits\n"
       << "           or \n"
       << "     image_select -db -band i -s 'TOADOBJE ~= D3' -s 'SESEEING <= 1'  `dbls 2003-03-03`\n\n";

  exit(-1);
}

bool DecodeFitsExpression(const FitsHeader &Head, const string&Exp)
{
  istringstream iss(Exp);
  string keyname, op, keyval;
  iss >> keyname >> op >> keyval;
  if (!Head.HasKey(keyname)) return false;

  if (op == "==") 
    return string(Head.KeyVal(keyname)) == keyval;
  if (op == "~=")
    return string(Head.KeyVal(keyname)).find(keyval) != string::npos;
  double dval = atof(keyval.c_str());
  if (op == "<") 
    return double(Head.KeyVal(keyname)) <  dval;
  if (op == "<=") 
    return double(Head.KeyVal(keyname)) <= dval;
  if (op == ">") 
    return double(Head.KeyVal(keyname)) >  dval;
  if (op == ">=") 
    return double(Head.KeyVal(keyname)) >= dval;
  return true;
}

int main(int nargs, char **args)
{
  // if nothing is given
  if (nargs < 2){usage(args[0]);}

  // default stuff
  bool withDbImages = false;
  double radius=1, ra=-1, dec=-1, minarea = 20;
  string datestr="", refstr="";

  StringList toSelect, explist;

  // loop over arguments
  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      // images
      if (arg[0] != '-') 
	{
	  toSelect.push_back(arg);
	  continue;
	}

      // options
      arg++;
      sscanf(arg, "%s", arg);
      if (strcmp(arg,"db")==0)  { withDbImages = true; continue;}
      if (strcmp(arg,"band")==0){ ++i; explist.push_back(string("TOADBAND == ")+args[i]); continue;}
      if (strcmp(arg,"chip")==0){ ++i; explist.push_back(string("TOADCHIP == ")+args[i]); continue;}
      if (strcmp(arg,"s")==0) { ++i; explist.push_back(args[i]); continue;}
      if (strcmp(arg,"date")==0){ ++i; datestr = args[i]; continue;}
      if (strcmp(arg,"ref")==0){ ++i; refstr = args[i]; continue;}
      if (strcmp(arg,"area")==0) { ++i; minarea = atof(args[i]); continue;}
      if (strcmp(arg,"rad")==0) { ++i; radius = atof(args[i])/60./* convert to deg */; continue;}
      if (strcmp(arg,"ra")==0)  { ++i; ra = RaStringToDeg(args[i]); continue;}
      if (strcmp(arg,"dec")==0) { ++i; dec = DecStringToDeg(args[i]); continue;}

      // unrecognized option
      usage(args[0]);      
    }

  double julmin,julmax,ramin,ramax,decmin,decmax;

  bool setdate = (datestr.length() > 0);
  if (setdate)
    {	
      int d1,m1,y1,d2,m2,y2;
      if (sscanf(datestr.c_str(),"%d/%d/%d-%d/%d/%d", &d1,&m1,&y1,&d2,&m2,&y2) !=6)
	{
	  cerr << " Error parsing date format. Please use DD/MM/YYYY-DD/MM/YYYY \n";
	  exit(1);
	}
      julmin = JulianDay(d1,m1,y1);
      julmax = JulianDay(d2,m2,y2);
    }

  bool setref = (refstr.length() > 0);
  FitsHeader *refhdr = 0;
  if (setref)
    {
      if (withDbImages)
	{
	  refstr = DbImage(refstr).FitsImageName(Calibrated);
	  if (!FileExists(refstr)) refstr = DbImage(refstr).FitsImageName(Raw);
	}
      refhdr = new FitsHeader(refstr);      
    }

  bool setrad = (ra != -1 && dec !=-1);
  Frame sourceFrame;
  if (setrad)
    {
      double cosdec = cos(M_PI*dec/180);
      ramin = ra-radius*cosdec;
      ramax = ra+radius*cosdec;
      decmin = dec-radius;
      decmax = dec+radius;
      sourceFrame = Frame(ramin,decmin, ramax, decmax);
    }

  for (StringCIterator it=toSelect.begin(); it!=toSelect.end(); ++it)
    {
      string name(*it);
      if (withDbImages)
	{
	  if (!DbImage(*it).IsValid()) continue;
	  name = DbImage(*it).FitsImageName(Calibrated);
	  if (!FileExists(name)) name = DbImage(*it).FitsImageName(Raw);
	}
      if (!FileExists(name)) continue;

      FitsHeader head(name);
      if (!head.IsValid()) continue;

      // decode dates
      bool date_is_true = true;
      if (setdate) 
	{
	  double juldate = JulianDay(head);
	  date_is_true = (juldate >= julmin) && (juldate <= julmax);
	}
      
      // decode ref matching
      bool ref_is_true = true;
      if (setref)
	{	  
	  double overlap = Arcmin2Overlap(*refhdr, head);
	  // cout << overlap << endl;
	  ref_is_true = (overlap >= minarea);
	}

      // decode radius
      bool rad_is_true = true;
      if (setrad)
	{
	  Frame pixFrame(head);
	  Gtransfo *Pix2RaDec;
	  rad_is_true = false;
	  if (WCSFromHeader(head, Pix2RaDec))
	    {
	      Frame sidFrame = ApplyTransfo(pixFrame,*Pix2RaDec);
	      Frame overlap = sidFrame*sourceFrame;
	      cout << "sidFrame " << sidFrame << endl;
	      cout << "sourceFrame " << sourceFrame << endl;
	      rad_is_true = (overlap.Area()>0);
	    }
	  else
	    cerr << " cannot select on Ra,Dec : no WCS in " << head.FileName() << endl;
	}
      
      // decode expressions
      bool exp_is_true = true;
      for (StringCIterator ie=explist.begin(); ie!=explist.end(); ++ie)
	exp_is_true = DecodeFitsExpression(head, *ie) && exp_is_true;

      // print image if selected
      if ((exp_is_true) && (date_is_true) && (rad_is_true) && (ref_is_true)) cout << *it << endl;
    }

  if (refhdr) delete refhdr;

  return EXIT_SUCCESS;
}
