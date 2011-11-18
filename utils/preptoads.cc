#include <iostream>
#include <vector>

#include "imageutils.h"
#include "fitsimage.h"
#include "histo1d.h"

static void usage(const char *progName)
{
  cerr << progName << " [options] <FITS files>\n" 
       << "   options:  -g  multiply the image by the gain \n"
       << "          :  -s  set the sky and rms\n"
       << "          :  -f  set saturation level\n"
       << "          :  -b  write and scale image with BITPIX=16 \n"
       << "          :  -c  clean all TOADS FITS keys\n"
       << "          :  -l  prepare LBL calibrated images\n";
}

static void SetGain(FitsImage &Im)
{
  if (Im.HasKey("TOADGAIN"))
    {
      double gain = Im.KeyVal("TOADGAIN");
      if (gain == 1) 
	{
	  cout << "   is already gain multiplied \n";
	  return;
	}
      cout << " Multiplying by gain = " << gain << endl;
      Im *= gain;
      Im.AddOrModKey("OLDGAIN",gain,"original gain before TOADS (counts/ADU)");  
      Im.AddOrModKey("TOADGAIN",1,"Hopefully image is in photoelectrons");
    }
  else cerr << "   does not have TOADGAIN \n";  
}

static void SetLbl(FitsHeader &Im)
{
  if (Im.HasKey("TOADGAIN"))
    {
      double gain = Im.KeyVal("TOADGAIN");
      if (gain == 1) 
	{
	  cout << "   is already gain multiplied \n";
	  return;
	}
      Im.AddOrModKey("OLDGAIN",gain,"original gain before TOADS (counts/ADU)");
      Im.AddOrModKey("TOADGAIN",1,"Hopefully image is in photoelectrons");
    }
  else cerr << "   does not have TOADGAIN \n";  

  if (Im.HasKey("SKY") && Im.HasKey("SIGMA"))
    {
      double sky = Im.KeyVal("SKY");
      double sigma = Im.KeyVal("SIGMA");
      Im.AddOrModKey("SKYLEV",sky," Copy of SKY key (LBL style)");
      Im.AddOrModKey("SKYSIGEX",sigma," Copy of SIGMA key (LBL style)");
    }

  if (Im.HasKey("WELLDEPT"))
    {
      double satur = Im.KeyVal("WELLDEPT");
      Im.AddOrModKey("SATURLEV",satur," Copy of WELLDEPT key (LBL style)");
    }
      
  Im.AddCommentLine("Image prepared from LBL style stuff");
}

void SetSky(FitsImage &Im)
{
  Pixel sky;
  if (!Im.HasKey("SKYLEV") && !Im.HasKey("SKYSIGEX"))
    {
      Pixel sig;
      Im.SkyLevel(&sky, &sig);
      Im.AddOrModKey("SKYLEV",sky," Original sky level");
      Im.AddOrModKey("SKYSIGEX",sig," Original sky sigma");
      float gain = Im.KeyVal("TOADGAIN");
      float rdnoise = Im.KeyVal("TOADRDON");
      float sigth = sqrt(sky*gain+rdnoise*rdnoise)/gain;
      Im.AddOrModKey("SKYSIGTH",sigth," Theoretical sigma (sqrt(sky*gain+rdnoise^2)/gain");
      cout << " Setting sky = "<< sky << " sigma = " << sig << " th sigma = " << sigth << endl;
    }
  else sky = Im.KeyVal("SKYLEV");
  return;
}

double SetOverscan(FitsHeader &Head)
{
  double overscan = 100.;
  if (Head.HasKey("OVERSCAN")) 
    {
      string soverscan = Head.KeyVal("OVERSCAN");
      size_t pos = soverscan.find("mean"); // HC, iraf type reductions
      if (pos == string::npos) overscan = Head.KeyVal("OVERSCAN");
      else 
	{
	  soverscan.erase(0,pos+5);
	  overscan = atof(soverscan.c_str());
	  Head.AddOrModKey("OVERSCAN", overscan,"mean subtracted overscan value");
	}      
    }
  cout << " Setting overscan = " << overscan << endl;
  return overscan;
}


static void SetSaturation(FitsImage &Im)
{
  double saturation = ComputeSaturation(Im);
  // WARNING WELLDEPTH CAN CHANGE
  double welldepth = 65536.;
  if (Im.HasKey("WELLDEPT")) welldepth = Im.KeyVal("WELLDEPT");
  double overscan = SetOverscan(Im);
  double factor = 1;
  double gain = Im.KeyVal("TOADGAIN");
  if ((gain == 1) && Im.HasKey("OLDGAIN")) factor = Im.KeyVal("OLDGAIN"); 
  double osatur = (welldepth-overscan)*factor;
  if (fabs(osatur-saturation) > 0.1*osatur) 
    {
      cout << " Found satur= " << saturation << " Forcing saturation from welldepth.\n";
      saturation = osatur;
    }
  
  cout << " Setting saturation level = "<< saturation <<endl;
  Im.AddOrModKey("SATURLEV",saturation," Current saturation level");
}


void SetBitpix16(FitsImage &Im)
{
  Pixel minIm=0,maxIm;
  if (!Im.HasKey("SATURLEV")) SetSaturation(Im);
  maxIm = Im.KeyVal("SATURLEV");
  // cout << " min = " << minIm<< " max = " << maxIm << endl,
  Im.EnforceMinMax(minIm,maxIm);
  Im.AddOrModKey("BITPIX",16," forced by TOADS");
}

void CleanToadsKeys(FitsHeader &Head)
{
  const int nkeys = 39;
  static string keyList[nkeys] = 
    {"TOADGAIN", "TOADEXPO", "TOADINST", "TOADFILT", 
     "TOADRDON", "TOADEQUI", "TOADDECL", "TOADRASC",
     "TOADPIXS", "TOADBAND", "TOADAIRM", "TOADDATE",
     "TOADUTIM", "TOADTYPE", "TOADOBJE", "TOADCHIP",
     "TOADNAMP", "TOADPZPT", "TOADJULI", 
     "SEXSKY", "SEXSIGMA", "SESEEING", "BACK_SUB",
     "SATURLEV", "SKYLEV", "SKYSIGEX", "SKYSIGTH",
     "OLDGAIN", "ZEROUSNO", "USNOSB23", "SCALFACT", 
     "KERNCHI2", "BACKLEV", "ORIGSATU", "FLATNOIS", 
     "SIGPSF", "SESIGX", "SESIGY", "SERHO"
    };
  for (int i=0; i<nkeys; ++i)
    {
      string keyname = keyList[i];
      if (Head.HasKey(keyname))
	{
	  string keyval = Head.KeyVal(keyname);
	  cout << " Removing " << keyname << "=" << keyval << endl;
	  Head.RmKey(keyname);
	}
    }      
}


int main(int nargs, char **args)
{
  if (nargs <=1) {usage(args[0]); exit(-1);}
  vector<string> toDo;
  bool set_gain=false, set_sky=false, set_satur=false, 
       set_b16=false,set_clean=false, set_lbl=false;
  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] != '-') 
	{
	  toDo.push_back(arg);
	  continue;
	}
      switch (arg[1])
	{
	case 'g': set_gain = true; break;
	case 's': set_sky = true; break;
	case 'f': set_satur = true; break;
	case 'b': set_b16 = true; break;
	case 'c': set_clean = true; break;
	case 'l': set_lbl = true; break;
	default : usage(args[0]); exit(1);
	}
    }

  int nim = toDo.size();
  cout << " " << nim << " images to process \n";
  for (int i=0; i<nim; ++i)
    {
      string name = toDo[i];
      cout << " Processing " << name << endl;
      if (set_lbl) 
	{
	  FitsHeader imFits(name,RW);
	  SetLbl(imFits);
	  continue;
	}

      FitsImage imFits(name,RW);
      if (set_clean) CleanToadsKeys(imFits);
      if (set_gain) SetGain(imFits); 
      if (set_satur) SetSaturation(imFits);
      if (set_sky) SetSky(imFits);
      if (set_b16) 
	{
	  SetBitpix16(imFits);
	  imFits.Write();
	}
    }
  return EXIT_SUCCESS;
}

