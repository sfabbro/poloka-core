#include <iostream>
#include <string>
#include <cstdio>

#include "fitsimage.h"
#include "fitstoad.h"

#include "fileutils.h"
void split_name(const string &FileName, string &Dir, string &Base, string &Type)
{
Dir = DirName(FileName);
Base = BaseName(FileName);
char name[256];
sprintf(name,"%s",Base.c_str());
char *gz = strstr(name,".gz");
if (gz) *gz = '\0';
char *dot = strrchr(name,'.');
if (dot)
  {
  Type = dot+1;
  *dot='\0';
  Base = name;
  }
else
  {
  Type = "";
  }
}


static void usage()
{
  cout << " usage : split_fits file [file ...] [output_directory] (input files can ge gzipped)" << endl;
}

int main(int nargs, char **args)
{
string outDirName = "./";
int nIn;
 if (nargs <=1) {usage(); exit(-1);}
if (nargs >=3 && IsDirectory(args[nargs-1]))
  {
  outDirName = AddSlash(args[nargs-1]);
  nIn = nargs-1;
  }
else
  {
  nIn = nargs; 
  }
cout << " writing output files in " << outDirName << endl;  


for (int i=1; i<nIn; ++i)
  {
    string fileName = args[i];
    
    FitsHeader inFile(fileName);
    if (!inFile.IsValid())
      {
	cerr << " file " << fileName << " is not a valid fits file " << endl;
        usage();
        continue;
      }

    if (!inFile.HasKey("EXTEND") || bool(inFile.KeyVal("EXTEND")) == false )
      {
	cerr << " The file " << fileName << " has no fits extensions " << endl;
	exit (1);
      }
    //    int next = inFile.KeyVal("NEXTEND");
    //    cout << " the file " << fileName << " has " << next << " extensions" << endl;


    // THE FOLLOWING LINE IS NECESSARY FOR SPLITING INT-WFC IMAGES :
    // THE EXTENTION HEADER DON'T HAVE ANY KEYWORD ALLOWING TO KNOW YOU
    // ARE DEALING WITH THE WFC. AND THE FOLLOWING LINE AFFECT INSTRUMENT
    // KEYWORD.
    cout << " splitting " << fileName << " from " << TelInstName(inFile) << endl;

    string dir,base,ext;
    split_name(fileName, dir, base, ext);
    //  store the main header to copy it to extensions.
    FitsHeader mainHeader(inFile, outDirName+base+"_mh."+ext);    

    int nhdu = inFile.NHDU();
    for (int i= 2; i<= nhdu; ++i)
      {
	if (inFile.MoveHDU() == 0) break;
	//	cout << inFile ;
	int chip = inFile.KeyVal("TOADCHIP",true);
	char schip[256];
        sprintf(schip,"%d",chip);
	//	string outFileName = outDirName + base + "_" + string(schip) + "." + ext;
	char outFileName[512];
        const char *file_name_format;
        if (IsOfKind<Cfht12K>(mainHeader)||IsOfKind<Megacam>(mainHeader))
          {
          file_name_format = "%s%s%02d.%s";
          }
        else
          {
          file_name_format = "%s%s_%d.%s";
          }
        sprintf(outFileName,file_name_format,outDirName.c_str(),
		                            base.c_str(),
                                            chip, ext.c_str());
	cout << "   chip " << chip << " -> " << outFileName << endl;
	FitsHeader outFile(inFile, outFileName);
        outFile.Append_LowPriority(mainHeader);
        outFile.AddHistoryLine("This header was patched from a main header and an image extension one");
        if (outFile.HasKey("EXTEND")) outFile.ModKey("EXTEND",false);
	inFile.CopyDataTo(outFile);
      }
  }
return 1;
}
