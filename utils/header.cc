#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>



#include "fitsimage.h"
#include "dbimage.h"
#include "fileutils.h"


static void usage()
{
cout << " header  [-k <keyname> ... ]  <fitsFile..>" << endl;
cout << " or \n header -<what> [-k <keyname> ... ] <dbImage..> " << endl;
exit(EXIT_FAILURE);
}



static void fits_header_process(const string &FileName, vector<char *> requested_keys, const string &line_start)
{
  if (!FileExists(FileName)) { cerr << FileName << " does not exist" << endl; return;}
  //  if (!IsFits(FileName))     { cerr << FileName << " is not a fits file" << endl; return;}
  FitsHeader header(FileName);
  if (!header.IsValid()) return;
  int nkeys = requested_keys.size();
  if (nkeys ==0)
    {
      cout << header;
    }
  else 
    {
      cout << line_start;
      for (int j=0; j<nkeys; j++)
	{
	  if (header.HasKey(requested_keys[j])) 
	    cout << ' ' << header.KeyVal(requested_keys[j]);
	  else
	    cout << ' ' << requested_keys[j] << ": absent" ;
	}
      cout << endl;
    }
}

int main(int argc, char**argv)
{
vector<char*> requested_keys;
vector<string> fits_files; 
vector<string> names;
vector<string> which_info;


for (int i=1; i< argc; i++)
  {
  char *arg = argv[i];
  if (arg[0] != '-') 
    {
    names.push_back(arg);
    }
  else
    {
    if (strlen(arg)>2)   /* search for possible : -raw, -dead, -flat, -cal */
      {
      which_info.push_back(arg+1);
      }
    else
      {
      switch (arg[1])
        {
        case 'k' : i++ ; requested_keys.push_back(argv[i]); break;
        default : cerr << " do not understand " << arg << endl; usage();
        }
      }
    }
  }
if (argc <=1)
  {
  usage();
  }


int ninfo = which_info.size();

for (int in = 0 ; in < int(names.size()); in++)
  {
  if (ninfo == 0)
    {
    fits_header_process(names[in], requested_keys, names[in]);
    continue;
    }
  DbImage a_dbimage(names[in]);
  DbImageList list;
  if (!a_dbimage.IsValid())
    {
    list.Collect(names[in].c_str());
    }
  else
    {
    list.push_back(a_dbimage);   
    }
  for (DbImageCIterator dbi = list.begin(); dbi != list.end(); ++dbi)
    {
    const DbImage &dbimage = *dbi;
    for (int ii=0; ii<ninfo ; ++ii)
      {
	string filename = dbimage.GetFileName(which_info[ii].c_str());
	if (IsFits(filename.c_str()))
	  {
	    string tag = dbimage.Name();
	    if (ninfo != 1) tag = tag + "("+which_info[ii]+")";
	    fits_header_process(filename, requested_keys, tag);
	  }
      }
    }
  }
return EXIT_SUCCESS;
}




