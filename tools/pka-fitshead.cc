#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

#include <poloka/fitsimage.h>
#include <poloka/dbimage.h>
#include <poloka/fileutils.h>
#include <poloka/polokaexception.h>


static void usage(const char* progname)
{
  cerr << "Usage: " << progname << " [OPTION]... FITS...\n"
       << "Usage: " << progname << " -what [OPTION]... DBIMAGE...\n"
       << "Display FITS header or FITS keys of a file\n\n"
       << "    -k FITSKEY: print FITSKEY\n"
       << "    -m STRING : print STRING when the key is missing (default: absent)\n"
       << "    -n        : do no print the filename\n"
       << "    -what     : for DBIMAGE uses (raw,cal,dead,weight,...)\n";
  exit(EXIT_FAILURE);
}



static void fits_header_process(const string &FileName, 
				vector<char *> requested_keys, 
				const string &line_start,
				const string &missing_key_placeholder)
{
  /* checking here that the file exists forbids to use "file[1]",
     which is a pity. If there is a problem FitsHeader::IsValid
     returns false.*/
  //  if (!FileExists(FileName)) { cerr << FileName << " does not exist" << endl; return;}
  //  if (!IsFits(FileName))     { cerr << FileName << " is not a fits file" << endl; return;}
  FitsHeader header(FileName);
  if (!header.IsValid()) return;
  size_t nkeys = requested_keys.size();
  if (nkeys ==0)
    {
      cout << header;
    }
  else 
    {
      cout << line_start;
      for (size_t j=0; j<nkeys; j++)
	{
	  if (header.HasKey(requested_keys[j])) 
	    cout << ' ' << header.KeyVal(requested_keys[j]);
	  else
	    if (missing_key_placeholder.empty())
	      cout << ' ' << requested_keys[j] << ": absent" ;
	    else
	      cout << ' ' << missing_key_placeholder << ' ';
	}
      cout << endl;
    }
}

int main(int argc, char**argv)
{
  if (argc <=1) usage(argv[0]);
  
  vector<char*> requested_keys;
  vector<string> fits_files; 
  vector<string> names;
  vector<char*> which_info;
  bool no_file_names = false; 
  string missing_key_placeholder;

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
		case 'n' : no_file_names = true; break;
		case 'm' : missing_key_placeholder = argv[++i]; break;
		default : cerr << argv[0] << ": do not understand " << arg << endl; usage(argv[0]);
		}
	    }
	}
    }
  
  size_t ninfo = which_info.size();
  bool ok = true;
  for (size_t in = 0 ; in < names.size(); in++)
    {
      try {

	if (ninfo == 0)
	  {
	    fits_header_process(names[in], requested_keys, 
				no_file_names ? " " : names[in],
				missing_key_placeholder);
	    continue;
	  }
	DbImage a_dbimage(names[in]);
	DbImageList imList;
	if (!a_dbimage.IsValid())
	  {
	    imList.Collect(names[in].c_str());
	  }
	else
	  {
	    imList.push_back(a_dbimage);   
	  }
	for (DbImageCIterator dbi = imList.begin(); dbi != imList.end(); ++dbi)
	  {
	    for (size_t ii=0; ii<ninfo ; ++ii)
	      {
		string filename = dbi->GetFileName(which_info[ii]);
		if (IsFits(filename))
		  {
		    string tag = dbi->Name();
		    if (ninfo != 1) tag = tag + "(" + string(which_info[ii]) + ")";
		    if (no_file_names) tag = "";
		    fits_header_process(filename, requested_keys, tag, 
					missing_key_placeholder);
		  }
	      }
	  }
	
      } catch (PolokaException p){
	p.PrintMessage(cerr);
	ok = false;
      }
    }
  
  return ok? EXIT_SUCCESS : EXIT_FAILURE;
}




