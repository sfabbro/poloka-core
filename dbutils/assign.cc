#include <iostream>
#include <list>

#include "fileutils.h"
#include "dbimage.h"
#include "fitsimage.h"
#include "stringlist.h"

struct Flat {
  int chip;
  string filter;
  string filename;
  Flat(const string file)
    {
      FitsHeader header(file);
      chip = header.KeyVal("TOADCHIP");
      filter = header.KeyVal("TOADBAND");
      filename = file;
    };
};


struct FlatList : public list<Flat> {
  Flat* locate(const int chip, const string filter)
    {
      for (iterator i=begin(); i!=end(); ++i)
	{
	  if ((*i).chip == chip && (*i).filter == filter)
	  return &(*i);
	}
      return NULL;
    }

  Flat* locate(const int chip)
    {
      for (iterator i=begin(); i!=end(); ++i)
	{
	  if ((*i).chip == chip)
	  return &(*i);
	}
      return NULL;
    }

  int addEntry(const string filename);
};

int FlatList::addEntry(const string filename)
{
  if (!FileExists(filename))
    {
      cerr << " Can not find file " << filename << " " << endl;
      return 0;
    }
  Flat flat(filename);
  Flat *p = locate(flat.chip, flat.filter);
  if (p)
    {
      cerr << " files " << filename << " and " << p->filename 
	   << " have the same filter and chip " << endl;
	exit(-1);
      return 0;
    }
  push_back(flat);
  return 1;
}

int process_one_image(const DbImage& dbimage, FlatList flatlist, char *WhichInfo)
{
  if (!dbimage.IsValid())
    {
      cerr << " cannot find image " << dbimage.Name() << " " << endl;
      return 0;
    }

  string image_file = dbimage.FitsImageName(Raw);
  if (!FileExists(image_file)) image_file = dbimage.FitsImageName(Calibrated);
  if (!FileExists(image_file)) image_file = dbimage.FitsImageName(Elixir);
  FitsHeader header(image_file);
  
  if (!header.IsValid())
    {
      cerr << " could not find an image file for DbImage " << dbimage.Name() << endl;
      return 0;
    }

  if (!header.HasKey("TOADCHIP"))
    {
      cerr << " file " << dbimage.Name() 
	   << " does not have a TOADCHIP Key in its header: nothing done " << endl;
      return 0;
    }

  const int chip = header.KeyVal("TOADCHIP");
  bool isbias = (strcmp(WhichInfo,"bias") == 0) || 
    (strcmp(WhichInfo,"dark") == 0) || (strcmp(WhichInfo,"dead") == 0);

  if (!header.HasKey("TOADBAND") && !isbias)
    {
      cerr << " file " << dbimage.Name() 
	   << " does not have a TOADBAND Key in its header: nothing done " << endl;
      return 0;
    }

  const string filter = header.KeyVal("TOADBAND");

  Flat* flat;
  if (isbias) flat = flatlist.locate(chip);
  else flat = flatlist.locate(chip, filter);
  
  
  if (!flat)
    {
      cerr << " NO " << WhichInfo << " found for image " << dbimage.Name() << endl;
      return 0;
    }
  return AssignInfo(dbimage, flat->filename, WhichInfo);
}



#ifdef OLD
/* there is one thing (at least which does not work on Linux with maps:
  there is no operator[] for const maps */
typedef map<int,string> IntStringMap;

int process_one_image(const DbImage& dbimage, IntStringMap &FlatMap, char *WhichInfo)
{
if (!dbimage.IsValid())
  {
    cerr << " cannot find image " << dbimage.Name() << " " << endl;
    return 0;
  }
const char *image_raw_file = dbimage.FitsImageName(Raw).c_str(); 
FitsHeader header(image_raw_file);
if (!header.IsValid())
  {
    cerr << " file " << image_raw_file << " does not seem to be a fits file " << endl;
    return 0;
  }

if (!header.HasKey("TOADCHIP"))
  {
    cerr << " file " << dbimage.Name() << " does not have a TOADCHIP Key in its header: nothing done " << endl;
    return 0;
  }
const int ccd = header.KeyVal("TOADCHIP");
string flatname = FlatMap[ccd];
if (FlatMap[ccd] == "")
  {
    cerr << " NO " << WhichInfo << " found for image " << dbimage.Name() << endl;
    return 0;
  }
return AssignInfo(dbimage, FlatMap[ccd], WhichInfo);
}


int IdentifyFlats(const StringList &FlatNames, IntStringMap& FlatMap)
{
for (StringCIterator fi = FlatNames.begin(); fi != FlatNames.end(); ++fi)
  {
  FitsHeader flat((*fi).c_str());
  if (!flat.IsValid())
    {
      cerr << " fits file " << *fi << " could not be open " << endl;
      FlatMap.clear();
      return 0;
    }
  int ccd = flat.KeyVal("TOADCHIP");
  if (FlatMap[ccd] != "")
    {
      cerr << *fi << " and " << FlatMap[ccd] << " both apply to ccd " << ccd <<  endl;
      FlatMap.clear();
      return 0;
    }
  FlatMap[ccd] = *fi;
  }
return 1;
}

int ReadList(char *fileName, StringList &Out)
{
ifstream list(fileName);
while (!list.eof())
  {
  string a_string;
  list >> a_string;
  Out.push_back(a_string);
  }
return Out.size();
}
#endif

static void usage()
{
  cerr << " Usage: assign -<what> -i <dbimage(s)> -f <fitsfile(s)> \n"
       << "  link a file to a dbimage, check filter and chip.\n"
       << "      <what>     : type of file to assign (in {raw,flat, bias,dead,fringe,calibrated})\n"
       << "      <dbimage(s)> : name of dbimages to assign \n"
       << "      <fitsfile(s)>    : the full path name of the FITS files to assign \n"
       << endl;
}

int main(int argc, char **argv)
{
if (argc < 6)
  {
    usage();
    exit(1);
  }


char *whichInfo = (argv[1]+1); /* suppress the leading '-' without check... */
  

StringList images;
StringList flats;
FlatList flatlist;
for (int i=2; i< argc; ++i)
  {
  char *arg = argv[i];
  if (arg[0] != '-') {cerr << " bad argument " << arg << endl; usage(); exit(1); }
  switch (arg[1])
    {
    case 'f' : 
      {
      for (++i ; i < argc ; ++i) 
        {
	if (argv[i][0] == '-') {--i; break;}
	flatlist.addEntry(argv[i]);
        flats.push_back(argv[i]);
        }
      break;
      }
    case 'i' :
      {
      for (++i ; i < argc ; ++i) 
        {
	if (argv[i][0] == '-') {--i; break;}
        images.push_back(argv[i]);
        }
      break;
      }
    default :
      {
	cerr << " bad argument " << arg << endl; usage(); exit(1);
      }
    }
  }

for (StringCIterator ci = images.begin(); ci != images.end(); ++ci)
  {
  DbImage dbimage((*ci).c_str());
  if (dbimage.IsValid()) 
    {
    process_one_image(dbimage, flatlist, whichInfo);
    continue;
    }
  /* try to interpret it as a path */
  DbImageList list((*ci).c_str());
  if (list.size() == 0)
    {
    cerr << " cannot interpret " << *ci << " as an image name nor a path " << endl;
    continue;
    }
  for (DbImageCIterator ii = list.begin(); ii != list.end(); ++ii)
    {
#ifdef OLD
    process_one_image(*ii, flatMap, whichInfo);
#endif
    process_one_image(*ii, flatlist, whichInfo);
    }
  }
return 0;
}
