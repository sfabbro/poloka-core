#include <iostream.h>
#include <stdio.h>
#include <string>
#include <map.h>
#include <list.h>

#include "dbimage.h"
#include "fitsimage.h"
#include "fileutils.h"


static void dump_im_list(const DbImageList &List)
{
for ( DbImageCIterator it = List.begin(); it != List.end(); ++it)
  {
    cout << (*it).Name() << ' ' ;
  }
 cout << endl;
}

struct Mariage
{
string name;
DbImageList Ref;
DbImageList New;
  Mariage(const string &Name) : name(Name) {};
  Mariage() {};

  void dump() const { cout << name << endl;
  cout << " ref " ; dump_im_list(Ref); 
  cout << " new " ; dump_im_list(New);};
};


string FieldCamName(const string &Object, int CcdId)
{
  /* pick up numbers at the end of the CAT-NAME field  (Object here) */
const char *pstart = Object.c_str() + Object.length() - 1;
RemovePattern(Object," ");
while (pstart > Object.c_str() && isdigit(*pstart)) pstart--;
char toto[12]; sprintf(toto,"%d",CcdId); 
//return "field"+ string(pstart+1) + "ccd" + string(toto);
return "field"+ Object+ "ccd" + string(toto);
}

// what comes just after is really specific for each run.
// rather than using complicated scripting we resort to
// hand modifications of the 2 next routines.

string Translate(const string &RawName)
{
if (strstr(RawName.c_str(),"CFDF 0300+00")) return "P1";
if (strstr(RawName.c_str(),"CFDF 2215+00")) return "P2";
return RawName;
}

string BuildName(const DbImage &Image)
{
string fileName = Image.FitsImageName(Calibrated);
if (!FileExists(fileName)) return "NotFlatFielded";
FitsHeader h(fileName);
if (!h.IsValid()) return "NotFlatFielded";

//string object = h.KeyVal("TOADOBJE");
//object = Translate(object);

string ra = h.KeyVal("TOADRASC");
string dec = h.KeyVal("TOADDECL");
string object = ra + dec;
int ccd = h.KeyVal("TOADCHIP");
return FieldCamName(object, ccd);
}


typedef map<string,Mariage> MariageMap;
typedef map<string,Mariage>::iterator MariageIterator;
MariageMap Mariages;


void Weddings(DbImageList &Ref, DbImageList &New)
{
  cout << "NREf = " << Ref.size() << " NNew =  " << New.size() << endl;
for (DbImageIterator itRef = Ref.begin(); itRef != Ref.end(); ++itRef)
  {
  string a_name = BuildName(*itRef);
  Mariages[a_name].Ref.push_back(*itRef);
  Mariages[a_name].name = a_name;
  }
for (DbImageIterator itNew = New.begin(); itNew != New.end(); ++itNew)
  {
  string a_name = BuildName(*itNew);
  Mariages[a_name].New.push_back(*itNew);
  Mariages[a_name].name = a_name;
  }
 cout << " covered  fields ================================================== " << endl;
for (MariageIterator mit = Mariages.begin(); mit != Mariages.end(); mit++)
  {
  (*mit).second.dump();
  cout << " --------------------------" << endl;
  }
}


int ImageReadyForWedding(const DbImage &image)
{
int ok = 1;
if (!FileExists(image.FitsImageName(Calibrated)))
  {
   cerr << " image " << image.Name() << " misses the fits calibrated image file -> ignored " << endl;
   ok = 0;
  }
if (!FileExists(image.ImageCatalogName(SExtractor)))
  {
   cerr << " image " << image.Name() 
	<< " misses the catalog -> ignored " 
	<< image.ImageCatalogName(SExtractor)
	<< endl;
   ok = 0;
  }
return ok;
}
    



int ImageListReady(const DbImageList &List)
{
for (DbImageCIterator it = List.begin(); it != List.end(); ++it)
  {
  if (!ImageReadyForWedding(*it)) return 0;
  }
return 1;
}


int MinRefSize = 1;
int MinNewSize = 2;


int MariageReady(const Mariage &AMariage)
{
int ok = 1;
if (int(AMariage.Ref.size()) < MinRefSize || int(AMariage.New.size()) < MinNewSize )
  {
  return 0;
  }
ok *= ImageListReady(AMariage.Ref) * ImageListReady(AMariage.New);
return ok;
}

void MariageListFilter()
{
for (MariageIterator mit = Mariages.begin(); mit != Mariages.end(); )
  {
  Mariage &mariage =(*mit).second ;
  if (!(MariageReady(mariage)))
    {
    MariageIterator next = mit;
    ++next;
    Mariages.erase(mit); /* le erase tout seul est bugge me semble-t-il ... P.A. */
    mit = next;
    }
  else ++mit;
  }
 cout << "\n\n ============ "<< Mariages.size() << " VIABLES SUBTRACTIONS =============================== " << endl;
for (MariageIterator mit = Mariages.begin(); mit != Mariages.end(); ++mit )
  {
  (*mit).second.dump();
  cout << " --------------------------" << endl;
  }  
}


static bool IsInDbImageList(const DbImageList &List, const string &DbImageName)
{
for (DbImageCIterator it = List.begin() ; it != List.end(); ++it)
  {
  if ((*it).Name() == DbImageName) return true;
  }
return false;
}

static void ImageListDump(const DbImageList &List, ostream &stream = cout)
{
for (DbImageCIterator it = List.begin(); it != List.end(); ++it)
  {
    stream << (*it).Name() << ' ' << BuildName(*it) << endl;
  }
}

Mariage* WhichMariage(DbImage & image)
{
for (MariageIterator mit = Mariages.begin(); mit != Mariages.end(); ++mit )
  {
  Mariage &mariage = (*mit).second;
  if (IsInDbImageList(mariage.Ref, image.Name()) || IsInDbImageList(mariage.New, image.Name()))
    {
    return &mariage;
    }
  }
return NULL;
}


void PendingImages(DbImageList Ref, DbImageList New) 
   /* no const nor & on purpose, we are going to delete entries... */
{
for (DbImageIterator itr = Ref.begin(); itr != Ref.end(); )
  {
    if (WhichMariage(*itr)) {itr = Ref.erase(itr);}
    else ++itr;
  }
for (DbImageIterator itr = New.begin(); itr != New.end(); )
  {
    if (WhichMariage(*itr)) {itr = New.erase(itr);}
    else ++itr;
  }
 cout << "\n\n ====== PENDING IMAGES ============" << endl;
 cout << " references : (" << Ref.size() << " images ) " << endl;
 ImageListDump(Ref);
 cout << " \n\n new images (" << New.size() << " images)" << endl;
 ImageListDump(New);
}

  


void MariageListMakeSubFiles(bool Overwrite, bool OnlyOneSub)
{
  cout << " \n\n ============== ATTEMPTING to create subtraction files =========================== " << endl;
char *sub =getenv("SUB");
if (!sub) sub = getenv("sub"); 
if (!sub)
  {
    cerr << " I was planning to write subtraction files in $SUB (or $sub) which is undefined " << endl;
    return;
  }
for (MariageIterator mit = Mariages.begin(); mit != Mariages.end(); ++mit )
  {
  Mariage &mariage = (*mit).second;
  string DirectoryName = AddSlash(string(sub))+mariage.name;
  string FileName = DirectoryName + "/subfile";
  if (!FileExists(DirectoryName))
    { 
    string command = "mkdir " + DirectoryName;
    if (system(command.c_str()) != 0)
      {
        cerr << " could not create directory " << DirectoryName << endl;
        continue;
      }
    }
  if (!Overwrite && FileExists(FileName)) continue;
  FILE *file = fopen(FileName.c_str(), "w");
  if (!file)
    {
      cerr << " cannot open " << FileName << " for writing " << endl;
      continue;
    }
  cout << " creating " << FileName << endl;
  fprintf(file,"# subtraction for %s\n",mariage.name.c_str());
  fprintf(file,"REF\n");
  for (DbImageIterator dbi = mariage.Ref.begin(); dbi != mariage.Ref.end(); ++dbi)
    {
     fprintf(file,"%s\n",(*dbi).Name().c_str());
    }
  fprintf(file,"NEW1\n");
  int half =  mariage.New.size()/2;
  int count =0;
  for (DbImageIterator dbi = mariage.New.begin(); dbi != mariage.New.end(); ++dbi)
    {
    if (count == half) fprintf(file,"NEW2\n"); 
    fprintf(file,"%s\n",(*dbi).Name().c_str());
    count++;
    }
  if (OnlyOneSub || getenv("ONLY_ONE_SUB"))
     fprintf(file,"ONLY ONE SUBTRACTION\n");
  fclose(file);
  }
}
  


static void usage()
{
  cerr << " mariages -R <references> -N <new1> -n (no subtraction files created)\n"
       << " where references and new are to tags in dbconfig " << endl 
       << " example : mariages -R '199903*' -N 19990410r 19990411r " << endl
       << "  -1 : generates subfiles for one sub scheme" << endl
       << "  -o : enable overwriting existing subfiles. " << endl
       ;
exit(1);
}



int main(int argc, char **argv)
{
bool no_sub_files = false;
if (argc < 5) usage();
string refName = "";
string newName = "";
DbImageList references;
DbImageList news;
bool onlyOneSub = false;
bool overwrite = false;
for (int i=1; i<argc; ++i)
  {
  char *arg = argv[i];
  switch (arg[1])
    {
    case 'R': i++; 
         while (argv[i][0] != '-' && i< argc) 
           {
           if (references.Collect(argv[i]) == 0)
             {
	       cerr << " No Images in " << arg[i] << endl;
             }
           ++i;
           }; 
         --i; break;
    case 'N': i++; 
         while (i<argc && argv[i][0] != '-') 
           {
           if (news.Collect(argv[i]) == 0)
             {
	       cerr << " No Images in " << arg[i] << endl;
             }
           ++i;
           }; 
         --i; break;
    case 'n': no_sub_files = true; break;
    case '1': onlyOneSub = true; break;
    case 'o': overwrite = true;
    default : cerr << " do not understand " << arg << endl; usage();
    }
  }
if (references.size() == 0) 
  {
    cerr << " no images in Reference list " << endl;
    usage();
  }

if (news.size() == 0) 
  {
    cerr << " no images in new list " << endl;
    usage();
  }
Weddings(references,news);
MariageListFilter();
PendingImages(references,news);
if (!no_sub_files) MariageListMakeSubFiles(overwrite, onlyOneSub);
return 1;
}
