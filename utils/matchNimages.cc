#include <iostream>


#include "transformedimage.h"
#include "sestar.h"
#include "imagematch.h"
#include "listmatch.h"
#include "fileutils.h"
#include "fitsimage.h"
#include "imagesum.h"
#include "nstarmatch.h"

#include <vector>

void usage(char *progName)
{
  cerr << progName << " <images> " << endl;
}

#ifdef STORAGE
int ImagesStarMatch(const vector<string> &toMatch, const string &filename)
{
  
  vector<SEStarList> seStarListToMatch;
  
  DbImage refimage(toMatch[0]);
  SEStarList slref(refimage.ImageCatalogName());
  
  int n = toMatch.size();
  NStarMatchList nStrMtchList(n);
  
  for (int i=0; i<n ; ++i)
    {
      DbImage dbimage(toMatch[i]);
      SEStarList* seStrList = new SEStarList(dbimage.ImageCatalogName());

      BaseStarList* bsl2 = SE2Base(seStrList);

      CountedRef<Gtransfo> guess1,guess2;
      ImageListMatch(refimage,dbimage, guess1, guess2);
      nStrMtchList.SetTransfo(guess1,0,i);
      nStrMtchList.SetTransfo(guess2,i,0);
      nStrMtchList.MatchAnotherList(*bsl2,i,2.0,new SEStar());
    }
  
  //    int ndump = nStrMtchList.Nlists();
  int ndump = 200;
  int count = 0;
  for (NStarMatchIterator si = nStrMtchList.begin(); si != nStrMtchList.end(); ++si)
  {    
    for (int i =0; i<nStrMtchList.Nlists(); i++)
      {
	NStarMatch starMatch = *si;
	if (!(starMatch.StarExist(i)))
	  {
	    cout << "----  " ;
	    continue;
	  }
	//	SEStar *sestar = (SEStar *) starMatch.star[i];
	//	cout << sestar->flux <<"=";
        if (i==0) cout << starMatch.star[i]->x << ' ' << starMatch.star[i]->y << ' ' ;
	cout << starMatch.flux(i)<<"  ";
      }
    cout << endl;
    count++;
    if (count==ndump) break;
  }
  cout << " filename " << filename << endl;
  if (strlen(filename.c_str()) != 0) nStrMtchList.write(filename);
 return nStrMtchList.size();
}
#endif

int ImagesStarMatch(const vector<string> &toMatch, const string &filename)
{
  
  DbImage refimage(toMatch[0]);
  int n=toMatch.size();
  NStarMatchList  nStrMtchList(n);
  
  vector<SEStarList*> lists(n); 
   // it is extremely important that the lists to match
  //  live at least as long as the NStarMatchList

  for (int i=0; i<n ; ++i)
    {
      DbImage dbimage(toMatch[i]);
      string sexCatName = dbimage.ImageCatalogName(SExtractor);
      lists[i] = new SEStarList(sexCatName);
      BaseStarList* bsl2 = SE2Base(lists[i]);
  
      CountedRef<Gtransfo> guess1,guess2;
      ImageListMatch(refimage,dbimage, guess1, guess2);
      nStrMtchList.SetTransfo(guess1,0,i);
      nStrMtchList.SetTransfo(guess2,i,0);
      nStrMtchList.MatchAnotherList(*bsl2,i,2.0,new SEStar());
    }
  
  //    int ndump = nStrMtchList.Nlists();
  int ndump = 200;
  int count = 0;
  for (NStarMatchIterator si = nStrMtchList.begin(); si != nStrMtchList.end(); ++si)
  {    
    for (int i =0; i<nStrMtchList.Nlists(); i++)
      {
	NStarMatch starMatch = *si;
	if (!(starMatch.StarExist(i)))
	  {
	    cout << "----  " ;
	    continue;
	  }
	//	SEStar *sestar = (SEStar *) starMatch.star[i];
	//	cout << sestar->flux <<"=";
        if (i==0) cout << starMatch.star[i]->x << ' ' << starMatch.star[i]->y << ' ' ;
	cout << starMatch.flux(i)<<"  ";
      }
    cout << endl;
    count++;
    if (count==ndump) break;
  }
  cout << " filename " << filename << endl;
  if (strlen(filename.c_str()) != 0) nStrMtchList.write(filename);
 return nStrMtchList.size();
}



int main(int nargs, char **args)
{
string filename;
vector<string> toMatch;
for (int i=1; i<nargs; ++i)
  {
  char *arg = args[i];
  if (arg[0] != '-') 
    {
      DbImage dbimage(arg);
      if (!dbimage.IsValid())
	{
	  cerr << " Be careful ! " << arg << " must be an image name ! " << endl;
	}
      string sexCatName = dbimage.ImageCatalogName(SExtractor);
      if (!FileExists(sexCatName.c_str()))
	{
	  cerr << "The SExtractor Catalog associated to image " << arg << " doesn't exist !! " << endl;
	  continue;
	}
      toMatch.push_back(arg);
      continue;
    }
  else if (arg[1]=='f')
    {
      i++; filename = args[i];
    }
  }
ImagesStarMatch(toMatch,filename);
return EXIT_SUCCESS;
}
