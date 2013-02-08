#include <cstdlib>
#include <string>
#include <iostream>


#include <poloka/fileutils.h>
#include <poloka/toadscards.h>
#include <poloka/polokaexception.h>

using namespace std;



static string FileName;

void SetDatacardsFileName(const string &NewFileName)
{
  if (!FileExists(NewFileName))
    throw(PolokaException("SetDatacardsFileName : cannot find "+NewFileName));
  cout << " INFO : Setting the default DataCards file name to '" 
       << NewFileName << "'" << endl;
  FileName = NewFileName;
}


#define A_NAME "sub.datacard"

string DefaultDatacards()
{
  if (FileName != "") return FileName;
  if (FileExists(A_NAME)) return string(A_NAME);
  char *where = getenv("TOADSCARDS");
  if (where)
    {
      return (AddSlash(where)+A_NAME);
    }
  else // expect troubles
    return "";
}

    
