#include <stdlib.h>
#include <string>


#include "fileutils.h" // for FileExists()

#include "toadscards.h"


using namespace std;



static string FileName;

void SetDatacardsFileName(const string &NewFileName)
{
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

    
