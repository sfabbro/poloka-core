#include <iostream>
#include <fstream>

#include "fitsset.h"
#include "fileutils.h"
//#include "dbimage.h"
#include "fitsimage.h"
#include "fitstoad.h"
#include "frame.h"
FitsSet::FitsSet(const string &ListName, const bool CheckFilter, const bool CheckCCD)
{
vector<string> names;

nx =0;
ny =0;
nxtot = 0;
nytot = 0;

  
if (!FileExists(ListName)) 
  {
    cout << "file : " << ListName << " does not exist " << endl;
    return;
  }

ifstream list(ListName.c_str());

string firstFileName;
string filtRef;
string ccdRef;
string objRef;


while (!list.eof())
  {
    string fileName;
    list >> fileName;
    if (fileName.length() == 0) continue;
    FitsHeader loop(fileName);
    if (!loop.IsValid())
      {
        cout << " FitsSet :: cannot open " << fileName << endl;
        continue;
      }
    if (fileName.length()==0) continue;
    Frame Illu = TotalIlluRegion(loop);
    int nx_cur = int(Illu.Nx());
    int ny_cur = int(Illu.Ny());
    //    int x0_cur, y0_cur, nx_cur, ny_cur;
    //loop.IlluParameters(x0_cur, y0_cur, nx_cur, ny_cur);
    int nxtot_cur = loop.KeyVal("NAXIS1");
    int nytot_cur =loop.KeyVal("NAXIS2");
    string filt_cur = loop.KeyVal("TOADFILT");
    string ccd_cur = loop.KeyVal("TOADCHIP");
    string obj_cur = loop.KeyVal("TOADTYPE");

    if ((nx==0) && (ny==0) && (filtRef=="") && (ccdRef=="")) 
      {
        nx = nx_cur;
        ny = ny_cur;
        nxtot = nxtot_cur;
        nytot = nytot_cur;
        filtRef = filt_cur;
	ccdRef = ccd_cur;
	objRef = obj_cur;
	firstFileName = fileName;
      }
    if ((nx_cur != nx)||(ny_cur != ny))
      {
	cout << "FitsSet : Image " << fileName << " doesn't have same nx or ny as Ref Image "<< firstFileName << endl; 
	continue;
      }

    if (obj_cur != objRef)
      {
	cout << "FitsSet : Image " << fileName << " doesn't have the same object type as Ref Image "<< firstFileName << endl;
	//	continue;
      }

    if (CheckFilter && (filt_cur != filtRef))
      {
	cout << "FitsSet : Image " << fileName << " doesn't have same filter as Ref Image "<< firstFileName << endl; 
	continue;
      }

    if (CheckCCD && (ccd_cur != ccdRef))
      {
	cout << "FitsSet : Image " << fileName << " doesn't refer to the same CCD as Ref Image "<< firstFileName << endl; 
	continue;
      }

    if ((nxtot_cur != nxtot)||(nytot_cur != nytot))
      {
	cout << "FitsSet : Image " << fileName << " doesn't have same overscan as Ref Image "<< firstFileName << endl; 
        continue;
      }
  fitsNames.push_back(fileName);
  // could have something more clever depending on what the read line refers to (actual file or DbImage) 
  all_names += fileName;
  }

if ( ny==0 || nx==0)
  cout << " FitsSet::FitsSet: No Image opened from " << ListName << " --> empty set " << endl;


}


