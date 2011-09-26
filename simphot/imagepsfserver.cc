//#include "imagepsfserver.h"



#include <map>
#include <string>

class ImagePSF;

#include "reducedimage.h"
#include "imagepsf.h"

class ImagePSFServer
{
  std::map<std::string, const ImagePSF*> psfMap;
  typedef std::map<std::string, const ImagePSF*>::const_iterator const_iterator;
  typedef std::map<std::string, const ImagePSF*>::iterator iterator;

 public:
  const ImagePSF* FindPSF(const ReducedImage* RI)
  {
    string key = RI->Dir();
    const_iterator i = psfMap.find(key);
    if (i != psfMap.end()) return i->second;
    return NULL;
  }

  void StorePSF(const ReducedImage* RI, const ImagePSF* Ipsf)
  {
    psfMap.insert(pair<string,const ImagePSF*>(RI->Dir(), Ipsf));
  }

  ~ImagePSFServer()
  {
    //    cout << " deleting ImagePSF's" << endl;
    for (iterator i = psfMap.begin(); i != psfMap.end(); ++i) delete i->second;
  }

};



static ImagePSFServer TheImagePSFServer;

const ImagePSF* FindImagePSF(const ReducedImage *Ri)
{
  const ImagePSF* i = TheImagePSFServer.FindPSF(Ri);
  if (!i) 
    {
      i = new ImagePSF(*Ri, false);
      TheImagePSFServer.StorePSF(Ri, i);
    }
  return i;
}

