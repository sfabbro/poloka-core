class PixCoord {
public:
  PixCoord() {
    x=y=0;
  }
  PixCoord(int nx,int ny) {
    x = nx;
    y = ny;
  }
  int x;
  int y;
};

class mycluster : public std::list<PixCoord> {
public:
  mycluster(int x,int y,Image &image) {
    addpoint(x,y,image);
  }
  void setinimage(Image &image, Pixel value) const {
    std::list<PixCoord>::const_iterator point = begin();
    std::list<PixCoord>::const_iterator endpoint = end();
    for(;point!=endpoint;++point) {
      image(point->x,point->y)=value;
    }
  }
  void enlarge(Image &image) {
    std::list<PixCoord>::iterator point = begin();
    std::list<PixCoord>::iterator endpoint = end();
    for(;point!=endpoint;++point) {
      addpoint(point->x,point->y,image,false);
    }   
  }
  
private:
  
  void addpoint(int x,int y,Image &image,bool addthispoint = true) {
    if(size()>10000) {
      cerr << "cluster overflow" << endl;
      return;
    }
    if(image(x,y)==0)
      return;
    image(x,y) = 0;
    if(addthispoint)
      push_back(PixCoord(x,y));
    if(y>0) {
      addpoint(x,y-1,image);
      if(x>0) addpoint(x-1,y-1,image);
      if(x<image.Nx()-1)  addpoint(x+1,y-1,image);
    }
    if(x>0) addpoint(x-1,y,image);
    if(x<image.Nx()-1)  addpoint(x+1,y,image);
    if(y<image.Ny()-1) {
      addpoint(x,y+1,image);
      if(x>0) addpoint(x-1,y+1,image);
      if(x<image.Nx()-1)  addpoint(x+1,y+1,image);
    }
  }
};

static void findclusters(const Image &image, std::list<mycluster> &clusters) {
  Image newimage = image;
  for (int ix=0;ix< newimage.Nx();ix++)
    for (int iy=0;iy< newimage.Ny();iy++)
      if(newimage(ix,iy)>0) {
	//cout << "New cluster" << endl;
	mycluster cluster(ix,iy,newimage);
	clusters.push_back(cluster);
      }
}


static void saveclustersinimage( const std::list<mycluster> &clusters , Image &image) {
  for (int ix=0;ix< image.Nx();ix++)
     for (int iy=0;iy< image.Ny();iy++)
       image(ix,iy)=0;
  std::list<mycluster>::const_iterator cluster = clusters.begin();
   std::list<mycluster>::const_iterator endcluster = clusters.end();
   for(;cluster!=endcluster;++cluster) {
     cluster->setinimage(image,1);
   }
}
