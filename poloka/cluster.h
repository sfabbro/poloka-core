#ifndef CLUSTER__H
#define CLUSTER__H

#include <list>
#include <poloka/image.h>



class ReducedImage;


class Cluster {

protected:
 public:
  long color ;
  long size ;
  double xsum, ysum ;
  double x2sum, xysum, y2sum ;


  const double Color() const {return color;};
  

  //! constructor
  Cluster(long c);

  long browse_and_color(const Image& src, long x0, long y0, 
			double threshold, Image& map);
  
} ;


class ClusterList : public list<Cluster>
{
 private:
  Image *colors;
  int debuglevel;
  bool trackfound;

 public:
  ClusterList(ReducedImage & inrim, int debuglevel=0);
  ~ClusterList(){delete colors;}
  bool Cut();
  Image Mask();
};

typedef ClusterList::const_iterator ClusterCIterator;
typedef ClusterList::iterator ClusterIterator;


#endif /*CLUSTER__H*/
