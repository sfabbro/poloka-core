#include <iostream>
#include <reducedimage.h>
#include <imagematch.h>
#include <gtransfo.h>


static void usage(const char* pg) {
  cerr << pg << " coordx coordy dbimage_from dbimage_to" << endl;
  exit(1);
}

int main (int argc, char **argv) {
  if (argc < 5)  usage(argv[0]);
  double coord_x = atof(argv[1]);
  double coord_y = atof(argv[2]);
  ReducedImage image_from(argv[3]);
  ReducedImage image_to(argv[4]);
  
  CountedRef<Gtransfo> direct,reverse;
  ImageListMatch(image_from,image_to,direct, reverse);
  
  
  Point P_from(coord_x,coord_y);
  Point P_to = direct->apply(P_from);
  
  cout << P_to.x << " " << P_to.y << endl;

  return 0;
}
