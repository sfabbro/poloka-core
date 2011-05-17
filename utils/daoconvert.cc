#include "daophotio.h"

static void usage(char *progName) {
  cerr << progName
       << " <dbimage> <daofile>: transform <dbimage> star list to a daophot star list with proper header<\n"
       << "  -s <sexfile> <daofile>: transform <sexfile> to <daofile> (does not need dbimage))\n"
       << "  -r <daofile> <sexfile>: reverse (does not need dbimage)\n"
       << "  -m <daofile> <sexfile>: merge daophot star list into a sextractor one (with star match)\n";
}

int main(int argc, char **argv) {

  if (argc < 3) {  usage(argv[0]); return EXIT_FAILURE; }

  string src,dest;

  if (argv[1][0] != '-') {
    src = argv[1];
    dest = argv[2];
    Sex2Dao(src, dest);
    return EXIT_SUCCESS;
  }

  if (argv[1][0] == 's') {
    src = argv[2];
    dest = argv[3];
    Sex2Dao(src, dest);
    return EXIT_SUCCESS;
  }

  if (argv[1][1] == 'r') {
    src = argv[2];
    dest = argv[3];
    Dao2Sex(src, dest);
    return EXIT_SUCCESS;
  }

  if (argv[1][1] == 'm') {
    src = argv[2];
    dest = argv[3];
    MergeSexDao(dest, src);
    return EXIT_SUCCESS;
  }

  usage(argv[0]);
  return EXIT_FAILURE;
}
