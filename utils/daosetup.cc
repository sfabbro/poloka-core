#include "allreducedimage.h"
#include "daophotio.h"

static void usage(char *progName) {
  cerr << progName << "[OPTION] DBIMAGE...\n"
       << "Setup a DBIMAGE to run DAOPHOT:\n"
       << " -d DIR: write in directory DIR instead of the dbimage directory\n"
       << " - write <dbimage>.fits with sky level added\n"
       << " - write daophot.opt, allstar.opt, photo.opt standard daophot option files\n"
       << " - write <dbimage>.als from se.list\n";
}

struct doDaoSetup {

  doDaoSetup() : Dir("") {}

  string Dir;

  void operator () (const string& DbImageName) const {
    ReducedImageRef im = ReducedImageNew(DbImageName);
    if (!im || !im->IsValid()) { 
      cerr << " not a valid dbimage: " << DbImageName << endl;
      return;
    }
    DaoSetup(*im, Dir);
  }
};

int main(int argc, char **argv) {

  if (argc < 2) {  usage(argv[0]); return EXIT_FAILURE; }

  list<string> imlist;
  doDaoSetup daosetup;

  for (int i=1; i<argc; ++i) {
    char *arg = argv[i];
    if (arg[0] != '-') {
      imlist.push_back(arg);
      continue;
    }
    if (arg[1] == 'd') { daosetup.Dir = argv[++i] ; continue; }
    cerr << " Error: unknown option: " << arg << endl;
    usage(argv[0]); 
    return EXIT_FAILURE;
  }

  for_each(imlist.begin(), imlist.end(), daosetup);

  return EXIT_SUCCESS;
}
