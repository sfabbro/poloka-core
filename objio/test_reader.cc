// -*- C++ -*-
// 
// file test_reader.cc
// 
// 
#include <string>

#include "xmlstream.h"
#include "reader.h"
#include "readermgr.h"

struct Star {
  double x;
  double y;
  double flux;
};


int 
main()
{
  reader_mgr<Star,xmlistream> mgr;
  mgr.getClassReader(0);
}

