#ifndef TOADSCARDS__H
#define TOADSCARDS__H

#include <string>

//! Where to read the datacards. if is either what was set by SetDatacardsFileName, or it is "sub.datcards" in the cwd, or $TOADSCARDS/sub.datacards. 
std::string DefaultDatacards();
void SetDatacardsFileName(const std::string &NewFileName);


#endif /* TOADSCARDS__H */ 
