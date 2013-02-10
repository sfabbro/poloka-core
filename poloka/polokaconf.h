#ifndef POLOKACONF__H
#define POLOKACONF__H

#include <string>

const string POLOKA_CONF_DEFAULT_NAME="poloka.conf";

//! Where to read the datacards. if is either:
//! 1. what was set by SetDatacardsFileName
//! 2. $PWD/poloka.conf
//! 3. $TOADSCARDS/poloka.conf
string DefaultDatacards(const string& DefaultName=POLOKA_CONF_DEFAULT_NAME);

//! Set a new file name for datacards
void SetDatacardsFileName(const string& NewFileName);


#endif // POLOKACONF__H
