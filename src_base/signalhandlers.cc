#include <signal.h>
#include <stdlib.h> // for abort

#include "polokaexception.h"


/* convert (some) signals into exceptions to allow for:
   - post mortem diagnostic (both in logs and get a core if needed).
   - exception catching that does some proper cleanup in destructors
      (especially fits files)
*/



class SignalException : public PolokaException
{
 private:
  int copy_count; // required because exceptions are copied without notice.
  bool should_abort; // triggers abort to get a core.

 public :
  SignalException(const string &Mess, const bool SA) 
    : PolokaException(Mess), copy_count(0), should_abort(SA) 
  {};
  
  SignalException(const SignalException &E) : PolokaException(E) 
  {
    copy_count = E.copy_count+1;
    should_abort = E.should_abort;
  }

  ~SignalException()
  {
    if (should_abort && (copy_count == 0)) abort();
  }
};




static void my_signal_handler(int Signum)
{
  switch (Signum)
    {
    case SIGUSR1 :
      throw(SignalException(" signal USR1 received: probably BQS Time Limit Reached ", false));
      break;
    }
}


// install handler function

static int install_signal_handler()
{
  signal(SIGUSR1, my_signal_handler);
  return 1;
}

// call it... implicitely

static int toto = install_signal_handler();
