// This may look like C code, but it is really -*- C++ -*-
#ifndef EXCEPTIONS_SEEN
#define EXCEPTIONS_SEEN


#define EXC_ABORT_NEG(_x)
#define EXC_ABORT_ALL(_x)

#define EXC_AWARE

#define END_CONSTRUCTOR

#define TRY try
#define CATCH(_var) catch(long _var)
#define CATCHALL catch(...)
#define ENDTRY

#define THROW(_i) throw((long) _i);

#define THROW_ THROW(0)

#define THROW_SAME throw;

#define ASSERT(_a_) if (!(_a_)) { \
     cerr << "Assertion failed " #_a_ " file " __FILE__ " line " << __LINE__ \
          << endl; \
     THROW_ }

#define RETURN(x) return(x)
#define RETURN_ return

void InitFailNewHandler();

#endif
