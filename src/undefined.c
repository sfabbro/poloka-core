/*
I got
/cern/pro/lib/libpacklib.a(qnexte.o): In function `qnexte_':
qnexte.o(.text+0x32): undefined reference to `__sigjmp_save'
qnexte.o(.text+0x3f): undefined reference to `__setjmp'
/cern/pro/lib/libpacklib.a(cfstati.o): In function `cfstati_':
cfstati.o(.text+0x65): undefined reference to `_xstat'
*/

/* in principle qnexte is useless for us.
   in redhat 6.0 _xstat is a weak reference of libc.a */


#if (HOST != lpnp54)
void qnexte_() {};


void _xstat() {}; 
#endif
