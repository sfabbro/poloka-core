*
* $Id: dfinv.F,v 1.1 2004/02/20 10:48:16 nrl Exp $
*
* $Log: dfinv.F,v $
* Revision 1.1  2004/02/20 10:48:16  nrl
* Initial revision
*
* Revision 1.1.1.2  2002/10/03 11:35:46  seb
*
* first release directly imported from lpnhe paris
*
* Revision 1.1.1.1  1999/09/13 14:57:00  astier
*
*
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*
#include "kernnum/pilot.h"
          SUBROUTINE          DFINV(N,A,IDIM,IR)
          INTEGER             IR(*)
          DOUBLE PRECISION    A(IDIM,*),ZERO,     X, Y, TI
          DOUBLE PRECISION    S31, S32, S33, S34, DOTF
          CHARACTER*6         HNAME
          DATA      HNAME               /  ' DFINV'  /
          DOTF(X,Y,S31)  =  X*Y + S31
          DATA      ZERO      /  0.D0  /
#include "finv.inc"
          RETURN
          END
