*
* $Id: dsinv.F,v 1.1 2004/02/20 10:48:16 nrl Exp $
*
* $Log: dsinv.F,v $
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
* Revision 1.1.1.1  1996/02/15 17:49:05  mclareni
* Kernlib
*
*
#include "kernnum/pilot.h"
          SUBROUTINE          DSINV(N,A,IDIM,IFAIL)
          DOUBLE PRECISION    A(IDIM,*),  ZERO,  ONE,  X, Y
          REAL                PIVOTF
          CHARACTER*6         HNAME
          DOUBLE PRECISION    S1, S31, S32, S33,  DOTF
          PIVOTF(X)    =  SNGL(X)
          DOTF(X,Y,S1)  =  X * Y + S1
          DATA      HNAME               /  'DSINV '  /
          DATA      ZERO, ONE           /  0.D0, 1.D0 /
          IF(IDIM .LT. N  .OR.  N .LE. 0)  GOTO 900
#include "sfact.inc"
#include "sfinv.inc"
          RETURN
 900      CALL TMPRNT(HNAME,N,IDIM,0)
          RETURN
          END
