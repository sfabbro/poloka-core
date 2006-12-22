*
* $Id: eisrs1d.F,v 1.2 2006/12/22 13:35:40 guy Exp $
*
* $Log: eisrs1d.F,v $
* Revision 1.2  2006/12/22 13:35:40  guy
* merged aperdev branch (where the tag is prod07-0-in-aperdev)
*
* Revision 1.1.2.1  2006/02/21 14:27:29  astier
*
*   These routines are originally single precision and were translated manually to double precision.
*   To forbid to resolve them in the genuine cernlib (because it just does not work), I change their
*   names here.
*
*  Added Files:
*   Tag: aperdev
*  	eisrs1d.F tred2d.F tql2d.F
*  Removed Files:
*   Tag: aperdev
*  	eisrs1.F tred2.F tql2.F
*
* Revision 1.1.1.1  2004/02/20 10:48:16  nrl
* First version.
*
* Revision 1.1  2003/10/13 14:31:32  guy
* CERNLIB functions for matrix diagonalization
*
* Revision 1.1.1.1  1996/04/01 15:02:39  mclareni
* Mathlib gen
*
*
* cernlib routines translated to double precision
#include "kernnum/pilot.h"
      SUBROUTINE EISRS1D(NM,N,AR,WR,ZR,IERR,WORK)
C     ALL EIGENVALUES AND CORRESPONDING EIGENVECTORS OF A REAL
C     SYMMETRIC MATRIX
C
      DOUBLE PRECISION AR(NM,NM),WR(N),ZR(NM,NM),WORK(1)
      CALL TRED2D(NM,N,AR,WR,WORK,ZR)
      CALL TQL2D(NM,N,WR,WORK,ZR,IERR)
      RETURN
      END
