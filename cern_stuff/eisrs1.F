*
* $Id: eisrs1.F,v 1.1 2004/02/20 10:48:16 nrl Exp $
*
* $Log: eisrs1.F,v $
* Revision 1.1  2004/02/20 10:48:16  nrl
* Initial revision
*
* Revision 1.1  2003/10/13 14:31:32  guy
* CERNLIB functions for matrix diagonalization
*
* Revision 1.1.1.1  1996/04/01 15:02:39  mclareni
* Mathlib gen
*
*
#include "kernnum/pilot.h"
      SUBROUTINE EISRS1(NM,N,AR,WR,ZR,IERR,WORK)
C     ALL EIGENVALUES AND CORRESPONDING EIGENVECTORS OF A REAL
C     SYMMETRIC MATRIX
C
      DIMENSION AR(NM,NM),WR(N),ZR(NM,NM),WORK(1)
      CALL TRED2(NM,N,AR,WR,WORK,ZR)
      CALL TQL2(NM,N,WR,WORK,ZR,IERR)
      RETURN
      END
