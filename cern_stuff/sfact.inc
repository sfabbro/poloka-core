*
* $Id: sfact.inc,v 1.3 2010/09/08 22:09:11 seb Exp $
*
* $Log: sfact.inc,v $
* Revision 1.3  2010/09/08 22:09:11  seb
* Revert blindly all changes previous to Aug 29, 2010, which break for paris users
*
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
* Revision 1.1.1.1  1996/02/15 17:49:04  mclareni
* Kernlib
*
*
*
* sfact.inc
*
          IFAIL  =  0
          DO 144    J  =  1, N
             IF(PIVOTF(A(J,J)) .LE. 0.)  GOTO 150
             A(J,J)  =  ONE / A(J,J)
             IF(J .EQ. N)  GOTO 199
 140         JP1  =  J+1
             DO 143   L  =  JP1, N
                A(J,L)  =  A(J,J)*A(L,J)
                S1      =  -A(L,J+1)
                DO 141  I  =  1, J
                   S1  =  DOTF(A(L,I),A(I,J+1),S1)
 141               CONTINUE
                A(L,J+1)  =  -S1
 143            CONTINUE
 144         CONTINUE
 150      IFAIL  =  -1
          RETURN
 199      CONTINUE
