      SUBROUTINE  PCKPSF  (ID, X, Y, M, S, INDEX, MAXN, FITRAD, PSFRAD,
     .     VARY, MAGFIL, STRFIL, NREQ, MAGLIM)
C
C=======================================================================
C
C Subroutine to read in any of the data files created by DAOPHOT and
C to select reasonable candidates for PSF stars:
C
C   (1) More than a critical distance from the edge of the frame, and
C   (2) Having no brighter star within a critical distance.
C
C              OFFICIAL DAO VERSION:  1991 April 18
C
C=======================================================================
C
      IMPLICIT NONE
      INTEGER MAXN
C
C MAXN is the maximum number of stars permitted in a data file.
C
      REAL X(MAXN), Y(MAXN), M(MAXN), S(MAXN)
      INTEGER ID(MAXN), INDEX(MAXN)
C
      REAL ABS
C
      CHARACTER*256 MAGFIL, STRFIL
      CHARACTER CASE*4
      REAL PSFRAD, FITRAD, RADSQ, DY, RADIUS, LOBAD, HIBAD, THRESH
      REAL AP1, PHPADU, READNS, FRAD, XYMIN, XMAX, YMAX, MAGLIM
      REAL VARY
      INTEGER I, J, N, NREQ, ISTAR, JSTAR, ISTAT, NL, NCOL, NROW
      INTEGER ITEMS, NSTAR, NVARY, NMIN
C
C      COMMON /FILNAM/ COOFIL, MAGFIL, PSFFIL, PROFIL, GRPFIL
C
      FRAD = 0.0
C
C-----------------------------------------------------------------------
C
C SECTION 1
C
C Get input file name, open the file, and read its header.
C
      NVARY = NINT(VARY)
      NVARY = (NVARY+2)*(NVARY+1)/2
      NMIN = MAX(1,NVARY)
C      CALL TBLANK
C      CALL GETNAM ('Input file name:', MAGFIL)
C      IF ((MAGFIL .EQ. 'END-OF-FILE') .OR. 
C     .     (MAGFIL .EQ. 'EXIT')) THEN
C         MAGFIL = ' '
C         RETURN
C      END IF
C
      CALL INFILE (2, MAGFIL, ISTAT)
      IF (ISTAT .NE. 0) THEN
         CALL STUPID ('Error opening input file '//MAGFIL)
         RETURN
      END IF
C
C      CALL GETDAT ('Desired number of stars, faintest magnitude:', X, 2)
C      IF (X(1) .LE. 0.5) THEN
C         CALL CLFILE (2)
C         RETURN
C      END IF
C      NREQ = NINT(X(1))
C      MAGLIM = AMIN1(X(2), 90.)
C
C Generate output file name and open the file.
C
C      FILE = SWITCH (MAGFIL, CASE('.lst'))
C      CALL GETNAM ('Output file name:', FILE)
      CALL OUTFIL (3, STRFIL, ISTAT)
      IF (ISTAT .NE. 0) THEN
         CALL STUPID ('Error opening output file '//STRFIL)
         CALL CLFILE (2)
         RETURN
      END IF
C
      NL = 0
      CALL RDHEAD (2, NL, NCOL, NROW, LOBAD, HIBAD, THRESH, AP1,
     .     PHPADU, READNS, FRAD)
      IF (NL .EQ. 0) NL = 1
      IF (NL .GT. 3) NL = 1
C
C Copy input file header to output file.
C
      ITEMS = 6
      IF (FRAD .GT. 0.) ITEMS = 7
      CALL WRHEAD (3, 3, NCOL, NROW, ITEMS, LOBAD, HIBAD, THRESH,
     .     AP1, PHPADU, READNS, FRAD)
      XYMIN = FITRAD + 1.
      XMAX = REAL(NCOL) - FITRAD
      YMAX = REAL(NROW) - FITRAD
      RADIUS = PSFRAD+FITRAD+2.
      RADSQ = RADIUS**2
C
C-----------------------------------------------------------------------
C
C SECTION 2
C
C Read the input file in star by star.
C
      I=0
 2000 I=I+1                                   ! Begin loop over stars
C
 2010 CALL RDSTAR (2, NL, ID(I), X(I), Y(I), M(I), S(I))
      IF (ID(I) .LT. 0) GO TO 2100            ! END-OF-FILE encountered
      IF (ID(I) .EQ. 0) GO TO 2010            ! Blank line encountered
      IF (M(I) .GT. 96.) M(I) = -M(I)
      IF (I .LT. MAXN) GO TO 2000
C
      CALL STUPID ('*** WARNING ***  Too many stars in input file.')
      WRITE (6,6) MAXN
    6 FORMAT (I10, ' stars have been read.  I will work with these.')
      I = I+1
C
C Perform the selection.
C
 2100 NSTAR=I-1                                      ! Number of stars
      CLOSE (2)
      CALL QUICK (M, NSTAR, INDEX)
      I = INDEX(1)
      IF ((M(1) .GT. -90.) .AND. (X(I) .GE. XYMIN) .AND.
     .     (Y(I) .GE. XYMIN) .AND. (X(I) .LE. XMAX) .AND.
     .     (Y(I) .LE. YMAX)) THEN
         IF (S(I) .LT. 9999.999) THEN
            WRITE (3,321) ID(I), X(I), Y(I), M(1), S(I)
  321       FORMAT (I6, 4F9.3)
         ELSE
            WRITE (3,322) ID(I), X(I), Y(I), M(1), S(I)
  322       FORMAT (I6, 3F9.3, F9.2)
         END IF
         N = 1
         IF (N .GE. NREQ) GO TO 2200
      ELSE
         N = 0
      END IF
C
      DO 2195 ISTAR=2,NSTAR
         I = INDEX(ISTAR)
         IF ( (M(ISTAR) .LT. -90.) .OR. (X(I) .LT. XYMIN) .OR.
     .         (Y(I) .LT. XYMIN) .OR. (X(I) .GT. XMAX) .OR.
     .         (Y(I) .GT. YMAX) .OR. 
     .         ((M(ISTAR) .GE. MAGLIM) .AND. (N .GE. NMIN))
     .         ) GO TO 2195
C
         DO 2190 JSTAR=1,ISTAR-1
            J = INDEX(JSTAR)
            DY = ABS(Y(J)-Y(I))
            IF (DY .GE. RADIUS) GO TO 2190
            DY = DY**2 + (X(J)-X(I))**2
            IF (DY .LT. RADSQ) GO TO 2195
 2190    CONTINUE
C
         IF (S(I) .LT. 9999.999) THEN
            WRITE (3,321) ID(I), X(I), Y(I), M(ISTAR), S(I)
         ELSE
            WRITE (3,322) ID(I), X(I), Y(I), M(ISTAR), S(I)
         END IF
         N = N+1
         IF (N .GE. NREQ) GO TO 2200
 2195 CONTINUE
C
 2200 CLOSE (3)
      IF (N .EQ. 1) THEN
         WRITE (6,8) N
    8    FORMAT (/I10, ' suitable candidate was found.'/)
      ELSE
         WRITE (6,7) N
    7    FORMAT (/I10, ' suitable candidates were found.'/)
      END IF
      RETURN
C
      END!
