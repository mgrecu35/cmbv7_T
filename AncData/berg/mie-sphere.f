      SUBROUTINE MIE_SPHERE (X, MIN, QSCAT, QEXTI, ASYM, QBSCAT)
**
      implicit    none
      SAVE
**
**    Mie Routine P. Bauer 
**
      integer    limitx
      PARAMETER (LIMITX = 1500)

**
      REAL        X
      REAL        MR, MI, N1, N2
      REAL        QSCAT, QEXTI, QABSO, ASYM, QBSCAT
**
      REAL        RFAC1, RFAC2
      REAL        RHELP1(2), RHELP2(2)
**
      COMPLEX     M, MX, MIN
      COMPLEX     CHELP1, CHELP2, CFAC1, CFAC2, CBSCAT
**
      COMPLEX     DN(0:LIMITX), WN(-1:LIMITX)
      COMPLEX     AN(LIMITX), BN(LIMITX)
**
      INTEGER     NEND
      INTEGER     I100, I101
**
      EQUIVALENCE (CHELP1, RHELP1 (1))
      EQUIVALENCE (CHELP2, RHELP2 (1))
**
************************************************************************
**
      M      = CONJG (MIN)
      CHELP1 = M
      MR     =        RHELP1 (1)
      MI     = -1.0 * RHELP1 (2)
**      
      MX   = M  * X
      N1   = MR * X
      N2   = MI * X
**
      IF (X .LE. 20000.0) NEND = X + 4.00 * X ** (1.0 / 3.0) + 2.0
      IF (X .LE.  4200.0) NEND = X + 4.05 * X ** (1.0 / 3.0) + 2.0
      IF (X .LE.     8.0) NEND = X + 4.00 * X ** (1.0 / 3.0) + 1.0
      IF (NEND .LE.    5) NEND = 5
      IF (NEND .GT. LIMITX) NEND = LIMITX
**
      RFAC1      = SIN  (N1) * SIN  (N1) + SINH (N2) * SINH (N2)
      RHELP1 (1) = SIN  (N1) * COS  (N1) / RFAC1
      RHELP1 (2) = SINH (N2) * COSH (N2) / RFAC1
**
      DN (0) = CHELP1
**
      RHELP1 (1) =             COS (X)
      RHELP1 (2) = -1.0 E+00 * SIN (X)
      RHELP2 (1) =             SIN (X)
      RHELP2 (2) =             COS (X)
**
      WN (-1) = CHELP1
      WN ( 0) = CHELP2
**
      QEXTI  = 0.0
      QSCAT  = 0.0
      QBSCAT = 0.0
      QABSO  = 0.0
      ASYM   = 0.0 
      CBSCAT = CMPLX (0.0,0.0)
**
      DO 100 I100 = 1, NEND
         DN (I100) = -1.0 * I100 / MX
     +             +  1.0 / (I100 / MX - DN (I100 - 1))
         WN (I100) = WN (I100 - 1) * (2.0 * I100 - 1.0) / X
     +             - WN (I100 - 2)
**
         CFAC1 = DN (I100) / M + I100 / X
         CFAC2 = M * DN (I100) + I100 / X
**
         CHELP1 = WN (I100)
         CHELP2 = WN (I100 - 1)
**
         AN (I100) = (CFAC1 * RHELP1 (1) - RHELP2 (1))
     +             / (CFAC1 * CHELP1     - CHELP2    )
         BN (I100) = (CFAC2 * RHELP1 (1) - RHELP2 (1))
     +             / (CFAC2 * CHELP1     - CHELP2    )
**
         CHELP1 = AN (I100)
         CHELP2 = BN (I100)
**
         RFAC1 = RHELP1 (1) + RHELP2 (1)
         RFAC2 = CABS (AN (I100)) * CABS (AN (I100))
     +         + CABS (BN (I100)) * CABS (BN (I100))
**
         QEXTI  = QEXTI  + (2.0 * I100 + 1.0) * RFAC1
         QSCAT  = QSCAT  + (2.0 * I100 + 1.0) * RFAC2
         CBSCAT = CBSCAT + (2.0 * I100 + 1.0) * (-1.0) ** I100
     +          * (AN (I100) - BN (I100))
**
         IF (I100 .EQ. 1) GO TO 100
**
         CHELP1 = AN (I100 - 1) * CONJG (AN (I100))
     +          + BN (I100 - 1) * CONJG (BN (I100))
         CHELP2 = AN (I100 - 1) * CONJG (BN (I100 - 1))
**
         I101 = I100 - 1
         RFAC1  = I101 * (I101 + 2) / (I101 + 1.0)
         RFAC2  = (2.0 * I101 + 1.0) / (I101 * (I101 + 1.0))
**
         ASYM = ASYM + RFAC1 * RHELP1 (1) + RFAC2 * RHELP2 (1)
100   CONTINUE
**
      QEXTI  = QEXTI * 2.0 / (X * X)
      QSCAT  = QSCAT * 2.0 / (X * X)
      ASYM   = ASYM  * 4.0 / (X * X * QSCAT)
      QBSCAT = CABS (CBSCAT) * CABS (CBSCAT) / (X * X)
      IF (QSCAT .GT. QEXTI) QSCAT = QEXTI
**
      RETURN
      END
