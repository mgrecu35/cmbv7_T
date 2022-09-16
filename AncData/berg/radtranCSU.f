      SUBROUTINE RADTRANCSU(nlyr, stype, view_angle, Tsfc, lyrhgt,  
     >                   lyrtemp, kabs, emis, refl, Tb, Tbdown)
C
C     Based on Kummerow PHD THESIS PLANE PARALLEL ATMOS, BUT
C     A) ASSUMING NO SCATTERING  
C     B) DOWNWELLING RADIANCE RECALCULATED AND CHANGED
C        BY GREG ELSAESSER (MAR 2005)

      implicit  none
      include 'parameters.inc'

      integer   nlyr, stype
      real      view_angle
      real      Tsfc
      real      lyrhgt(0:mxlyr), lyrtemp(0:mxlyr)
      real      kabs(mxlyr)
      real      emis, refl
      real      Tb, Tbdown

      real      umu, umu_down
      real      B0(mxlyr), B1(mxlyr)
      real      Z(0:mxlyr), Iout(0:mxlyr), Iin(mxlyr+1)

      real      XA, XB, YA, YB, dz
      real      TERM1, TERM2, TERM3, TERM4, TERM5
      integer   j
C                              
      umu = cos(view_angle*3.14159/180.0)
      Z(0) = lyrhgt(0)
      do j = 1, nlyr
        Z(j)  = lyrhgt(j)
        B0(j) = lyrtemp(j-1)
        B1(j) = (lyrtemp(j) - lyrtemp(j-1))/(lyrhgt(j) - lyrhgt(j-1))
        !print*, j, Z(j), B0(j), B1(j), Kabs(j)
      end do
      !stop
C     CALCULATE THE DOWNWELLING Tb THROUGH THE ATMOSPHERE AT UMU
      Iin(nlyr+1) = 2.7     ! Cosmic background radiation
C     LOOP THROUGH THE REMAINING LAYERS from the top down
      
      if(stype .gt. 1) then
        umu_down = 0.6
      else
        umu_down = umu
      endif
      !print*, view_angle, umu, umu_down
      do j = nlyr, 1, -1
C       CALCULATE RADIANCE FROM TOP OF LAYER "J"
        XA = B0(j) 
        XB = B1(j)
        YA = Kabs(j)/-umu
        dz = Z(j) - Z(j-1)

        TERM1 = Iin(j+1)*exp(YA*dz)
        TERM2 = XA*(1. - exp(YA*dz))
        if ( abs(ya*dz) .lt. 1.e-5 ) then
          TERM3 = 0.
        else
          TERM3 = XB/YA*(EXP(YA*dz) - YA*DZ - 1.)
        end if
        Iin(j) = TERM1 + TERM2 + TERM3
        !print*, j, Iin(j), TERM1, TERM2, TERM3
      end do
      Tbdown = Iin(1)
      
C
c     Do reflection from sfc and begin computing the upwelling through layers
      IOUT(0) = EMIS * Tsfc + REFL * Iin(1)
c      write(6,'("TB_down = ",f8.2,", emis = ",f8.2,", refl = ",f8.2)')
c     >          Iin(1),EMIS*Tsfc,Refl*Iin(1)

      DO J = 1, nlyr
C       CALCULATE THE UPWELLING RADIANCES AT THE TOP OF EACH LAYER J
        XA = B0(J) 
        XB = B1(J)
        YA = Kabs(j)/UMU
        DZ = Z(J) - Z(J-1)
 
        TERM1 = IOUT(J-1)*EXP(-YA*DZ)
        TERM2 = XA*(1. - EXP(-YA*DZ))
        if(abs(ya*dz) .lt. 1.e-5) then
          TERM3 = 0.
        else
          TERM3 = XB/YA*(YA*DZ - 1. + EXP(-YA*DZ))
        end if
        Iout(J) = TERM1 + TERM2 + TERM3 
      end do
      TB = Iout(nlyr)

      RETURN
      END
