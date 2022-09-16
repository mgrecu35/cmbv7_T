      subroutine monortm_lut(freq_index, pavg, tavg, ravg, kext)

      USE define_vars

!     Look-up table version of MonoRTM

      implicit  none

      integer :: freq_index
      real    :: pavg
      real    :: tavg
      real    :: ravg
      real    :: kext

      integer :: i,j,k,n
      integer :: p1,p2
      integer :: t1,t2
      integer :: r1,r2
      real    :: pw1,pw2
      real    :: tw1,tw2
      real    :: rw1,rw2

      p1 = 1
      do n=1,npres
        if (pres(n) .gt. pavg) then
          p2 = n
          exit
        else
          p1 = n
          p2 = n
        endif
      enddo
      if (p1 .eq. p2) then
        pw1 = 0.5
        pw2 = 0.5
      else
        pw1 = (pres(p2) - pavg) / (pres(p2) - pres(p1))
        pw2 = (pavg - pres(p1)) / (pres(p2) - pres(p1))
      endif

      t1 = 1
      do n=1,ntemp
        if (temp(n) .gt. tavg) then
          t2 = n
          exit
        else
          t1 = n
          t2 = n
        endif
      enddo
      if (t1 .eq. t2) then
        tw1 = 0.5
        tw2 = 0.5
      else
        tw1 = (temp(t2) - tavg) / (temp(t2) - temp(t1))
        tw2 = (tavg - temp(t1)) / (temp(t2) - temp(t1))
      endif

      r1 = 1
      do n=1,nrmix
        if (rmix(n) .gt. ravg) then
          r2 = n
          exit
        else
          r1 = n
          r2 = n
        endif
      enddo
      if (r1 .eq. r2) then
        rw1 = 0.5
        rw2 = 0.5
      else
        rw1 = (rmix(r2) - ravg) / (rmix(r2) - rmix(r1))
        rw2 = (ravg - rmix(r1)) / (rmix(r2) - rmix(r1))
      endif

      kext = (rw1 * tw1 * pw1 * kabs(p1,t1,r1,freq_index)) + &
             (rw1 * tw1 * pw2 * kabs(p2,t1,r1,freq_index)) + &
             (rw1 * tw2 * pw1 * kabs(p1,t2,r1,freq_index)) + &
             (rw1 * tw2 * pw2 * kabs(p2,t2,r1,freq_index)) + &
             (rw2 * tw1 * pw1 * kabs(p1,t1,r2,freq_index)) + &
             (rw2 * tw1 * pw2 * kabs(p2,t1,r2,freq_index)) + &
             (rw2 * tw2 * pw1 * kabs(p1,t2,r2,freq_index)) + &
             (rw2 * tw2 * pw2 * kabs(p2,t2,r2,freq_index))

      return
      end
