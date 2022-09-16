      program test_read_melt_percentages

    ! test the read_melt_percentages subroutine for reading tables
    ! of lwc melt and rate fractions as functions of depth in the
    ! melting layer and d0 value

    ! Bill Olson  July, 2013

      use Tables_frac

      implicit none

      integer :: i, j, k, kk

      real :: heightBB, heightML
      real :: mu, d0, height, mlwc_frac, mrate_frac

      call read_melt_percentages

    ! write out fractions
      k = 1
      write(6, '(/,"for mu: ", f10.4)') mu_init + float(k - 1) * mu_interval
      write(6, '(/, "lwc melt fraction: ")')
      do i = 1, d0_nx
        write(6, '(f10.6, 2x, i5, 2x, 51f5.2)') &
         d0_init + float(i - 1) * d0_interval, max_j(i, k), (mlwc_fraction(i, j, k), j = 1, height_ny)
      end do

      write(6, '(/, "rate melt fraction: ")')
      do i = 1, d0_nx
        write(6, '(f10.6, 2x, i5, 2x, 51f5.2)') &
         d0_init + float(i - 1) * d0_interval, max_j(i, k), (mrate_fraction(i, j, k), j = 1, height_ny)
      end do

    ! test values
      heightBB = 500.
      heightML = 1000.
      mu = -2.
      d0 = 0.00
      height = 100.

!      do kk = 1, 700

!        d0 = 0.0 + float(kk - 1) * d0_interval

        call interp_melt_percentages(heightBB, heightML, mu, d0, height, mlwc_frac, mrate_frac)

        write(6, '("input heightBB: ", f10.4, "  heightML: ", f10.4, "  mu: ", f10.4, "  d0: ", f10.4, &
         " and  height: ", f10.4, "  output mlwc_frac: ", f10.4, "  mrate_frac: ", f10.4)') &
         heightBB, heightML, mu, d0, height, mlwc_frac, mrate_frac

!      end do



      stop
      end
