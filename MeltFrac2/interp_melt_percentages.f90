      subroutine interp_melt_percentages(heightBB, heightML, &
       mu, d0, height, mlwc_frac, mrate_frac)

    ! Bill Olson  July, 2013

      use Tables_frac

      implicit none

      logical :: printflag
      integer :: ilayer
      integer :: int_mu, int_d0_lo, int_d0_hi, int_height

      real :: heightBB, heightML, mu, d0, height, mlwc_frac, mrate_frac
      real :: frac_layer, height_table
      real :: heightBB_table, heightML_table
      real :: d0_low, d0_high, height_low, height_high

      printflag = .false.

    ! find mu index (assumes integral value)
      int_mu = nint((mu - mu_init) / mu_interval) + 1

    ! find d0 interval for interpolation
      int_d0_lo = int((d0 - d0_init) / d0_interval) + 1
      int_d0_hi = int_d0_lo + 1

      if(int_mu > mu_nz) int_mu = mu_nz
      if(int_mu < 1) int_mu = 1
      if(int_d0_lo > d0_nx - 1) then
        int_d0_lo = d0_nx
        int_d0_hi = d0_nx
      end if
      if(int_d0_lo <  1) then
        int_d0_lo = 1
        int_d0_hi = 1
      endif

    ! depth relative to bright band maximum, and fractional depth
      if(printflag) write(6, '("ny: ", i5, "  height_init: ", f10.4, &
       "  height_int: ", f10.4, "  int_mu: ", i5, "  int_d0_lo: ", i5, &
       "  int_d0_hi: ", i5, "  maxj: ", i5)') &
       height_ny, height_init, height_interval, int_mu, int_d0_lo, &
       int_d0_hi, max_j(int_d0_lo, int_mu)

      heightBB_table = height_init + float(max_j(int_d0_lo, int_mu) - 1) * height_interval
      heightML_table = height_init + float(height_ny - 1) * height_interval
      if(height < heightBB) then
        ilayer = 0
        frac_layer = height / heightBB
        height_table = frac_layer * heightBB_table
      else
        ilayer = 1
        frac_layer = (height - heightBB) / (heightML - heightBB)
        height_table = heightBB_table + frac_layer * (heightML_table - heightBB_table)
      endif
      int_height = int(height_table / height_interval) + 1

      if(int_height > height_ny - 1) int_height = height_ny - 1

      if(printflag) write(6, '("mu: ", f7.2, "  d0: ", f8.4, &
       "  height: ", f10.4, "  int_mu: ", i5, "  int_d0_lo: ", i5, &
       "  int_d0_hi: ", i5, "  ilayer: ", i5, "  heightBB_table: ", f10.4,  &
       "  heightML_table: ", f10.4, "  frac_layer: ", f8.4, &
       "  height_table: ", f10.6, "  int_height: ", i5)') mu, d0, height, int_mu, &
       int_d0_lo, int_d0_hi, ilayer, heightBB_table, heightML_table, &
       frac_layer, height_table, int_height

    ! bilinearly interpolate melt lwc and melt rate fractions
    ! height interpolation
      d0_low = d0_init + float(int_d0_lo - 1) * d0_interval
      d0_high = d0_init + float(int_d0_hi - 1) * d0_interval
      height_low = height_init + float(int_height - 1) * height_interval
      height_high = height_init + float(int_height) * height_interval

      if(printflag) write(6, '("d0_init: ", f8.4, &
       "  d0_interval: ", f10.6)') d0_init, d0_interval
      if(printflag) write(6, '("height_init: ", f8.4, &
       "  height_interval: ", f8.4)') height_init, height_interval
      if(printflag) write(6, '("d0_low: ", f10.6, &
       "  d0_high: ", f10.6)') d0_low, d0_high
      if(printflag) write(6, '("height_low: ", f10.6, &
       "  height_high: ", f10.6)') height_low, height_high

      if(printflag) write(6, '("mlwc_fractions: ", 4f10.4)') &
       mlwc_fraction(int_d0_lo,     int_height, int_mu), &
       mlwc_fraction(int_d0_hi,     int_height, int_mu), &
       mlwc_fraction(int_d0_lo, int_height + 1, int_mu), &
       mlwc_fraction(int_d0_hi, int_height + 1, int_mu)

      if(int_d0_hi > int_d0_lo) then

        mlwc_frac = &
         (mlwc_fraction(int_d0_lo,     int_height, int_mu) * (height_high - height_table) * (d0_high - d0) + &
          mlwc_fraction(int_d0_hi,     int_height, int_mu) * (height_high - height_table) * (d0 - d0_low) + &
          mlwc_fraction(int_d0_lo, int_height + 1, int_mu) * (height_table - height_low) * (d0_high - d0) + &
          mlwc_fraction(int_d0_hi, int_height + 1, int_mu) * (height_table - height_low) * (d0 - d0_low)) / &
         ((height_high - height_low) * (d0_high - d0_low))

      else

        mlwc_frac = &
         (mlwc_fraction(int_d0_lo,     int_height, int_mu) * (height_high - height_table)  + &
          mlwc_fraction(int_d0_lo, int_height + 1, int_mu) * (height_table - height_low)) / &
         (height_high - height_low)

      endif

      if(printflag) write(6, '("mrate_fractions: ", 4f10.4)') &
       mrate_fraction(int_d0_lo,     int_height, int_mu), &
       mrate_fraction(int_d0_hi,     int_height, int_mu), &
       mrate_fraction(int_d0_lo, int_height + 1, int_mu), &
       mrate_fraction(int_d0_hi, int_height + 1, int_mu)

      if(int_d0_hi > int_d0_lo) then

        mrate_frac = &
         (mrate_fraction(int_d0_lo,     int_height, int_mu) * (height_high - height_table) * (d0_high - d0) + &
          mrate_fraction(int_d0_hi,     int_height, int_mu) * (height_high - height_table) * (d0 - d0_low) + &
          mrate_fraction(int_d0_lo, int_height + 1, int_mu) * (height_table - height_low) * (d0_high - d0) + &
          mrate_fraction(int_d0_hi, int_height + 1, int_mu) * (height_table - height_low) * (d0 - d0_low)) / &
         ((height_high - height_low) * (d0_high - d0_low))

      else

        mrate_frac = &
         (mrate_fraction(int_d0_lo,     int_height, int_mu) * (height_high - height_table)  + &
          mrate_fraction(int_d0_lo, int_height + 1, int_mu) * (height_table - height_low)) / &
         (height_high - height_low)


      end if


    ! round up very high melt fractions to avoid extremely low ice concentrations
      if(mlwc_frac .ge. 0.99) mlwc_frac = 1.00
      if(mrate_frac .ge. 0.99) mrate_frac = 1.00

      return
      end subroutine interp_melt_percentages
