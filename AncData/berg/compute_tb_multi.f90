      program compute_tb

      USE define_vars
      USE RSS_RTM

!     Code to compute Tb for a layered atmosphere including clouds and
!     hydrometeors as specified by user in MODEL.DAT

      implicit   none
      include   'parameters.inc'
    
      ! Input atmospheric model variables
 
      character(len=100) :: lut_file
      character(len=100) :: input_file
      character(len=10)  :: sat
      character(len=10)  :: sensor
      character(len=1)   :: cpol(2) = (/'v','h'/)

      logical, parameter :: diffuse   = .TRUE.
      integer, parameter :: verbose   = 0
      real,    parameter :: salinity  = 35.0
      real,    parameter :: mis_val = -999.9

      real               :: emis
      real               :: ebar
      real               :: refl
      real               :: view_angle
      real               :: azim
      real               :: e0(2), ewind(2), eharm(2,4), land_emis(13)
      real               :: omega(2), trans, optdepth

      integer :: n_arguments, iargc 
      integer :: nlyr, sfc_type
      logical :: file_test
      real    :: tskin, sfc_wind
      real    :: hgt_lev(0:mxlyr)
      real    :: press_lev(0:mxlyr)
      real    :: temp_lev(0:mxlyr)
      real    :: mix_ratio(mxlyr)
      real    :: cloud_water(mxlyr)
      
      ! Internally assigned/computed cloud parameters

      real    :: tavg, pavg, qavg, atm_ext, kext_clw, salb_clw, asym_clw
      real    :: kext(mxlyr)
      real    :: tb, tbdown
      real    :: tb_out(mxchan)
      real    :: tb_down(mxchan)
      real    :: cfreq(mxchan)

      ! Loop counters 

      integer :: ichan, ilyr, ifile, npol, i, k, n
      

      ! Read input file (atmospheric profile)

      !n_arguments = iargc()
      !if (n_arguments .eq. 2) then
      !  call getarg(1, input_file)
      !  call getarg(2, sensor)
      !else
      !  write(6,*) 'Error in command line arguments! '
      !  write(6,*) 'compute_tb <input_file> <sensor_id>'
      !  stop
      !endif
      lut_file = 'MonoRTM_v5.3-GMI.tbl'
      call read_lut(lut_file)
      cfreq(1:nchan) = freq(ifreq)
      !write(6,'(" Satellite: ",a)') sat
      !write(6,'(" Sensor:    ",a)') sensor
      !write(6,'(" Channels = ",15(f7.2,a1))') (cfreq(n),cpol(ipol(n)+1),n=1,nchan)
      !write(6,*)
      !loop over input files
      do ifile=0,360
        write(input_file,'(A,I3.3,A4)') '/PANFS/user/home/smunchak/data/gmi_input/merra2/MODEL.',ifile,'.DAT'
        INQUIRE(FILE=input_file, EXIST=file_test)
        if(file_test) then
          !print*, input_file
          open(unit=14, file=input_file, form='formatted')
          read(14,*) nlyr, sfc_type, hgt_lev(0), press_lev(0), temp_lev(0), tskin, sfc_wind, view_angle, azim
          read(14,*) land_emis
          !print '(13F8.5)',land_emis
          if (verbose .eq. 1) write(6,501) hgt_lev(0),press_lev(0),temp_lev(0),tskin,sfc_wind,view_angle,azim
          do k = 1, nlyr
            read(14,*) i, hgt_lev(k), press_lev(k), temp_lev(k), mix_ratio(k), cloud_water(k)
            if (verbose .eq. 1) write(6,502) k, hgt_lev(k), press_lev(k), temp_lev(k), mix_ratio(k), cloud_water(k)
          end do
          close(14)
          
          do ichan = 1, nchan
            if(ichan <= 9) then
              view_angle = 52.8
            else
              view_angle = 49.1
            endif
            optdepth = 0.0
            
            do ilyr = 1, nlyr   
              pavg = (press_lev(ilyr) - press_lev(ilyr-1)) / log(press_lev(ilyr)/press_lev(ilyr-1))
              tavg = (temp_lev(ilyr) + temp_lev(ilyr-1)) / 2.0
              if(ilyr > 1) then
                qavg = (mix_ratio(ilyr)+mix_ratio(ilyr-1)) / 2.0
              else
                qavg = mix_ratio(ilyr)
              endif
              call monortm_lut(ifreq(ichan), pavg, tavg, qavg, atm_ext)
              !call mie_clw(cfreq(ichan), tavg, cloud_water(ilyr), kext_clw, salb_clw, asym_clw)
          
              kext(ilyr) = atm_ext + kext_clw
              optdepth = optdepth + kext(ilyr) * (hgt_lev(ilyr) - hgt_lev(ilyr-1))
              !print*, pavg, tavg, mix_ratio(ilyr), atm_ext, hgt_lev(ilyr), optdepth
            end do 
            trans = exp(-optdepth)
          
            if((sfc_type .eq. 1) .or. (minval(land_emis) .le. 0)) then !ocean simulation
              call radtran(nlyr, sfc_type, view_angle, tskin, hgt_lev, temp_lev, kext, 0.5, 0.5, tb, tbdown)
              !write(6,*) cfreq(ichan)!, tbdown, tb
              call find_surface_tb(freq=cfreq(ichan), surtep=tskin, ssws=sfc_wind, tht=view_angle, &
                             phir=azim, sal=salinity, e0=e0, ewind=ewind, eharm=eharm, &
                             tran=trans, tbdw=tbdown, omega=omega)
              emis = e0(ipol(ichan)+1) + ewind(ipol(ichan)+1)
              if (diffuse) then
                refl = (1.0 - emis) * (1.0 + omega(ipol(ichan)+1))
              else
                refl = 1.0 - emis
              endif
              if (refl .lt. (1.0 - emis)) refl = 1.0 - emis ! Reflectivity cannot be less than 1.0 - emis
              if (refl .gt. 1.0) refl = 1.0                 ! Reflectivity cannot be greater than 1.0
            else
              emis = land_emis(ichan)
              refl = 1.0-emis
            endif
            call radtran(nlyr, sfc_type, view_angle, tskin, hgt_lev, temp_lev, kext, emis, refl, tb, tbdown)
            if ((tb .gt. 50.0) .and. (tb .lt. 350.0)) then
              tb_out(ichan)  = tb
              tb_down(ichan) = tbdown
            else
             write(6,*) ' Tb is outside physical range'
             print*, tb, tbdown, emis
             stop
            endif              
          end do  ! end loop over nchan

          ! write output file
          print '(I5,13F8.2)', ifile, tb_out(1:13)
          !write(6,*) 'Eddington Tb (MonoRTM-LUT)'
          !write(6,*) 
          !do ichan = 1, nchan
          !  write(6,512) cfreq(ichan), cpol(ipol(ichan)+1), tb_out(ichan)
          !end do
  501 format(f6.2,f8.2,f7.2,f8.2,3(f7.2))
  502 format(i3,f7.2,f8.2,2f7.2,f7.4)
 512  format(f6.2,a,': ',f9.2)
        endif
     end do
     stop
   end
