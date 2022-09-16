      subroutine read_lut(input_file)

      USE define_vars

!     Read routine for MonoRTM look-up table

      implicit  none

      character(len=100) :: input_file

      open(unit=10, file=trim(input_file), access='stream')

      read(10) nfreq
      allocate(freq(nfreq))
      read(10) freq
      !print*, nfreq
      read(10) nchan
      !print*, nchan
      allocate(ifreq(nchan))
      allocate(ipol(nchan))
      read(10) ifreq
      !print*, ifreq
      read(10) ipol
      !print*, ipol

      read(10) npres
      !print*, npres
      allocate(pres(npres))
      read(10) pres
      !print*, pres
      
      read(10) ntemp
      !print*, ntemp
      allocate(temp(ntemp))
      read(10) temp
      !print*, temp
      read(10) nrmix
      !print*, nrmix
      allocate(rmix(nrmix))
      read(10) rmix
      !print*, rmix
      allocate(kabs(npres,ntemp,nrmix,nfreq))
      read(10) kabs
      close(10)

      return
      end
