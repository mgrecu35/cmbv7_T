      MODULE define_vars

      implicit none

      integer              :: nfreq
      integer              :: nchan
      integer              :: npres
      integer              :: ntemp
      integer              :: nrmix

      integer, allocatable :: ifreq(:)
      integer, allocatable :: ipol(:)
      real, allocatable    :: freq(:)
      real, allocatable    :: pres(:)
      real, allocatable    :: temp(:)
      real, allocatable    :: rmix(:)
      real, allocatable    :: kabs(:,:,:,:)

      END MODULE define_vars
