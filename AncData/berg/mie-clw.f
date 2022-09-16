      subroutine mie_clw(freqy, temp, lwc, ksca, asca, gsca)
      

c     Compute the extinction, absorption, asymmetry parameter and
c     backscatter for a given water content of cloud water in [g/m^3]. 
c     Code assumes that cloud water drops are mono-disperse with an 
c     effective radius reff supplied in the parameter file.

c     Input:
c     freqy		frequency of radiation [GHz]
c     temp		temperature of particles [K]
c     lwc		water content of cloud water distribution [g/m**3]
c     reff              effective radius of particles [mm]

c     Output:
c     ksca		extinction coefficient [1/km]
c     asca		single-scatter albedo []
c     gsca		asymmetry factor []
c     pbck		backscatter phase function/(4*pi) []
c

      implicit none

      real    freqy, temp, lwc, reff_cloudwater
      real    ksca, asca, gsca, pbck

      real    pi, wavel, x
      real    densliq, density
      real    rad, dropmass
      real    num
      real    qext, qsca, asym, qbsca
      real    bext, bsca, bsym, bq11
      
      real     epsreal, epsimag
      complex  ewat, cref

c      
c     Assign some useful constants
      data       pi /3.14159265/
      data       densliq /1.0e+3/                 ! kg/m^3
      data       reff_cloudwater  / 0.10 /        ! [mm]
  
      wavel = 300./freqy
c
c     Begin by checking if hydrometeors of this species are present.
c     If not, set scattering parameters to zero and return.
c
      if(lwc .lt. 0.001) then
        ksca=0.
        asca=0.
        gsca=0.
        pbck=0.
        return
      endif
c
c     Initialize the scattering parameters
c
      bext=0.
      bsca=0.
      bsym=0.
      bq11=0.
c
c     Compute scattering properties
c
      rad = reff_cloudwater
      density = densliq
      dropmass = 4./3.*pi*density*rad*rad*rad*1.E-09
      num = lwc/dropmass      
c
c     Get complex refractive index of liquid water
c
      call watoptic(freqy,temp,0.0,epsreal,epsimag)
      ewat = cmplx(epsreal,epsimag)
      cref = csqrt(ewat)
c
c     call Mie program
c
      x = 2.*pi*rad/wavel 
      call mie_sphere(x,cref,qsca,qext,asym,qbsca)

      bext=num*qext*pi*rad*rad*1.e-6
      bsca=bsca+num*qsca*pi*rad*rad*1.e-6
      bsym=bsym+num*qsca*asym*pi*rad*1.e-6
      bq11=bq11+num*qbsca*pi*rad*rad*1.e-6
c
c     check for distribution with very small extinction;
c     set parameters to zero to avoid numerical problems
c
      if( bext .gt. 1.e-6) then
        ksca=bext
        asca=bsca/bext
        gsca=bsym/bsca
        pbck=bq11/bsca
      else
        ksca=0.
        asca=0.
        gsca=0.
        pbck=0.
      end if

      return
      end
