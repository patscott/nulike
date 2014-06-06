***********************************************************************
*** nulike_aeinit initialises the IceCube effective area/volume and 
*** its errors, as well as the angular resolution.
***
*** input: filename       name of effective area file
***        nbins          number of energy bins in file
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: April 8, 2011
*** Modified: Jun 3, 6, 2014
***********************************************************************

      subroutine nulike_eainit(filename,nbins)

      implicit none
      include 'nulike.h'

      character (len=*) filename
      character (len=1) instring
      integer nbins, i, IER
      real*8 totalSystematic(max_nBinsEA), working(2*nbins-2)
      logical clusterdebug
      parameter(clusterdebug = .false.)

      !Open effective area file for reading
      open(lun,file=filename, ACTION='READ')

      !Skip header
      instring = '#'
      do while (instring .eq. '#')
        read(lun, fmt='(A1)'), instring
      enddo

      !Read in effective neutrino and anti-neutrino areas, uncertainties on eff area and angular resolutions
      do i = 1,nbins
        read(lun, *) instring, 
     &   effArea_logE(1,i,analysis), effArea_logE(2,i,analysis)
        !Convert to log(GeV) scale 
        effArea_logE(1,i,analysis) = dlog10(effArea_logE(1,i,analysis))
        effArea_logE(2,i,analysis) = dlog10(effArea_logE(2,i,analysis))
        !Find log-weighted bin centres 
        effArea_logEcentres(i,analysis) = 0.5d0*(effArea_logE(1,i,analysis)
     &   +effArea_logE(2,i,analysis))
        read(lun, *) instring, 
     &   effArea_nu(i), effArea_nubar(i)
        read(lun, *) instring, 
     &   effArea_syserr(i), effArea_staterr(i)
        read(lun, *) instring, 
     &   effArea_AngRes(i)
        read(lun, fmt='(A1)') instring
        if (i .ne. nbins) read(lun, fmt='(A1)'), instring
      enddo

      close(lun)

      !Calculate the percentage total systematic error from the effective 
      !area estimation in each bin, as the quadrature sum of the true
      !systematic and the Monte Carlo statistical error.  (Both contributions
      !have the character of a pure systematic when considered in the 
      !context of using a fixed effective area to do an analysis of real events.)
      totalSystematic = dsqrt(effArea_syserr*effArea_syserr + 
     & effArea_staterr*effArea_staterr) * 0.01d0      
      EAErr = maxval(totalSystematic)


      !Now need to init the interpolators in effective area and angular resolution.

      !Set up interpolation in neutrino effective area
      call TSPSI(nbins,effArea_logEcentres(:,analysis),effArea_nu,
     & 2,0,.false.,.false.,2*nbins-2,working,effArea_nuderivs,
     & effArea_nusigma,IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_eainit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up neutrino eff area.'
        stop
      endif

      !Set up interpolation in anti-neutrino effective area
      call TSPSI(nbins,effArea_logEcentres(:,analysis),effArea_nubar,
     & 2,0,.false.,.false.,2*nbins-2,working,effArea_nubarderivs,
     & effArea_nubarsigma,IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_eainit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up anti-nu eff area.'
        stop
      endif

      !Set up interpolation in angular resolution
      call TSPSI(nbins,effArea_logEcentres(:,analysis),effArea_AngRes,
     & 2,0,.false.,.false.,2*nbins-2,working,effArea_AngResderivs,
     & effArea_AngRessigma,IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_eainit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up angular resolution.'
        stop
      endif

      end subroutine nulike_eainit
