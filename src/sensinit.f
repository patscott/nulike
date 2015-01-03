***********************************************************************
*** nulike_sensinit initialises the IceCube effective area/volume and 
*** its errors, as well as the angular resolution.
***
*** input: filename       name of effective area/volume file
***        nbins          number of energy bins in file
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: April 8, 2011
*** Modified: Jun 3, 6, 15 2014
***********************************************************************

      subroutine nulike_sensinit(filename,nbins)

      implicit none
      include 'nulike.h'

      character (len=*) filename
      character (len=1) instring
      integer nbins, i, IER
      real*8 totalSystematic(max_nSensBins), working(2*nbins-2)
      logical clusterdebug
      parameter(clusterdebug = .false.)

      !Save bin number for other routines
      nSensBins(analysis) = nbins

      !Open effective area/volume file for reading
      open(lun,file=filename, ACTION='READ')

      !Skip header
      instring = '#'
      do while (instring .eq. '#')
        read(lun, fmt='(A1)'), instring
      enddo
      !Skip density block if it exists
      if (instring .eq. 'r') then
        read(lun, fmt='(A1)'), instring
        read(lun, fmt='(A1)'), instring
      endif

      !Read in effective neutrino and anti-neutrino areas/vols, uncertainties on eff area/vol and angular resolutions
      do i = 1,nbins
        read(lun, *) instring, 
     &   sens_logE(1,i,analysis), sens_logE(2,i,analysis)
        !Convert to log(GeV) scale 
        sens_logE(1,i,analysis) = dlog10(sens_logE(1,i,analysis))
        sens_logE(2,i,analysis) = dlog10(sens_logE(2,i,analysis))
        !Find log-weighted bin centres 
        sens_logEcentres(i,analysis) = 0.5d0*(sens_logE(1,i,analysis)
     &   +sens_logE(2,i,analysis))
        read(lun, *) instring, 
     &   sens_nu(i,analysis), sens_nubar(i,analysis)
        read(lun, *) instring, 
     &   sens_syserr(i,analysis), sens_staterr(i,analysis)
        read(lun, *) instring, 
     &   sens_AngRes(i,analysis)
        read(lun, fmt='(A1)') instring
        if (i .ne. nbins) read(lun, fmt='(A1)'), instring
      enddo

      close(lun)

      !Fix up the end bins
      if (nbins .lt. max_nSensBins) sens_AngRes(nbins+1:,analysis) = sens_AngRes(nbins,analysis)

      if (any(sens_AngRes(:,analysis) .lt. 1.d2*epsilon(0.d0))) then
        write(*,*)
        write(*,*) "Error: your effective volume file contains mean angular"
        write(*,*) "error values approximately equal to or less than zero."
        write(*,*) "You cannot have zero mean angular errors.  Please fix this."
        write(*,*)
        stop
      endif

      !Set the minimum detectable energy (only used in 2014 partial likelihood calculation).
      min_detectable_logE = sens_logE(1,1,analysis)

      !Calculate the percentage total systematic error from the effective 
      !area/volume estimation in each bin, as the quadrature sum of the true
      !systematic and the Monte Carlo statistical error.  (Both contributions
      !have the character of a pure systematic when considered in the 
      !context of using a fixed effective area/volume to do an analysis of real events.)
      totalSystematic = dsqrt(sens_syserr(:,analysis)**2 + 
     & sens_staterr(:,analysis)**2) * 0.01d0      
      EAErr(analysis) = maxval(totalSystematic)


      !Now need to init the interpolators in effective area/volume and angular resolution.

      !Set up interpolation in neutrino effective area/volume
      call TSPSI(nbins,sens_logEcentres(:,analysis),sens_nu(:,analysis),
     & 2,0,.false.,.false.,2*nbins-2,working,sens_nuderivs(:,analysis),
     & sens_nusigma(:,analysis),IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_sensinit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up nu/lepton effective area/vol.'
        stop
      endif

      !Set up interpolation in anti-neutrino effective area/volume
      call TSPSI(nbins,sens_logEcentres(:,analysis),sens_nubar(:,analysis),
     & 2,0,.false.,.false.,2*nbins-2,working,sens_nubarderivs(:,analysis),
     & sens_nubarsigma(:,analysis),IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_sensinit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up anti-nu eff area/vol.'
        stop
      endif

      !Set up interpolation in angular resolution
      call TSPSI(nbins,sens_logEcentres(:,analysis),sens_AngRes(:,analysis),
     & 2,0,.false.,.false.,2*nbins-2,working,sens_AngResderivs(:,analysis),
     & sens_AngRessigma(:,analysis),IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_sensinit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up angular resolution.'
        stop
      endif

      end subroutine nulike_sensinit
