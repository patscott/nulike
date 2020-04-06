***********************************************************************
*** nulike_sensinit initialises the IceCube effective area/volume and
*** its errors, as well as the angular resolution.
***
*** input: filename       name of effective area/volume file
***        nbins          number of energy bins in file
***        like           likelihood version (2012 or 2015)
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: April 8, 2011
*** Modified: Jun 3, 6, 15 2014
***********************************************************************

      subroutine nulike_sensinit(filename,nbins,like)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      character (len=*) filename
      character (len=1) instring
      integer nbins, like, i, IER
      real*8 totalSystematic(max_nSensBins + 1), working(2*nbins)
      logical clusterdebug
      parameter(clusterdebug = .false.)

      !Make sure we understand what sort of file to expect
      if (like .ne. 2012 .and. like .ne. 2015) stop 'nulike_sensinit: unrecognised likelihood version.'

      !Save bin number for other routines.  Extra 'bin' is at low-E edge of first bin.
      nSensBins(analysis) = nbins + 1

      !Open effective area/volume file for reading
      open(lun,file=filename, ACTION='READ')

      !Skip header
      instring = '#'
      do while (instring .eq. '#')
        read(lun, fmt='(A1)') instring
      enddo
      !Skip density block if it exists
      if (instring .eq. 'r') then
        read(lun, fmt='(A1)') instring
        read(lun, fmt='(A1)') instring
      endif

      !Read in effective neutrino and anti-neutrino areas/vols, uncertainties on eff area/vol and angular resolutions
      do i = 2, nbins + 1
        read(lun, *) instring,
     &   sens_logE(1,i,analysis), sens_logE(2,i,analysis)
        !Convert to log(GeV) scale
        sens_logE(1,i,analysis) = dlog10(sens_logE(1,i,analysis))
        sens_logE(2,i,analysis) = dlog10(sens_logE(2,i,analysis))
        !Find log-weighted bin centres
        sens_logEcentres(i,analysis) = 0.5d0*(sens_logE(1,i,analysis)+sens_logE(2,i,analysis))
        if (like .eq. 2012) then
          read(lun, *) instring, sens_nu(i,analysis), sens_nubar(i,analysis)
        else
          read(lun, *) instring, sens_nu(i,analysis)
        endif
        read(lun, *) instring, sens_syserr(i,analysis), sens_staterr(i,analysis)
        read(lun, *) instring, sens_AngRes(i,analysis)
        read(lun, fmt='(A1)') instring
        if (i .ne. nbins + 1) read(lun, fmt='(A1)') instring
      enddo

      close(lun)

      !Make the effective volume go to zero at the lower edge of the first bin.
      sens_logE(1,1,analysis) = sens_logE(1,2,analysis)
      sens_logE(2,1,analysis) = sens_logE(1,2,analysis)
      sens_logEcentres(1,analysis) = sens_logE(1,2,analysis)
      sens_nu(1,analysis) = 0.d0
      if (like .eq. 2012) sens_nubar(i,analysis) = 0.d0
      !Set the first entry to have the same errors as the first bin.
      sens_syserr(1,analysis) = sens_syserr(2,analysis)
      sens_staterr(1,analysis) = sens_staterr(2,analysis)
      sens_AngRes(1,analysis) = sens_AngRes(2,analysis)
      !Set the entries above the uppermost bin to have the same errors as the edge bins.
      if (nbins .lt. max_nSensBins) then
        sens_syserr(nbins+2:,analysis) = sens_syserr(nbins+1,analysis)
        sens_staterr(nbins+2:,analysis) = sens_staterr(nbins+1,analysis)
        sens_AngRes(nbins+2:,analysis) = sens_AngRes(nbins+1,analysis)
      endif

      if (any(sens_AngRes(:,analysis) .lt. 1.d2*epsilon(0.d0))) then
        write(*,*)
        write(*,*) "Error: your effective volume file contains mean angular"
        write(*,*) "error values approximately equal to or less than zero."
        write(*,*) "You cannot have zero mean angular errors.  Please fix this."
        write(*,*)
        stop
      endif

      !Set the minimum detectable energy (only used in 2015 partial likelihood calculation).
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

      !Set up interpolation in neutrino effective area/CP-invariant effective volume
      call TSPSI(nbins+1,sens_logEcentres(:,analysis),sens_nu(:,analysis),
     & 2,0,.false.,.false.,2*nbins,working,sens_nuderivs(:,analysis),
     & sens_nusigma(:,analysis),IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_sensinit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up nu/lepton effective area/vol.'
        stop
      endif

      if (like .eq. 2012) then
        !Set up interpolation in anti-neutrino effective area
        call TSPSI(nbins+1,sens_logEcentres(:,analysis),sens_nubar(:,analysis),
     &   2,0,.false.,.false.,2*nbins,working,sens_nubarderivs(:,analysis),
     &   sens_nubarsigma(:,analysis),IER)
        if (IER .lt. 0) then
          write(*,*) 'Error in nulike_sensinit: TSPSI failed with error'
          write(*,*) 'code',IER,' in setting up anti-nu eff area/vol.'
          stop
        endif
      endif

      !Set up interpolation in angular resolution
      call TSPSI(nbins+1,sens_logEcentres(:,analysis),sens_AngRes(:,analysis),
     & 2,0,.false.,.false.,2*nbins,working,sens_AngResderivs(:,analysis),
     & sens_AngRessigma(:,analysis),IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_sensinit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up angular resolution.'
        stop
      endif

      end subroutine nulike_sensinit
