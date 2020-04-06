***********************************************************************
*** nulike_bginit initialises the pdfs which describe the IceCube
*** background spectrum and angular distribution, and reads in
*** the total number of background events.
***
*** input:   filename       name of file containing background
***                          distributions
***          nbins_angular  number of bins for angular distribution
***          nbins_ee       number of bins for energy estimator distribution
***          first, second  indexes between 1 and 3 indicating
***                          the identities of the first and second
***                          blocks in the file (angular, energy
***                          estimator or number of events).
***                          See nucommon.h for the key.
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: April 8, 2011
*** Modified: Jun 6, 2014
***********************************************************************

      subroutine nulike_bginit(filename, nbins_angular, nbins_ee,
     & first, second, like)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      character (len=*) filename
      character (len=200) instring
      character (len=15) headerstring
      integer nbins_angular, nbins_ee, counts(3), IER
      integer i, j, dummyint, first, second, indices(3), like
      real*8 dummyfloat1, dummyfloat2
      real*8 BGangdist_phi_temp(max_nBinsBGAng)
      real*8 BGangdist_prob_temp(max_nBinsBGAng)
      real*8 angworking(2*nbins_angular-2)
      real*8 eeworking(2*nbins_ee-2), TSINTL


      !Open background file for reading
      open(lun,file=filename,ACTION='READ')

      !Set record ordering parameters
      headerstring = hstring(first)
      indices = (/first, second, 6 - first - second/)
      do i = 1,3
        select case (indices(i))
          case(angular)
            counts(i) = nbins_angular
          case(enrgyest)
            counts(i) = nbins_ee
          case(events)
            counts(i) = 1
        end select
      enddo

      !Skip over header lines
      instring = '#'
      do while (instring .ne. headerstring)
        read(lun, *) instring
      enddo

      !Read in background angular and energy (estimator) distributions, as well
      !as total number of BG events across the whole sky
      do i = 1, 3
        do j = 1, counts(i)
          read(lun, *) instring, dummyint
          if (indices(i).ne.events) then
            read(lun, *) instring, dummyfloat1, dummyfloat2
            if (indices(i).eq.angular) then
              !Read in observed angular distribution of background events
              BGangdist_phi_temp(j) = dummyfloat1
              BGangdist_prob_temp(j) = dummyfloat2
            else
              !Read in observed distribution of energy estimators (e.g. nchan)
              BGeedist_ee(j,analysis) = dummyfloat1
              BGeedist_prob(j,analysis) = dummyfloat2
            endif
          else
            FullSkyBG(analysis) = dummyint
          endif
          read(lun, *) instring
       enddo
       if (i .ne. 3) read(lun, fmt=*) instring
      enddo

      close(lun)

      ! Switch according to likelihood version.
      select case (like)

      ! 2012 likelihood, as per arXiv:1207.0810
      case (2012)
        !If the energy dispersion files don't go high or low enough in nchan
        !to cover the whole tabulated range in the BG file, mark them for extension.
        if (ee_min(analysis) .gt. BGeedist_ee(1,analysis))
     &   ee_min(analysis) = BGeedist_ee(1,analysis)
        if (ee_max(analysis) .lt. BGeedist_ee(nbins_ee,analysis))
     &   ee_max(analysis) = BGeedist_ee(nbins_ee,analysis)
        !Reset nnchan_total
        nnchan_total(analysis) = nint(ee_max(analysis) - ee_min(analysis)) + 1
        !Make sure we didn't break everything
        if (nnchan_total(analysis) .gt. max_ncols) then
          write(*,*)
          write(*,*) 'Extension of histograms gives more nchan values'
          write(*,*) 'than nulike has been configured to handle.  '
          write(*,*) 'Increase max_ncols in nucommon.h and recompile.'
          write(*,*)
          call exit(0)
        endif

      !2015 likelihood, as per arXiv:1601.00653 (Set up interpolation in distribution of the energy estimator.)
      case (2015)

        !Set up interpolation in energy distribution.
        call TSPSI(nbins_ee,BGeedist_ee(:,analysis),BGeedist_prob(:,analysis),2,0,.false.,.false.,
     &   2*nbins_ee-2,eeworking,BGeedist_derivs(:,analysis),BGeedist_sigma(:,analysis),IER)
        if (IER .lt. 0) then
          write(*,*) 'TSPSI error from spectral distribution of'
          write(*,*) 'background in nulike_bginit, code: ', IER
          stop
        endif

      case default
        write(*,*) "Unrecognised likelihood version in nulike_bginit."
        write(*,*) "Quitting..."
        stop

      end select


      !Set up interpolation in angular distribution

      do i = 1, nbins_angular - 1
        !Take bin centres for angular values
        BGangdist_phi_temp(i) = 0.5d0 * (BGangdist_phi_temp(i)
     &                                + BGangdist_phi_temp(i+1))
      enddo
      !Do the same for the last bin, assuming its upper limit is 180 degrees
      BGangdist_phi_temp(nbins_angular) =
     & 0.5d0 * (BGangdist_phi_temp(nbins_angular) + 180.d0)

      do i = 1, nbins_angular
        !Convert angles to radians and flip em (for fussy interpolator)
        BGangdist_phi(i,analysis) =
     &   BGangdist_phi_temp(nbins_angular+1-i)/180.d0 * pi
        !Convert probabilities from dP/dphi to dP/dcos(phi) and flip em
        BGangdist_prob(i,analysis) = BGangdist_prob_temp(nbins_angular+1-i)/
     &   dsin(BGangdist_phi(i,analysis)) * 180.d0 / pi
      enddo
      BGangdist_phi(:,analysis) = dcos(BGangdist_phi(:,analysis))

      !Initialise interpolator
      call TSPSI(nbins_angular,BGangdist_phi(:,analysis),BGangdist_prob(:,analysis),
     & 2,0,.false.,.false.,2*nbins_angular-2,angworking,BGangdist_derivs(:,analysis),
     & BGangdist_sigma(:,analysis),IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_bgnit: TSPSI failed with error'
        write(*,*) 'code',IER
        stop
      endif

      !Calculate renormalisation factor required to cancel any tiny normalisation
      !change introduced by interpolation (typically order 1e-3)
      BGangdist_norm(analysis) = TSINTL (-1.d0,1.d0,nbins_angular,
     & BGangdist_phi(:,analysis),BGangdist_prob(:,analysis),
     & BGangdist_derivs(:,analysis),BGangdist_sigma(:,analysis),IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_bgnit: TSINTL failed with error'
        write(*,*) 'code',IER
        stop
      endif

      !Make sure energy estimator histograms are properly normalised
      BGeedist_prob(:,analysis) = BGeedist_prob(:,analysis)/sum(BGeedist_prob(:,analysis))


      end subroutine nulike_bginit

