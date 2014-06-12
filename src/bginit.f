***********************************************************************
*** nulike_bginit initialises the pdfs which describe the IceCube
*** background spectrum and angular distribution, and reads in
*** the total number of background events.
***
*** input:   filename       name of file containing background 
***                          distributions
***          nbins_angular  number of bins for angular distribution
***          nbins_nchan    number of bins for nchan (energy) 
***                          distribution
***          first, second  indexes between 1 and 3 indicating
***                          the identities of the first and second
***                          blocks in the file (angular, nchan or 
***                          number of events).  See nulike.h for
***                          the key.
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: April 8, 2011
*** Modified: Jun 6, 2014
***********************************************************************

      subroutine nulike_bginit(filename, nbins_angular, nbins_nchan, 
     & first, second, like)

      implicit none
      include 'nulike.h'

      character (len=*) filename
      character (len=200) instring
      character (len=15) headerstring
      integer nbins_angular, nbins_nchan, counts(3), IER
      integer i, j, dummyint, first, second, indices(3), like
      real*8 dummyfloat1, dummyfloat2
      real*8 BGangdist_phi_temp(max_nBinsBGAng)
      real*8 BGangdist_prob_temp(max_nBinsBGAng)
      real*8 working(2*nbins_angular-2), TSINTL

      !Open background file for reading
      open(lun,file=filename,ACTION='READ')

      !Set record ordering parameters
      headerstring = hstring(first)
      indices = (/first, second, 6 - first - second/)
      do i = 1,3
        select case (indices(i))
          case(angular)
            counts(i) = nbins_angular
          case(nchannels)
            counts(i) = nbins_nchan
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
              !Read in observed distribution of nchan (energies)
              BGnchandist_nchan(j,analysis) = dummyfloat1
              BGnchandist_prob(j,analysis) = dummyfloat2
            endif
          else
            FullSkyBG(analysis) = dummyint
          endif
          read(lun, *) instring
       enddo
       if (i .ne. 3) read(lun, fmt=*), instring
      enddo 

      close(lun)

      ! Switch according to likelihood version.
      select case (like)

      ! 2012 likelihood, as per arXiv:1207.0810
      case (2012)
        !If the energy dispersion files don't go high or low enough in nchan
        !to cover the whole tabulated range in the BG file, mark them for extension.
        if (nchan_min(analysis) .gt. BGnchandist_nchan(1,analysis)) 
     &   nchan_min(analysis) = BGnchandist_nchan(1,analysis)
        if (nchan_max(analysis) .lt. BGnchandist_nchan(nbins_nchan,analysis))
     &   nchan_max(analysis) = BGnchandist_nchan(nbins_nchan,analysis)
        !Reset nnchan_total
        nnchan_total(analysis) = nint(nchan_max(analysis) - nchan_min(analysis)) + 1
        !Make sure we didn't break everything
        if (nnchan_total(analysis) .gt. max_nnchan) then
          write(*,*)
          write(*,*) 'Extension of histograms gives more nchan values'
          write(*,*) 'than nulike has been configured to handle.  '
          write(*,*) 'Increase max_nnchan in nulike.h and recompile.'
          write(*,*)
          call exit(0)
        endif

      !2014 likelihood, as per arXiv:141x.xxxx (Set up interpolation in distribution of the energy estimator.)
      case (2014)
        !FIXME

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
     & 2,0,.false.,.false.,2*nbins_angular-2,working,BGangdist_derivs(:,analysis),
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

      !Make sure nchan histograms are properly normalised
      BGnchandist_prob(:,analysis) = BGnchandist_prob(:,analysis)/
     &                               sum(BGnchandist_prob(:,analysis))


      end subroutine nulike_bginit

