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
***********************************************************************

      subroutine nulike_bginit(filename, nbins_angular, nbins_nchan, 
     & first, second)

      implicit none
      include 'nulike.h'

      character (len=*) filename
      character (len=200) instring
      character (len=15) headerstring
      integer nbins_angular, nbins_nchan, counts(3), IER
      integer i, j, dummyint, first, second, indices(3)
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
              BGnchandist_nchan(j) = int(dummyfloat1)
              BGnchandist_prob(j) = dummyfloat2
            endif
          else
            FullSkyBG = dummyint
          endif
          read(lun, *) instring
       enddo
       if (i .ne. 3) read(lun, fmt=*), instring
      enddo 

      close(lun)

      !Throw a warning if energy dispersion files don't go low enough
      !in nchan to cover whole tabulated range in BG file, and then mark
      !them for extension.
      if (nchan_min .gt. BGnchandist_nchan(1)) then
       ! write(*,*) 'Warning from nulike_bginit: nchan
     & !values in the observed background spectrum go below
     & !the range tabulated in the energy dispersion histograms.
     & !Assuming zeros for histograms entries outside given range.'
        nchan_min = BGnchandist_nchan(1)
      endif

      !Throw a warning if energy dispersion files don't go high enough
      !in nchan to cover whole tabulated range in BG file, and then mark
      !them for extension.
      if (nchan_max .lt. BGnchandist_nchan(nBinsBGE)) then
       ! write(*,*) 'Warning from nulike_bginit: nchan
     & !values in the observed background spectrum go above
     & !the range tabulated in the energy dispersion histograms.
     & !Assuming zeros for histograms entries outside given range.'
        nchan_max = BGnchandist_nchan(nBinsBGE)
      endif

      !Reset nnchan_total
      nnchan_total = nchan_max - nchan_min + 1

      !Make sure we didn't break everything
      if (nnchan_total .gt. max_nnchan) then
        write(*,*)
        write(*,*) 'Extension of histograms gives more nchan values'
        write(*,*) 'than DarkSUSY has been configured to handle.  '
        write(*,*) 'Increase max_nnchan in nulike.h and recompile.'
        write(*,*)
        call exit(0)
      endif
      

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
        BGangdist_phi(i) = 
     &   BGangdist_phi_temp(nbins_angular+1-i)/180.d0 * pi
        !Convert probabilities from dP/dphi to dP/dcos(phi) and flip em
        BGangdist_prob(i) = BGangdist_prob_temp(nbins_angular+1-i)/
     &   dsin(BGangdist_phi(i)) * 180.d0 / pi
      enddo

      !Initialise interpolator
      call TSPSI(nbins_angular,dcos(BGangdist_phi),BGangdist_prob,
     & 2,0,.false.,.false.,2*nbins_angular-2,working,BGangdist_derivs,
     & BGangdist_sigma,IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_bgnit: TSPSI failed with error'
        write(*,*) 'code',IER
        stop
      endif

      !Calculate renormalisation factor required to cancel any tiny normalisation
      !change introduced by interpolation (typically order 1e-3)
      BGangdist_norm = TSINTL (-1.d0,1.d0,nbins_angular,
     & dcos(BGangdist_phi),BGangdist_prob,BGangdist_derivs,
     & BGangdist_sigma,IER)     
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_bgnit: TSINTL failed with error'
        write(*,*) 'code',IER
        stop
      endif

      !Make sure nchan histograms are properly normalised
      BGnchandist_prob = BGnchandist_prob/sum(BGnchandist_prob)


      end subroutine nulike_bginit

