***********************************************************************
*** nulike_bgspec provides the interpolated probability distribution
*** function for the energy estimator due to background events.
***
*** Input:	ee	value of the energy estimator (e.g. nchan)
*** Output:             pdf (units of estimator^-1)
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Apr 24, 2011
*** Modified: Jun 6, 2014
***********************************************************************

      real*8 function nulike_bgspec(ee,like)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      real*8  ee, ee_a(1), nulike_bgspec_a(1)
      integer nchan_index, like, IER

      ! Switch according to likelihood version.
      select case (like)

      ! 2012 likelihood, as per arXiv:1207.0810
      case (2012)

        if (ee .lt. ee_min(analysis) .or. ee .gt. ee_max(analysis)) then
          write(*,*) 'Error in nulike_bgspec: energy estimator outside'
          write(*,*) 'tabulated range: ee=',ee,'.  Quitting...'
          stop
        endif

        nchan_index = nint(ee) - nint(ee_min(analysis)) + 1 - nchan_hist2BGoffset(analysis)
        if (BGeedist_ee(nchan_index,analysis) .ne. ee) then
          write(*,*) BGeedist_ee(nchan_index,analysis), ee, nchan_index
          stop 'Something is wrong with nchan_index in nulike_bgspec.'
        endif

        nulike_bgspec = BGeedist_prob(nchan_index,analysis)

      !2015 likelihood, as per arXiv:1601.00653
      case (2015)

        !If the value of the energy estimator is outside the range of this histogram, return zero.
        if (ee .lt. BGeedist_ee(1,analysis) .or. ee .gt. BGeedist_ee(nBinsBGE(analysis),analysis) ) then
          nulike_bgspec = 0.d0
        else !If the measured value of the energy estimator is in the range of this histogram, set the prob by interpolating.
          ee_a(1) = ee
          call TSVAL1(nBinsBGE(analysis),BGeedist_ee(:,analysis),
     &     BGeedist_prob(:,analysis),BGeedist_derivs(:,analysis),
     &     BGeedist_sigma(:,analysis),0,1,ee_a,nulike_bgspec_a,IER)
          if (IER .lt. 0) then
            write(*,*) 'TSVAL1 error from background spectral'
            write(*,*) 'distribution in nulike_bgspec, code:', IER
            stop
          endif
          nulike_bgspec = nulike_bgspec_a(1)
          if (nulike_bgspec .lt. 0.d0) then
            nulike_bgspec = 0.d0
          endif
        endif

      case default
        write(*,*) "Unrecognised likelihood version in nulike_bgspec."
        write(*,*) "Quitting..."
        stop

      end select

      end function nulike_bgspec
