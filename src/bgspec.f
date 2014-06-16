***********************************************************************
*** nulike_bgspec provides the interpolated probability distribution
*** function for the energy estimator due to background events. 
***
*** Input:	ee	value of the energy estimator (e.g. nchan)
*** Output:             pdf (units of estimator^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
*** Modified: Jun 6, 2014
***********************************************************************

      real*8 function nulike_bgspec(ee,like)

      implicit none
      include 'nulike.h'

      real*8  ee
      integer nchan_index, like

      if (ee .lt. ee_min(analysis) .or. ee .gt. ee_max(analysis)) then
        write(*,*) 'Error in nulike_bgspec: energy estimator outside'
        write(*,*) 'tabulated range: ee=',ee,'.  Quitting...'
        stop
      endif
  
      ! Switch according to likelihood version.
      select case (like)

      ! 2012 likelihood, as per arXiv:1207.0810
      case (2012)

        nchan_index = nint(ee) - nint(ee_min(analysis)) + 1 - nchan_hist2BGoffset(analysis)
        if (BGeedist_ee(nchan_index,analysis) .ne. ee) then
          write(*,*) BGeedist_ee(nchan_index,analysis), ee, nchan_index
          stop 'Something is wrong with nchan_index in nulike_bgspec.'
        endif

        nulike_bgspec = BGeedist_prob(nchan_index,analysis)

      !2014 likelihood, as per arXiv:141x.xxxx
      case (2014)

        !TSVAL1( 

        !nulike_bgspec, IER)

      case default
        write(*,*) "Unrecognised likelihood version in nulike_bgspec."
        write(*,*) "Quitting..."
        stop

      end select



      end function nulike_bgspec
