***********************************************************************
*** nulike_bgspec provides the interpolated probability distribution
*** function for the number of hit DOMs due to background events. 
***
*** Input:	nchan		number of hit DOMs
*** Output:                     pdf (chan^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
*** Modified: Jun 6, 2014
***********************************************************************

      real*8 function nulike_bgspec(nchan,like)

      implicit none
      include 'nulike.h'

      real*8  nchan
      integer nchan_index, like

      if (nchan .lt. nchan_min(analysis) .or. nchan .gt. nchan_max(analysis)) then
        write(*,*) 'Error in nulike_bgspec: nchan outside tabulated'
        write(*,*) 'range, nchan=',nchan,'.  Quitting...'
        stop
      endif
  
      ! Switch according to likelihood version.
      select case (like)

      ! 2012 likelihood, as per arXiv:1207.0810
      case (2012)

        nchan_index = nint(nchan) - nint(nchan_min(analysis)) + 1 - nchan_hist2BGoffset(analysis)
        if (BGnchandist_nchan(nchan_index,analysis) .ne. nchan) then
          write(*,*) BGnchandist_nchan(nchan_index,analysis), nchan, nchan_index
          stop 'Something is wrong with nchan_index in nulike_bgspec.'
        endif

        nulike_bgspec = BGnchandist_prob(nchan_index,analysis)

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
