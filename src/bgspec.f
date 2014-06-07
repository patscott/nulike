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

      real*8 function nulike_bgspec(nchan)

      implicit none
      include 'nulike.h'

      integer nchan, nchan_index

      if (nchan .lt. nchan_min(analysis) .or. nchan .gt. nchan_max(analysis)) then
        write(*,*) 'Error in nulike_bgspec: nchan outside tabulated'
        write(*,*) 'range, nchan=',nchan,'.  Quitting...'
        stop
      endif
  
      nchan_index = nchan - nchan_min(analysis) + 1 - nchan_hist2BGoffset(analysis)
      if (BGnchandist_nchan(nchan_index,analysis) .ne. nchan) then
        write(*,*) BGnchandist_nchan(nchan_index,analysis), nchan, nchan_index
        stop 'Something is wrong with nchan_index in nulike_bgspec.'
      endif

      nulike_bgspec = BGnchandist_prob(nchan_index,analysis)

      end function nulike_bgspec
