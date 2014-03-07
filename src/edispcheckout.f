***********************************************************************
*** nulike_edispcheckout sets up the correct interpolation vectors for
*** nulike_edisp, corresponding to the requested number of channels.
***
*** Input:	nchan		number of hit DOMs
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
***********************************************************************

      subroutine nulike_edispcheckout(nchan)

      implicit none
      include 'nulike.h'

      integer nchan, nchan_index, j

      if (nchan .lt. nchan_min .or. nchan .gt. nchan_max) then
        write(*,*) 'Error in nulike_edispcheckout: nchan outside'
        write(*,*) 'tabulated range, nchan=',nchan,'.  Quitting...'
        stop
      endif

      nchansaved = nchan
      nchan_index = nchan - nchan_min + 1
      if (hist_nchan(1,nchan_index) .ne. nchan) then
        stop'Something is wrong with nchan_index in nulike_edispcheckout'
      endif

      do j = 1, nHistograms
        edisp_prob(j) = hist_prob(j,nchan_index)
        edisp_derivs(j) = hist_derivs(j,nchan_index)
        edisp_sigma(j) = hist_sigma(j,nchan_index)
      enddo

      end subroutine nulike_edispcheckout

