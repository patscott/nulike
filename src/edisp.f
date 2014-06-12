***********************************************************************
*** nulike_edisp provides the interpolated IceCube energy dispersion
*** estimator for a given incoming neutrino energy and observed
*** number of hit DOMs. 
***
*** Input:	log10E	        log(neutrino energy/GeV)
*** 		nchan		number of hit DOMs
*** Output:                     energy dispersion (chan^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
*** Modified: Jun 7, 2014
***********************************************************************

      real*8 function nulike_edisp(log10E, nchan)

      implicit none
      include 'nulike.h'

      real*8 log10E, nchan
      integer nchan_index, IER

      if (nchan .lt. nchan_min(analysis) .or. nchan .gt. nchan_max(analysis)) then
        write(*,*) 'Error in nulike_edisp: nchan outside'
        write(*,*) 'tabulated range, nchan=',nchan,'.  Quitting...'
        stop
      endif

      nchan_index = nint(nchan) - nint(nchan_min(analysis)) + 1
      if (hist_nchan(1,nchan_index,analysis) .ne. nchan) then
        stop'Something is wrong with nchan_index in nulike_edisp'
      endif

      call TSVAL1(nHistograms(analysis),hist_logEcentres(:,analysis),
     & hist_prob(:,nchan_index,analysis),hist_derivs(:,nchan_index,analysis),
     & hist_sigma(:,nchan_index,analysis),0,1,log10E,nulike_edisp,IER)

      if (nulike_edisp .lt. 0.d0) nulike_edisp = 0.d0

      if (IER .lt. 0) then
        write(*,*) 'TSVAL1 error from energy dispersion'
        write(*,*) 'in nulike_edisp, code:', IER
        stop
      endif

      end function nulike_edisp
