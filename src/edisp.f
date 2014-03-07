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
***********************************************************************

      real*8 function nulike_edisp(log10E, nchan)

      implicit none
      include 'nulike.h'

      real*8 log10E
      integer nchan, IER

      if (nchan .ne. nchansaved) then
        call nulike_edispcheckout(nchan)
      endif

      call TSVAL1(nHistograms,hist_logEcentres,edisp_prob,
     & edisp_derivs,edisp_sigma,0,1,log10E,nulike_edisp,IER)

      if (nulike_edisp .lt. 0.d0) nulike_edisp = 0.d0

      if (IER .lt. 0) then
        write(*,*) 'TSVAL1 error from energy dispersion'
        write(*,*) 'in nulike_edisp, code:', IER
        stop
      endif

      end function nulike_edisp
