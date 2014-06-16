***********************************************************************
*** nulike_edisp provides the interpolated IceCube energy dispersion
*** estimator for a given incoming neutrino energy and observed
*** number of hit DOMs. 
***
*** Input:	log10E	        log(neutrino energy/GeV)
*** 		ee		value of the energy estimator (e.g. nchan)
*** Output:                     energy dispersion (units of ee^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
*** Modified: Jun 7, 15 2014
***********************************************************************

      real*8 function nulike_edisp(log10E, ee)

      implicit none
      include 'nulike.h'

      real*8 log10E, ee
      integer nchan_index, IER

      if (ee .lt. ee_min(analysis) .or. ee .gt. ee_max(analysis)) then
        write(*,*) 'Error in nulike_edisp: energy estimator outside'
        write(*,*) 'tabulated range: ee=',ee,'.  Quitting...'
        stop
      endif

      nchan_index = nint(ee) - nint(ee_min(analysis)) + 1
      if (hist_nchan(1,nchan_index,analysis) .ne. nint(ee)) then
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
