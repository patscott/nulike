***********************************************************************
*** nulike_edisp provides the interpolated IceCube energy dispersion
*** estimator for a given incoming neutrino energy and observed
*** number of hit DOMs.
***
*** Input:	log10E	        log(neutrino energy/GeV)
*** 		ee		value of the energy estimator (e.g. nchan)
*** Output:                     energy dispersion (units of ee^-1)
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Apr 24, 2011
*** Modified: Jun 7, 15 2014
***********************************************************************

      real*8 function nulike_edisp(log10E, ee, like)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'
      include 'nuprep.h'

      real*8 log10E, ee, nulike_edisp_a(1)
      integer nchan_index, IER, like

      !Switch according to likelihood version.
      select case (like)

      !2012 likelihood, as per arXiv:1207.0810
      case (2012)

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
     &   hist_prob(:,nchan_index,analysis),hist_derivs(:,nchan_index,analysis),
     &   hist_sigma(:,nchan_index,analysis),0,1,log10E,nulike_edisp_a,IER)
        nulike_edisp = nulike_edisp_a(1)

        if (log10E .lt. hist_logEcentres(1,analysis)) then
          nulike_edisp = hist_prob(1,nchan_index,analysis)
        else if (log10E .gt. hist_logEcentres(nHistograms(analysis),analysis)) then
          nulike_edisp = hist_prob(nHistograms(analysis),nchan_index,analysis)
        endif

      !2015 likelihood, as per arXiv:1601.00653
      case (2015)

        call TSVAL1(nhist,hist_logEnergies,hist_single_ee_prob,
     &   hist_single_ee_derivs,hist_single_ee_sigma,0,1,log10E,nulike_edisp_a,IER)
        nulike_edisp = nulike_edisp_a(1)

      case default
        write(*,*) "Unrecognised likelihood version in nulike_init."
        write(*,*) "Quitting..."
        stop

      end select

      if (nulike_edisp .lt. 0.d0) nulike_edisp = 0.d0

      if (IER .lt. 0) then
        write(*,*) 'TSVAL1 error from energy dispersion'
        write(*,*) 'in nulike_edisp, code:', IER
        stop
      endif

      end function nulike_edisp
