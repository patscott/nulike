***********************************************************************
*** nulike_bias provides the interpolated analysis bias factors of the 
*** IceCube experiment to either neutrinos or anti-neutrinos.
***
*** This routine is used only with the 2015 likelihood.
***
*** Input:  log10E  log(neutrino energy/GeV)
***         ptype   = 1 neutrinos
***                 = 2 anti-neutrinos
*** Output: bias factor (dimensionless)
***       
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Jul 27, 2015
***********************************************************************

      real*8 function nulike_bias(log10E, ptype)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      real*8 log10E, log10E_a(1), nulike_bias_a(1)
      integer ptype, IER

      !Abort if no bias factors are actually being considered
      if (no_bias(analysis)) then
        nulike_bias = 1.d0
        return
      endif
 
      !Abort if outside the valid energy range.
      if (log10E .lt. bias_logE(1,1,analysis) .or.
     &    log10E .gt. bias_logE(2,nBiasBins(analysis),analysis) ) then
        nulike_bias = 0.d0
        return
      endif

      log10E_a(1) = log10E

      !Choose relevant species
      if (ptype .eq. 1) then

        !neutrinos
        call TSVAL1(nBiasBins(analysis),bias_logEcentres(:,analysis),
     &   bias_nu(:,analysis),bias_nuderivs(:,analysis),
     &   bias_nusigma(:,analysis),0,1,log10E_a,nulike_bias_a,IER)
        nulike_bias = nulike_bias_a(1)

      else if (ptype .eq. 2) then

        !anti-neutrinos
        call TSVAL1(nBiasBins(analysis),bias_logEcentres(:,analysis),
     &   bias_nubar(:,analysis),bias_nubarderivs(:,analysis),
     &   bias_nubarsigma,0,1,log10E_a,nulike_bias_a,IER)
        nulike_bias = nulike_bias_a(1)

      else

        write(*,*) 'Error in nulike_bias:'
        write(*,*) 'unrecognised ptype; quitting...'
        stop

      endif

      if (IER .lt. 0) then
        write(*,*) 'TSVAL1 error from bias calculation'
        write(*,*) 'in nulike_bias, code:', IER
        stop
      endif

      end function nulike_bias
