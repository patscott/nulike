***********************************************************************
*** nulike_angres calculates the 1 sigma angular resolution of the IceCube
*** detector for neutrinos and anti-neutrinos with energy E, in degrees.
*** Here 1 sigma refers to a single dimension of a 2D Gaussian PSF;
*** in terms of absolute containment angle, this therefore corresponds to
*** 39.3% containment, not 68%.
***
*** Input:	log10E	        log(neutrino energy/GeV)
*** Output:			39.3% containment ang resolution (degrees)
***       
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Apr 24, 2011
*** Modified: Jun 6, 2014
***********************************************************************

      real*8 function nulike_angres(log10E)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      real*8 log10E, nulike_angres_a(1)
      integer IER
      
      call TSVAL1(nSensBins(analysis),sens_logEcentres(:,analysis),
     & sens_AngRes(:,analysis),sens_AngResderivs(:,analysis),
     & sens_AngRessigma(:,analysis),0,1,log10E,nulike_angres_a,IER)

      if (IER .lt. 0) then
        write(*,*) 'TSVAL1 error from angular'
        write(*,*) 'resolution in nulike_angres, code:', IER
        stop
      endif

      nulike_angres = nulike_angres_a(1)

      end function nulike_angres
