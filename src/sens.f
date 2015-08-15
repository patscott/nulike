***********************************************************************
*** nulike_sens provides the interpolated effective area/volume of a  
*** loaded experiment/analysis to either neutrinos or anti-neutrinos.
***
*** Input:  log10E  log(neutrino energy/GeV)
***         ptype   = 1 neutrino efefctive area
***                 = 2 anti-neutrinos effective area
***                 = 3 CP-invariant effective volume
*** Output: effective area (m^2) or volume (km^3)
***       
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Apr 24, 2011
*** Modified: Jun 6 2014
***********************************************************************

      real*8 function nulike_sens(log10E, ptype)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      real*8 log10E, log10E_a(1), nulike_sens_a(1)
      integer ptype, IER
 
      !Abort if outside the valid energy range.     
      if (log10E .lt. sens_logE(1,1,analysis) .or.
     &    log10E .gt. sens_logE(2,nSensBins(analysis),analysis) ) then
        nulike_sens = 0.d0
        return
      endif

      log10E_a(1) = log10E

      !Choose relevant species/quantity
      if (ptype .eq. 1 .or. ptype .eq. 3) then

        !neutrino effective area or CP-invariant effective volume 
        call TSVAL1(nSensBins(analysis),sens_logEcentres(:,analysis),
     &   sens_nu(:,analysis),sens_nuderivs(:,analysis),
     &   sens_nusigma(:,analysis),0,1,log10E_a,nulike_sens_a,IER)
        nulike_sens = nulike_sens_a(1)

      else if (ptype .eq. 2) then

        !anti-neutrinos
        call TSVAL1(nSensBins(analysis),sens_logEcentres(:,analysis),
     &   sens_nubar(:,analysis),sens_nubarderivs(:,analysis),
     &   sens_nubarsigma,0,1,log10E_a,nulike_sens_a,IER)
        nulike_sens = nulike_sens_a(1)

      else

        write(*,*) 'Error in nulike_sens:'
        write(*,*) 'unrecognised ptype; quitting...'
        stop

      endif

      if (IER .lt. 0) then
        write(*,*) 'TSVAL1 error from effective area'
        write(*,*) 'in nulike_sens, code:', IER
        stop
      endif

      end function nulike_sens
