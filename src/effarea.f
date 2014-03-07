***********************************************************************
*** nulike_effarea provides the interpolated effective area of the 
*** IceCube experiment to either neutrinos or anti-neutrinos.
***
*** Input:	log10E	        log(neutrino energy/GeV)
*** 		ptype	= 1	neutrinos
***			= 2	anti-neutrinos
*** Output:                     effective area (m^2)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
***********************************************************************

      real*8 function nulike_effarea(log10E, ptype)

      implicit none
      include 'nulike.h'

      real*8 log10E
      integer ptype, IER
      
      !Choose relevant species
      if (ptype .eq. 1) then

        !neutrinos
        call TSVAL1(nBinsEA,effArea_logEcentres,effArea_nu,
     &   effArea_nuderivs,effArea_nusigma,0,1,log10E,nulike_effarea,IER)

      else if (ptype .eq. 2) then

        !anti-neutrinos
        call TSVAL1(nBinsEA,effArea_logEcentres,effArea_nubar,
     &   effArea_nubarderivs,effArea_nubarsigma,0,1,
     &   log10E,nulike_effarea,IER)

      else

        write(*,*) 'Error in nulike_effarea:'
        write(*,*) 'unrecognised ptype; quitting...'

      endif

      if (IER .lt. 0) then
        write(*,*) 'TSVAL1 error from effective area'
        write(*,*) 'in nulike_effarea, code:', IER
        stop
      endif

      end function nulike_effarea
