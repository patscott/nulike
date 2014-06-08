***********************************************************************
*** nulike_tabulated_weight provides access to the the precomputed 
*** weights to be used in
*** a) the signal component of the combined angular-spectral likelihood
***    in nulike_specanglike
*** b) the total predicted number of signal events in nulike_signal.
*** This routine is used only with the 2014 likelihood.
***
*** Input:		log10E       log10(neutrino energy/GeV)
***                     ptype        1=nu, 2=nubar
*** 			eventnum     the unique index number of the event
***
*** Output:             weight       eventnum = 0 => effective area (m^2)
***                                  eventnum > 0 => weighting for event
***                                   eventnum (m^2 chan^-1 degrees^-1 )
***
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Jun 8, 2014
***********************************************************************


      real*8 function nulike_tabulated_weight(log10E,ptype,eventnum)

      implicit none
      include 'nulike.h'

      real*8 log10E
      integer ptype, eventnum

      nulike_tabulated_weight = 1.d0

      end function nulike_tabulated_weight
