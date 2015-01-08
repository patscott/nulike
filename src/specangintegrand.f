***********************************************************************
*** nulike_specangintegrand provides the integrand for computing
*** a) the signal component of the combined angular-spectral likelihood
***    in nulike_specanglike
*** b) the total predicted number of signal events in nulike_signal.
*** This routine is used only with the 2014 likelihood.
***
*** Input:		log10E         log_10(neutrino energy/GeV)
***                     nuyield        external double function that returns
***                                     the differential muon/neutrino flux
***                                     at the detector in units of m^-2 
***                                     GeV^-1 annihilation^-1
*** Hidden Input:	eventnumshare  the unique index number of this event
***
*** Output:             integrand      eventnum = 0 => (annihilation^-1)
***                                    eventnum > 0 => (annihilation^-1
***                                                     chan^-1 degrees^-1 )
***
*** Note that the factor of ln(10) in the logarithmic integral has
*** been left for post-multiplication in order to increase efficiency.      
***
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Jun 8, 2014
***********************************************************************


      real*8 function nulike_specangintegrand(log10E,nuyield,context)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      real*8 log10E, nuyield, nulike_tabulated_weight
      type(c_ptr) context
      external nuyield

      !Scale neutrino and anti-neutrino yields by saved weights.
      nulike_specangintegrand =  
     & nuyield(log10E,1,context)*nulike_tabulated_weight(log10E,1,eventnumshare) +
     & nuyield(log10E,2,context)*nulike_tabulated_weight(log10E,2,eventnumshare)

      !Weight by E to give full integrand
      nulike_specangintegrand = nulike_specangintegrand * 10.d0**log10E

      !write(*,*) 10.d0**log10E, nulike_specangintegrand

      end function nulike_specangintegrand
