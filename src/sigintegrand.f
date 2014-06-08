***********************************************************************
*** nulike_sigintegrad provides the integrand for computing theta_S,
*** the number of expected signal events at IceCube.
*** Input:		log10E		log10(neutrino energy/GeV)
***                     nuyield         external double function that returns
***                                     the differential neutrino flux
***                                     at the detector in units of m^-2 
***                                     GeV^-1 annihilation^-1
*** Hidden Input:	ptypeshare   1 => return integrand for neutrinos
***                                  2 => for anti-neutrinos 
*** Output:             integrand       dimensionless
***
*** Note that the factor of ln(10) in the logarithmic integral has
*** been left for post-multiplication in order to increase efficiency.      
***
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
***********************************************************************

      real*8 function nulike_sigintegrand(log10E,nuyield)

      implicit none
      include 'nulike.h'

      real*8 log10E, nulike_dPhi_SdE, nuyield
      external nuyield

      !Return either predicted differential neutrino or
      !anti-neutrino signal in IceCube, depending on 
      !ptypeshare
      nulike_sigintegrand = nulike_dPhi_SdE(log10E,ptypeshare,
     & nuyield)*10.d0**log10E

      end function nulike_sigintegrand
