***********************************************************************
*** nulike_dP_SdE provides the differential probability for any neutrino
*** or anti-neutrino signal event in IceCube to have the energy E.
***
*** Input:	log10E	  log10(neutrino energy/GeV)
***		tS_tot	  the total number of signal events expected,
***			  of any energy or type (nu or nubar)
***             muonyield external double function that returns
***                       the differential muon/neutrino flux
***                       at the detector in units of m^-2 GeV^-1
***
*** Output:               differential probability (GeV^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
*** Modified: March 6, 2014
***********************************************************************

      real*8 function nulike_dP_SdE(log10E,tS_tot,muonyield)

      implicit none
      include 'nulike.h'

      real*8 log10E, tS_tot, nulike_dPhi_SdE, muonyield
      external muonyield

      !Obtain the neutrino contribution in yr^-1 GeV^-1
      nulike_dP_SdE = nulike_dPhi_SdE(log10E,1,muonyield)
      !Obtain the anti-neutrino contribution in yr^-1 GeV^-1
      nulike_dP_SdE = nulike_dP_SdE + nulike_dPhi_SdE(log10E,2,muonyield)
      !Normalise the distribution to differential probability per GeV
      !nulike_dP_SdE in GeV^-1, exp_time in s, annrate in s^-1
      nulike_dP_SdE = nulike_dP_SdE * exp_time * annrate / tS_tot

      end function nulike_dP_SdE

