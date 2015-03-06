***********************************************************************
*** nulike_dP_SdE provides the differential probability for any neutrino
*** or anti-neutrino signal event in IceCube to have the energy E.
*** This routine is used only with the 2012 likelihood.
***
*** Input:	log10E	  log10(neutrino energy/GeV)
***		tS_tot	  the total number of signal events expected,
***			  of any energy or type (nu or nubar)
***             nuyield   external double function that returns
***                       the differential neutrino flux
***                       at the detector in units of m^-2 GeV^-1
***                       annihilation^-1
***             context   A c_ptr passed in to nuyield when it is called
***
*** Output:               differential probability (GeV^-1)
***       
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Apr 24, 2011
*** Modified: March 6, 2014
***********************************************************************

      real*8 function nulike_dP_SdE(log10E,tS_tot,nuyield, context)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      real*8 log10E, tS_tot, nulike_dPhi_SdE, nuyield
      type(c_ptr) context
      external nuyield

      !Obtain the neutrino contribution in GeV^-1 annihilation^-1
      nulike_dP_SdE = nulike_dPhi_SdE(log10E,1,nuyield,context)
      !Obtain the anti-neutrino contribution in GeV^-1 annihilation^-1
      nulike_dP_SdE = nulike_dP_SdE + nulike_dPhi_SdE(log10E,2,nuyield,context)
      !Normalise the distribution to differential probability per GeV
      !final nulike_dP_SdE in GeV^-1, exp_time in s, annrate in s^-1
      nulike_dP_SdE = nulike_dP_SdE * exp_time(analysis) * annrateshare / tS_tot

      end function nulike_dP_SdE

