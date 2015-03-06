***********************************************************************
*** nulike_dPhi_SdE provides the differential spectrum of neutrinos or
*** anti-neutrinos to be expected at IceCube, after weighting 
*** by the effective area and the angular loss factor (accounting for
*** the loss of neutrinos to regions beyond the analysis cone due to the
*** PSF).  Note that the annihilation rate in the Sun is not included.
*** The energy dispersion is also not taken into
*** account here, as this function is intended for use as either
*** a) input into an integral over all energies, which will give the
***    total number of epected events (essentially identical regardless
***    of whether energy dispersion is included or not)
*** b) a prefactor, to be weighted by the energy dispersion
*** This routine is used only with the 2012 likelihood.
*** 
*** Input:	log10E		log10(neutrino energy/GeV)
*** 		ptype	= 1	neutrinos
***			= 2	anti-neutrinos 
***             nuyield         external double function that returns
***                             the differential neutrino flux
***                             at the detector in units of m^-2 GeV^-1
***                             annihilation^-1
***             context         A c_ptr passed in to nuyield when it is called
***
*** Output:                     differential flux (GeV^-1 annihilation^-1)
***       
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Apr 24, 2011
*** Modified: March 6, 2014
***********************************************************************


      real*8 function nulike_dPhi_SdE(log10E,ptype,nuyield,context)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      real*8 log10E, effArea, angLossFac
      real*8 sigma, spec, nulike_sens, nulike_angres
      integer ptype
      type(c_ptr) context

      interface
        real(c_double) function nuyield(log10E,ptype,context) bind(c)
          use iso_c_binding, only: c_ptr, c_double, c_int
          implicit none
          real(c_double), intent(in) :: log10E
          integer(c_int), intent(in) :: ptype
          type(c_ptr), intent(inout) :: context
        end function
      end interface

      !Obtain differential neutrino or anti-neutrino 
      !flux spectrum as it arrives at the detector; spec in m^-2 GeV^-1 annihilation^-1
      spec = nuyield(log10E,ptype,context)

      !Obtain effective area for relevant species and energy; effArea in m^2
      effArea = nulike_sens(log10E, ptype)

      !Obtain angular resolution for nu and nubar with energy E; sigma in degrees
      sigma = nulike_angres(log10E)

      !Obtain angular loss factor for neutrinos and anti-neutrinos
      !with energy E, to account for events that will leak out of
      !the analysis cone.
      angLossFac = 1.d0-exp(phi_max_deg(analysis)*phi_max_deg(analysis)
     & /(-2.d0*sigma*sigma))

      !Put everything together to obtain the spectrum observed
      !by the neutrino telescope; nulike_dPhi_SdE in GeV^-1 annihilation^-1
      nulike_dPhi_SdE = spec * effArea * angLossFac

      end function nulike_dPhi_SdE


