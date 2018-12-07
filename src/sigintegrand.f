***********************************************************************
*** nulike_sigintegrad provides the integrand for computing theta_S,
*** the number of expected signal events at IceCube.
*** This routine is used only with the 2012 likelihood.
***
*** Input:		NumFun         =1
***               X(1)           log10(neutrino energy/GeV)
*** Hidden input: nuyield_ptr%f   external double function that returns
***                               the differential neutrino flux
***                               at the detector in units of m^-2
***                               GeV^-1 annihilation^-1
***               ptypeshare     1 => return integrand for neutrinos
***                              2 => for anti-neutrinos
***               context_shared A c_ptr passed in to nuyield when it is
***                               called
*** Output:       integrand      dimensionless
***
*** Note that the factor of ln(10) in the logarithmic integral has
*** been left for post-multiplication in order to increase efficiency.
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Apr 24, 2011
***********************************************************************

      function nulike_sigintegrand(NumFun,X) result(Value)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      integer, intent(in) :: NumFun
      real*8, intent(in) :: X(:)
      real*8 :: Value(NumFun)
      real*8 nulike_dPhi_SdE

      !Return either predicted differential neutrino or
      !anti-neutrino signal in IceCube, depending on
      !ptypeshare
      Value(1) = nulike_dPhi_SdE(X(1),ptypeshare,
     & nuyield_ptr%f,context_shared)*10.d0**X(1)

      end function nulike_sigintegrand
