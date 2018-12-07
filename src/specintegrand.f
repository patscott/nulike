***********************************************************************
*** nulike_specintegrand provides the integrand for computing the
*** signal component of the spectral likelihood in nulike_speclike.
*** This routine is used only with the 2012 likelihood.
***
*** Input:        NumFun         = 1
***               X(1)           log10(neutrino energy/GeV)
*** Hidden input: nuyield_ptr%f  external double function that returns
***                              the differential muon/neutrino flux
***                              at the detector in units of m^-2
***                              GeV^-1 annihilation^-1
***               context_shared A c_ptr passed in to nuyield when it is called
***               nchanshare     number of hit DOMs for this event
***               thetashare     total number of predicted signal
***                              events (nu + nubar)
*** Output:       integrand     (chan^-1)
***
*** Note that the factor of ln(10) in the logarithmic integral has
*** been left for post-multiplication in order to increase efficiency.
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Apr 2011
*** Modified: March 6, 2014
***********************************************************************

      function nulike_specintegrand(NumFun,X) result(Value)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      integer, intent(in) :: NumFun
      real*8, intent(in) :: X(:)
      real*8 :: Value(NumFun)
      real*8 nulike_dP_SdE, edisp, specpdf, nulike_edisp

      !Return energy dispersion
      edisp = nulike_edisp(X(1),nchanshare,likelihood_version(analysis))
      !Return spectral probability distribution function
      specpdf = nulike_dP_SdE(X(1),thetashare,nuyield_ptr%f,context_shared)
      !Put them together, weight by E to give full integrand
      Value(1) = specpdf * edisp * 10.d0**X(1)

      end function nulike_specintegrand
