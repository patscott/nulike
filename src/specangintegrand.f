***********************************************************************
*** nulike_specangintegrand provides the integrand for computing
*** a) the signal component of the combined angular-spectral likelihood
***    in nulike_specanglike
*** b) the total predicted number of signal events in nulike_signal.
*** This routine is used only with the 2015 likelihood.
***
*** Input:        NumFun         =1
***               X(1)           log10(neutrino energy/GeV)
*** Hidden input: nuyield_ptr%f   external double function that returns
***                               the differential neutrino flux
***                               at the detector in units of m^-2
***                               GeV^-1 annihilation^-1
***               eventnumshare  the unique index number of this event
***                               (indexed by the thread number)
***               context_shared A c_ptr passed in to nuyield when it is
***                               called
***
*** Output:       integrand      eventnum = 0 => (annihilation^-1)
***                              eventnum > 0 => (annihilation^-1
***                                               chan^-1 degrees^-1 )
***
*** Note that the factor of ln(10) in the logarithmic integral has
*** been left for post-multiplication in order to increase efficiency.
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Jun 8, 2014
***********************************************************************

      function nulike_specangintegrand(NumFun,X) result(Value)

      use iso_c_binding, only: c_ptr
      use omp_lib

      implicit none
      include 'nulike_internal.h'

      integer, intent(in) :: NumFun
      integer :: thread
      real*8, intent(in) :: X(:)
      real*8 :: Value(NumFun)
      real*8 nulike_tabulated_weight, nulike_bias, w1, w2

      !Get weights
      thread = omp_get_thread_num() + 1
      w1 = nulike_tabulated_weight(X(1),1,eventnumshare(thread))
      w2 = nulike_tabulated_weight(X(1),2,eventnumshare(thread))

      !Scale neutrino and anti-neutrino yields by saved weights.
      Value(1) = nuyield_ptr%f(X(1),1,context_shared)*nulike_bias(X(1),1)*w1 +
     &           nuyield_ptr%f(X(1),2,context_shared)*nulike_bias(X(1),2)*w2

      !Weight by E to give full integrand
      Value(1) = Value(1) * 10.d0**X(1)

      end function nulike_specangintegrand
