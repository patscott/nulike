***********************************************************************
*** nulike_specintegrand provides the integrand for computing the
*** signal component of the spectral likelihood in nulike_speclike.
***
*** Input:		log10E       log10(neutrino energy/GeV)
***                     muonyield    external double function that returns
***                                  the differential muon/neutrino flux
***                                  at the detector in units of m^-2 GeV^-1
*** Hidden Input:	nchanshare   number of hit DOMs for this event
***                     thetashare   total number of predicted signal 
***                                  events (nu + nubar) 
*** Output:             integrand    (chan^-1)
***
*** Note that the factor of ln(10) in the logarithmic integral has
*** been left for post-multiplication in order to increase efficiency.      
***
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 2011
*** Modified: March 6, 2014
***********************************************************************

      real*8 function nulike_specintegrand(log10E,muonyield)

      implicit none
      include 'nulike.h'

      real*8 log10E, nulike_dP_SdE, edisp, specpdf
      real*8 nulike_edisp, muonyield
      external muonyield

      !Return energy dispersion
      edisp = nulike_edisp(log10E,nchanshare)
      !Return spectral probability distribution function
      specpdf = nulike_dP_SdE(log10E,thetashare,muonyield)
      !Put them together, weight by E to give full integrand
      nulike_specintegrand = specpdf * edisp * 10.d0**log10E

      end function nulike_specintegrand
