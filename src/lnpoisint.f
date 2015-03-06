***********************************************************************
*** Wrapper function to perform the integral required for producing a      
*** Poissonian likelihood marginalised over a systematic
*** uncertainty in the sensitivity of the experiment.  See
*** hep-ex/0202013 or arXiv:0909.3300 for details.
***
*** input:  counts_tot   total predicted number of counts
***	    counts_sig   predicted number of counts due to signal
***         obscounts    observed number of counts
***         errorsq      squared fractional systematic error on 
***                       experimental sensitivity
***
*** output: nulike_lnpoisint ln(integral)
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: June, July 2011
***********************************************************************

      real*8 function nulike_lnpoisint(counts_tot,counts_sig,obscounts,
     & errorsq,logNormal)

      implicit none

      integer obscounts
      real*8 counts_tot, counts_sig, errorsq
      real*8 nulike_lnpiln, nulike_lnpin
      logical logNormal

      if (logNormal) then
        !Treat percentage error as log-normal distributed
        nulike_lnpoisint = nulike_lnpiln(obscounts,counts_tot-counts_sig,
     &   counts_sig,dsqrt(errorsq))
      else
        !Treat percentage error as Gaussianly-distributed
        nulike_lnpoisint = nulike_lnpin(obscounts,counts_tot-counts_sig,
     &   counts_sig,dsqrt(errorsq))
      endif


      end function nulike_lnpoisint

