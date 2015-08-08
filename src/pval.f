***********************************************************************
*** nulike_pval performs the calculation of the p-value 
*** in IceCube calculations based on Poissonian statistics.
*** 
*** input:  ntot      observed number of events
***         theta_tot total predicted number of events
***         theta_sig predicted number of signal events
***         err_th    fractional theoretical error
*** output:           p value
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: May, June, July, Dec 2011
*** Modified: Jun 3, 6 2014
***********************************************************************

      double precision function nulike_pval(ntot,theta_tot,theta_sig,err_th)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      integer ntot
      real*8 theta_tot, theta_sig, sigma, lnpval, lnpin, lngesum, err_th

      if (theta_tot.lt.0.d0.or.ntot.lt.0.or.theta_sig.lt.0.d0) then 
        write(*,*) theta_tot, theta_sig, ntot
        stop 'Error: something has gone negative in nulike_pval!'
      endif

      sigma = dsqrt(EAErr(analysis)*EAErr(analysis)+err_th*err_th)
      if (sysErrDist_logNorm(analysis)) then
        !Treat percentage error as log-normal distributed
        call nulike_lnpilnsum(ntot,theta_tot-theta_sig,theta_sig,
     &   sigma,lnpin,lnpval,lngesum)
      else
        !Treat percentage error as Gaussianly-distributed
        call nulike_lnpinsum(ntot,theta_tot-theta_sig,theta_sig,
     &   sigma,lnpin,lnpval,lngesum)
      endif
      nulike_pval = dexp(lnpval)

      end function nulike_pval

