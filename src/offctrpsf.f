***********************************************************************
*** nulike_offctrpsf returns the point spread function for the IceCube 
*** detector implied by an input angular error sigma, at some observed 
*** and predicted angles relative to some position on the sky.  The psf
*** takes into account the offset between the origin of the observing 
*** frame and the origin of the PSF frame.  It is normalised so as to 
*** produce unit integrated probability over the range phi_obs = [0,pi].
***
*** Sigma is interpreted as the 1 sigma error in a single direction of
*** a 2D Gaussian PSF, i.e. the 39.3% confidence value of the angular 
*** deviation (not 68%), and may be given by e.g. the parabaloid sigma
*** for an individual event, or the global mean angular error.
***
*** This routine is used only with the 2014 likelihood.
*** 
*** Input:	phi_obs	        observed angle (degrees)
***             phi_pred        predicted angle (degrees)
***             sigma           angular error corresponding to 39.3%
***                             containment angle (degrees)
***          
*** Output:			PSF (degrees^-1)
***       
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Jan 2, 2015 
***********************************************************************

      real*8 function nulike_offctrpsf(phi_obs, phi_pred, sigma)

      use MarcumQ

      implicit none
      include 'nulike.h'

      real*8 phi_obs, phi_pred, pred2, sigma, sigma2, expo 
      real*8 bess, p, q, besi0, arg3
      integer ierr

      sigma2 = sigma*sigma
      pred2 = phi_pred*phi_pred
      expo = exp(-(pred2 + phi_obs*phi_obs) / (2.d0 * sigma2))
      bess = besi0(phi_pred * phi_obs / sigma2)
      arg3 = 1.62d4/sigma2
      if (arg3 .lt. 1.d5) then
        call marcum(1.d0,0.5d0*pred2/sigma2,arg3,p,q,ierr)
      else
        p = 1.d0
      endif
      nulike_offctrpsf = phi_obs * expo * bess / p

      if (ierr .gt. 1) then
        write(*,*) 'phi_obs, phi_pred, sigma:',phi_obs, phi_pred, sigma
        write(*,*) 'ierr = ',ierr
        stop 'Catastrophic error when calling marcum in nulike_offctrpsf!'
      endif

      end function nulike_offctrpsf
