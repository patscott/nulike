***********************************************************************
*** nulike_psf returns the point spread function for the IceCube detector
*** implied by an input angular error sigma, at some observed and predicted
*** angles relative to some position on the sky.  The psf is rescaled so
*** as to produce unit integrated probability over the range [0,phimax]
*** Sigma is interpreted as the 1 sigma error in a single direction of
*** a 2D Gaussian PSF, i.e. the 39.3% confidence value of the angular 
*** deviation (not 68%), and may be given by e.g. the parabaloid sigma
*** for an individual event, or the global mean angular error.
*** This is the 2012 version, only used by the 2012 likelihood.
*** 
*** Input:	phi_obs	        observed angle (degrees)
***             phi_pred        predicted angle (degrees)
***             sigma           angular error corresponding to 39.3%
***                             containment angle (degrees)
***          
*** Output:			PSF (degrees^-1)
***       
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: May 6, 2011
***********************************************************************

      real*8 function nulike_psf(phi_obs, phi_pred, sigma)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      real*8 phi_obs, phi_pred, sigma, sigma2, expo, diff, prefac

      sigma2 = sigma*sigma
      diff = 180.d0 - phi_pred
      expo = exp(-phi_pred*phi_pred / (2.d0 * sigma2))
      expo = expo + exp(-diff*diff / (2.d0 * sigma2))
      prefac = (2.d0-expo)

      diff = phi_obs - phi_pred
      expo = exp(-diff*diff / (2.d0 * sigma2))
      nulike_psf = expo*abs(diff)/(sigma2*prefac) 

      end function nulike_psf
