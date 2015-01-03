***********************************************************************
*** nulike_psf returns the point spread function for the IceCube detector
*** implied by an input angular error sigma, at some observed and predicted
*** angles relative to some position on the sky.  The psf is rescaled so
*** as to produce unit integrated probability over the range 
***  phi_obs = [0,pi].
***
*** Sigma is interpreted as the 1 sigma error in a single direction of
*** a 2D Gaussian PSF, i.e. the 39.3% confidence value of the angular 
*** deviation (not 68%), and may be given by e.g. the parabaloid sigma
*** for an individual event, or the global mean angular error.
*** 
*** Input:	phi_obs	        observed angle (degrees)
***             phi_pred        predicted angle (degrees)
***             sigma           angular error corresponding to 39.3%
***                             containment angle (degrees)
***          
*** Output:			PSF (degrees^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: May 6, 2011
*** Updated: Jan 2, 2015 (taking into account co-ordinate shift and 
***                       better normalisation.)
***********************************************************************

      real*8 function nulike_psf(phi_obs, phi_pred, sigma)

      use MarcumQ

      implicit none
      include 'nulike.h'

      real*8 phi_obs, phi_pred, pred2, sigma, sigma2, expo 
      real*8 bess, p, q, halfpi2, besi0
      integer ierr
      parameter (halfpi2 = 0.5d0*pi*pi)

      sigma2 = sigma*sigma
      pred2 = phi_pred*phi_pred
      expo = exp(-(pred2 + phi_obs*phi_obs) / (2.d0 * sigma2))
      bess = besi0(phi_pred * phi_obs / sigma2)
      call marcum(1.d0,0.5d0*pred2/sigma2,halfpi2,p,q,ierr)
      nulike_psf = phi_obs * expo * bess / p

      if (ierr .ne. 0) then
        write(*,*) 'phi_obs, phi_pred, sigma:',phi_obs, phi_pred, sigma
        write(*,*) 'ierr = ',ierr
        stop 'Catastrophic error when calling marcum in nulike!'
      endif

      end function nulike_psf
