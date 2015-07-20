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
*** This routine is used only with the 2015 likelihood.
*** 
*** Input:	phi_obs	        observed angle (degrees)
***         phi_pred            predicted angle (degrees)
***         sigma               angular error corresponding to 39.3%
***                             containment angle (degrees)
***          
*** Output:	PSF (degrees^-1)
***       
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Jan 2, 2015 
***********************************************************************

      real*8 function nulike_offctrpsf(phi_obs, phi_pred, sigma)

      use iso_c_binding, only: c_ptr
      use MarcumQ

      implicit none
      include 'nulike_internal.h'

      real*8 phi_obs, phi_pred, sigma, sigma2, expbesi0, expbess  
      real*8 p, q, arg1, arg3
      integer ierr
      
      sigma2 = sigma*sigma
      expbess = expbesi0(phi_pred,phi_obs,sigma2)
      arg1 = 0.5d0*phi_pred*phi_pred/sigma2
      arg3 = 1.62d4/sigma2
      if (arg3 .gt. 1.d4 .or. (arg1 .eq. 0.d0 .and. arg3 .gt. 5.d2)) then
        p = 1.d0
      else
        call marcum(1.d0,arg1,arg3,p,q,ierr)
        if (ierr .gt. 1) then
          write(*,*) 'phi_obs, phi_pred, sigma:',phi_obs, phi_pred, sigma
          write(*,*) 'ierr = ',ierr
          stop 'Catastrophic error when calling marcum in nulike_offctrpsf!'
        endif
      endif

      nulike_offctrpsf = phi_obs * expbess / (p * sigma2)

      end function nulike_offctrpsf
