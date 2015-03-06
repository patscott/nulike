***********************************************************************
*** nulike_bgangpdf provides the interpolated probability distribution
*** function for the observed angle between the arrival direction
*** of background events and the Solar direction on the sky, normalised
*** to unit integrated probability over the range [0,phi_cut].  
***
*** Input:	cosphi		cos(arrival angle of event)
*** Output:                     pdf (degrees^-1)
***       
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Apr 24, 2011
***********************************************************************

      real*8 function nulike_bgangpdf(cosphi)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      real*8 cosphi, CosToInvDeg, nulike_bgangpdf_a(1)
      integer IER

      call TSVAL1(nBinsBGAng(analysis),BGangdist_phi(:,analysis),
     & BGangdist_prob(:,analysis),BGangdist_derivs(:,analysis),
     & BGangdist_sigma(:,analysis),0,1,cosphi,nulike_bgangpdf_a,IER)

      if (IER .lt. 0) then
        write(*,*) 'TSVAL1 error from background angular'
        write(*,*) 'distribution in nulike_bgangpdf, code:', IER
        stop
      endif

      CosToInvDeg = sqrt(1.d0 - cosphi*cosphi) * pi / 180.d0
      nulike_bgangpdf = nulike_bgangpdf_a(1) * CosToInvDeg / 
     &                  BGangdist_conenorm(analysis)

      end function nulike_bgangpdf
