***********************************************************************
*** nulike_bgangpdf provides the interpolated probability distribution
*** function for the observed angle between the arrival direction
*** of background events and the Solar direction on the sky, normalised
*** to unit integrated probability over the range [0,phi_cut].  
***
*** Input:	cosphi		cos(arrival angle of event)
*** Output:                     pdf (degrees^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
***********************************************************************

      real*8 function nulike_bgangpdf(cosphi)

      implicit none
      include 'nulike.h'

      real*8 cosphi, CosToInvDeg
      integer IER

      call TSVAL1(nBinsBGAng(analysis),BGangdist_phi(:,analysis),
     & BGangdist_prob(:,analysis),BGangdist_derivs(:,analysis),
     & BGangdist_sigma(:,analysis),0,1,cosphi,nulike_bgangpdf,IER)

      if (IER .lt. 0) then
        write(*,*) 'TSVAL1 error from background spectral'
        write(*,*) 'distribution in nulike_bgangpdf, code:', IER
        stop
      endif

      CosToInvDeg = sqrt(1.d0 - cosphi*cosphi) * pi / 180.d0
      nulike_bgangpdf = nulike_bgangpdf * CosToInvDeg / 
     &                  BGangdist_conenorm(analysis)

      end function nulike_bgangpdf
