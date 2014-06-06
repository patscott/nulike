***********************************************************************
*** nulike_bgpredinit calculates the expected number of background 
*** counts, based on the background spectrum defined in nulike_bginit.
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: April 8, 2011
*** Modified: Jun 3, 6 2014
***********************************************************************

      subroutine nulike_bgpredinit

      implicit none
      include 'nulike.h'

      integer IER
      real*8 TSINTL

      !Work out the total number of background events expected inside phi_cut, as integral 
      !of BG angular distribution from phi = 0 to phi_cut. Then multiply by number of events
      !across full sky, and dole out events into different bins.
      BGangdist_conenorm = TSINTL(dcos(phi_max_rad(analysis)),
     & 1.d0,nBinsBGAng(analysis),dcos(BGangdist_phi),BGangdist_prob,
     & BGangdist_derivs,BGangdist_sigma,IER)
      theta_BG(analysis) = BGangdist_conenorm / BGangdist_norm
     & * dble(FullSkyBG)

      end subroutine nulike_bgpredinit

