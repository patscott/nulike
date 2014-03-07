***********************************************************************
*** nulike_bgpredinit calculates the expected number of background 
*** counts in each superbin, based on the superbinning defined in
*** nulike_eainit and the background spectrum defined in nulike_bginit.
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: April 8, 2011
***********************************************************************

      subroutine nulike_bgpredinit

      implicit none
      include 'nulike.h'

      integer i, j, IER
      real*8 normFactor, TSINTL

      !Work out the fraction of BG events expected to fall into each superbin
      !This is given by the sum of partial fractions, which are the products of the
      !respective fractions of BG events with each nchan (BGnchandist_prob)
      !and the probability that an event with that nchan would fall into each superbin (relProb).
      theta_BG = 0.d0
      normFactor = 0.d0
      do j = 1, nBinsEAError
        do i = 1, nBinsBGE
          theta_BG(j) = theta_BG(j) + 
     &       relProb(i+nchan_hist2BGoffset,j) * BGnchandist_prob(i)
        enddo
        normFactor = normFactor + theta_BG(j)
      enddo

      !Make sure the final fractions are properly normalised 
      theta_BG = theta_BG / normFactor

      !Work out the total number of background events expected inside phi_cut, as integral 
      !of BG angular distribution from phi = 0 to phi_cut. Then multiply by number of events
      !across full sky, and dole out events into different bins.
      BGangdist_conenorm = TSINTL(dcos(phi_max_rad),1.d0,nBinsBGAng,
     & dcos(BGangdist_phi),BGangdist_prob,BGangdist_derivs,
     & BGangdist_sigma,IER)
      theta_BG = theta_BG * BGangdist_conenorm / BGangdist_norm
     & * dble(FullSkyBG)

      end subroutine nulike_bgpredinit

