***********************************************************************
*** nulike_bgpredinit calculates the expected number of background 
*** counts, based on the background spectrum defined in nulike_bginit.
***        
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: April 8, 2011
*** Modified: Jun 3, 6 2014
***********************************************************************

      subroutine nulike_bgpredinit(cosphimax)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      integer IER
      real*8 cosphimax, TSINTL

      !Die if the cut angle is so small it is fully contained in the first angular bin
      !of the background prediction.
      if (cosphimax .ge. BGangdist_phi(nBinsBGAng(analysis),analysis)) then
        write(*,*) "Analysis '"//trim(analysis_name_array(analysis))//"'"
        write(*,*) 'phi_cut = ',acos(cosphimax)*180.d0/pi,' deg'
        write(*,*) '1st entry of BGangdist_phi = ',acos(BGangdist_phi(nBinsBGAng(analysis),analysis))*180.d0/pi,' deg'
        write(*,*) 'Error: requested cut angle phi_cut is equal to '
        write(*,*) 'or smaller than the centre of the first bin in'
        write(*,*) 'which the angular distribution of the background'
        write(*,*) 'is specified.  Please use a larger value of phi_cut,'
        write(*,*) 'or provide a higher-resolution background datafile.'
        stop
      endif

      !Work out the total number of background events expected inside phi_cut, as integral 
      !of BG angular distribution from phi = 0 to phi_cut. Then multiply by number of events
      !across full sky, and dole out events into different bins.
      BGangdist_conenorm(analysis) = TSINTL(cosphimax,1.d0,
     & nBinsBGAng(analysis),BGangdist_phi(:,analysis),
     & BGangdist_prob(:,analysis),BGangdist_derivs(:,analysis),
     & BGangdist_sigma(:,analysis),IER)
      theta_BG(analysis) = BGangdist_conenorm(analysis) /
     & BGangdist_norm(analysis) * dble(FullSkyBG(analysis))

      end subroutine nulike_bgpredinit

