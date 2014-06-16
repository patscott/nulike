!*         -*- mode: fortran -*-
!*######################################################################*
!*                       i n c l u d e     f i l e                      *
!*######################################################################*

!************************************************************************
!***                           nuprep.h                               ***
!***         this piece of code is needed as a separate file          ***
!----------------------------------------------------------------------c
!  author: Pat Scott (patscott@physics.mcgill.ca), June 12, 2014


      real*8 phi_max
      real*8 events_ee(max_nEvents)
      real*8 hist_ee_flip(max_ncols,max_nHistograms)
      real*8 hist_prob_flip(max_ncols,max_nHistograms)
      real*8 hist_derivs_flip(max_ncols,max_nHistograms)
      real*8 hist_sigma_flip(max_ncols,max_nHistograms)
      real*8 hist_single_ee_prob(max_nHistograms)
      real*8 hist_single_ee_derivs(max_nHistograms)
      real*8 hist_single_ee_sigma(max_nHistograms)
      real*8 hist_logEnergies(max_nHistograms)

      integer nhist

      common /nu_prep_comm/ events_ee, hist_ee_flip, hist_prob_flip,
     & hist_derivs_flip, hist_sigma_flip, hist_single_ee_prob,
     & hist_single_ee_derivs, hist_single_ee_sigma, hist_logEnergies,
     & phi_max, nhist
      save /nu_prep_comm/


