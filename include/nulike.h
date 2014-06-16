!*         -*- mode: fortran -*-
!*######################################################################*
!*                       i n c l u d e     f i l e                      *
!*######################################################################*

!************************************************************************
!***                           nulike.h                               ***
!***         this piece of code is needed as a separate file          ***
!----------------------------------------------------------------------c
!  author: Pat Scott (patscott@physics.mcgill.ca), March 26 2011, June 3, 6 2014

      include 'nucommon.h'

      integer max_analyses, nchan_maxallowed
      parameter(max_analyses     = 15) !Maximum number of analyses that can be loaded up
      parameter(nchan_maxallowed =100) !Maximum value of nchan

      real*8 bigBadLike
      parameter(bigBadLike = -50.d0)

      character (len=100) analysis_name_array(max_analyses)
      integer likelihood_version(max_analyses)
      logical sysErrDist_logNorm(max_analyses)
      real*8 phi_max_deg(max_analyses)
      real*8 exp_time(max_analyses), theoryErr(max_analyses) 
      real*8 theta_BG(max_analyses), EAErr(max_analyses)

      integer nSensBins(max_analyses), nBinsBGAng(max_analyses)
      integer nEvents(max_analyses), nEvents_in_file(max_analyses)
      integer nHistograms(max_analyses), nnchan_total(max_analyses)
      real*8  ee_min(max_analyses), ee_max(max_analyses) 

      real*8 sens_logE(2,max_nSensBins,max_analyses)
      real*8 sens_logEcentres(max_nSensBins,max_analyses)
      real*8 sens_nu(max_nSensBins,max_analyses)
      real*8 sens_nubar(max_nSensBins,max_analyses)
      real*8 sens_syserr(max_nSensBins,max_analyses)
      real*8 sens_staterr(max_nSensBins,max_analyses)
      real*8 sens_AngRes(max_nSensBins,max_analyses)
      real*8 sens_nuderivs(max_nSensBins,max_analyses)
      real*8 sens_nubarderivs(max_nSensBins,max_analyses)
      real*8 sens_AngResderivs(max_nSensBins,max_analyses)
      real*8 sens_nusigma(max_nSensBins,max_analyses)
      real*8 sens_nubarsigma(max_nSensBins,max_analyses)
      real*8 sens_AngRessigma(max_nSensBins,max_analyses)

      real*8  BGangdist_phi(max_nBinsBGAng,max_analyses)
      real*8  BGangdist_prob(max_nBinsBGAng,max_analyses)
      real*8  BGangdist_derivs(max_nBinsBGAng,max_analyses)
      real*8  BGangdist_sigma(max_nBinsBGAng,max_analyses)
      real*8  BGangdist_norm(max_analyses), BGangdist_conenorm(max_analyses)
      real*8  BGeedist_ee(max_nBinsBGE,max_analyses)
      real*8  BGeedist_prob(max_nBinsBGE,max_analyses)
      integer FullSkyBG(max_analyses)

      real*8  hist_logE(2,max_nHistograms,max_analyses)
      real*8  hist_logEcentres(max_nHistograms,max_analyses)
      integer hist_nchan(max_nHistograms,max_ncols,max_analyses)
      real*8  hist_prob(max_nHistograms,max_ncols,max_analyses)
      real*8  hist_derivs(max_nHistograms,max_ncols,max_analyses)
      real*8  hist_sigma(max_nHistograms,max_ncols,max_analyses)
      integer nchan_hist2BGoffset(max_analyses)

      real*8  events_nchan(max_nEvents,max_analyses)
      real*8  events_cosphi(max_nEvents,max_analyses)
      real*8  events_cosphiErr(max_nEvents,max_analyses)

      logical pvalBGPoisComputed(max_analyses)
      real*8  BGpvalPoissonian(max_analyses)

      integer nPrecompE(max_nPrecompE)
      real*8 precomp_energies(max_nPrecompE,max_analyses)
      real*8 precomp_weights(max_nPrecompE,max_nEvents,2,max_analyses)
      real*8 precomp_derivs(max_nPrecompE,max_nEvents,2,max_analyses)
      real*8 precomp_sigma(max_nPrecompE,max_nEvents,2,max_analyses)
      real*8 precompEA_weights(max_nPrecompE,2,max_analyses)
      real*8 precompEA_derivs(max_nPrecompE,2,max_analyses)
      real*8 precompEA_sigma(max_nPrecompE,2,max_analyses)

      real*8  Eshare, thetashare, annrateshare, nchanshare
      integer ptypeshare, eventnumshare

      common /nulike_comm/ events_nchan,events_cosphi,events_cosphiErr,
     & sens_logE,sens_nu,theta_BG, BGeedist_prob,
     & sens_logEcentres, sens_nuderivs, sens_nubarderivs,
     & sens_AngResderivs, sens_nusigma, sens_nubarsigma,
     & sens_AngRessigma, BGangdist_derivs, BGangdist_sigma,
     & BGangdist_phi, BGangdist_prob, BGeedist_ee,
     & phi_max_deg, theoryErr, EAErr, BGpvalPoissonian, BGangdist_norm,
     & sens_AngRes, exp_time, Eshare, thetashare, annrateshare, 
     & BGangdist_conenorm, hist_LogE, hist_logEcentres, hist_nchan,
     & hist_prob, hist_derivs, hist_sigma, sens_nubar,sens_syserr,
     & sens_staterr, precomp_energies, precomp_weights, precomp_derivs,
     & precomp_sigma, precompEA_weights, precompEA_derivs, 
     & precompEA_sigma, nSensBins, nBinsBGAng, nEvents, nEvents_in_file,
     & nPrecompE, nHistograms, nnchan_total, ee_min, ee_max,
     & nchanshare, nchan_hist2BGoffset, FullSkyBG, ptypeshare,
     & pvalBGPoisComputed, sysErrDist_logNorm, eventnumshare,
     & analysis_name_array, likelihood_version
      save /nulike_comm/

      ! This parameter will be initialized in a block data routine
      ! (see flagblocks.f).  The external line below ensures the 
      ! block data routine is loaded at compile time.
      logical nulike_init_called     
      common /nulike_init_flag/ nulike_init_called
      save /nulike_init_flag/
      external nulike_init_flagblock

!*************************** end of nulike.h *****************************
