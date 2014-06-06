!*         -*- mode: fortran -*-
!*######################################################################*
!*                       i n c l u d e     f i l e                      *
!*######################################################################*

!************************************************************************
!***                           nulike.h                               ***
!***         this piece of code is needed as a separate file          ***
!----------------------------------------------------------------------c
!  author: Pat Scott (patscott@physics.mcgill.ca), March 26 2011, June 3, 6 2014

      integer lun, ridiculousNumberOfChannels
      integer max_nBinsEA,max_nBinsBGAng,max_nBinsBGE,max_nEvents
      integer max_nHistograms,max_nnchan,nchan_maxallowed,max_analyses
      integer angular, nchannels, events
      character (len=15) hstring(3)

      parameter(lun = 20)
      parameter(ridiculousNumberOfChannels = 1000)
      parameter(max_analyses    = 15)
      parameter(max_nBinsEA     = 10)
      parameter(max_nBinsBGAng  = 180)
      parameter(max_nBinsBGE    = 20)
      parameter(max_nEvents     = 2000)
      parameter(max_nHistograms = 20)
      parameter(max_nnchan      = 23) !Max number of different nchan values for 2012 likelihood
      parameter(nchan_maxallowed=100) !Maximum value of nchan
      parameter(angular = 1, nchannels = 2, events = 3)
      parameter(hstring = (/'####--Angular--',
     &                      '####--Nchan--  ',
     &                      '####--Nevents--'/))

      real*8 pi, bigBadLike
      parameter (pi=3.141592653589793238d0)
      parameter(bigBadLike = -50.d0)

      character (len=100) analysis_name_array(max_analyses)
      integer likelihood_version(max_analyses)
      logical sysErrDist_logNorm(max_analyses)
      real*8 phi_max_rad(max_analyses), phi_max_deg(max_analyses)
      real*8 exp_time(max_analyses), theoryErr(max_analyses) 
      real*8 theta_BG(max_analyses)

      integer nBinsEA(max_analyses)
      integer nEvents,nEvents_in_file
      integer nBinsBGAng,nBinsBGE,nHistograms
      integer nchan_min, nchan_max, nnchan_total
      integer analysis, nAnalyses

      real*8 effArea_logE(max_analyses,2,max_nBinsEA)
      real*8 effArea_logEcentres(max_analyses,max_nBinsEA)
      real*8 effArea_nu(max_nBinsEA)
      real*8 effArea_nubar(max_nBinsEA)
      real*8 effArea_syserr(max_nBinsEA)
      real*8 effArea_staterr(max_nBinsEA)
      real*8 effArea_AngRes(max_nBinsEA)
      real*8 effArea_nuderivs(max_nBinsEA)
      real*8 effArea_nubarderivs(max_nBinsEA)
      real*8 effArea_AngResderivs(max_nBinsEA)
      real*8 effArea_nusigma(max_nBinsEA)
      real*8 effArea_nubarsigma(max_nBinsEA)
      real*8 effArea_AngRessigma(max_nBinsEA)

      real*8  BGangdist_phi(max_nBinsBGAng)
      real*8  BGangdist_prob(max_nBinsBGAng)
      real*8  BGangdist_derivs(max_nBinsBGAng)
      real*8  BGangdist_sigma(max_nBinsBGAng)
      real*8  BGangdist_norm
      integer BGnchandist_nchan(max_nBinsBGE)
      real*8  BGnchandist_prob(max_nBinsBGE)

      real*8  hist_logE(2,max_nHistograms)
      real*8  hist_logEcentres(max_nHistograms)
      integer hist_nchan(max_nHistograms, max_nnchan)
      real*8  hist_prob(max_nHistograms, max_nnchan)
      real*8  hist_derivs(max_nHistograms, max_nnchan)
      real*8  hist_sigma(max_nHistograms, max_nnchan)
      integer nchan_hist2BGoffset
      real*8  edisp_prob(max_nHistograms)
      real*8  edisp_derivs(max_nHistograms)
      real*8  edisp_sigma(max_nHistograms)

      integer events_nchan(max_nEvents)
      real*8  events_cosphi(max_nEvents)
      real*8  events_cosphiErr(max_nEvents)

      real*8  theta_S,theta_Snu,theta_Snubar
      real*8  Eshare, thetashare, log10mwimp, EAErr
      real*8  BGpvalPoissonian, annrate, BGangdist_conenorm
      integer FullSkyBG, ptypeshare, nchanshare, nchansaved
      logical pvalBGPoisComputed

      common /nulike_comm/ events_nchan,events_cosphi,events_cosphiErr,
     & effArea_logE,effArea_nu,theta_BG, theta_S, BGnchandist_prob,
     & effArea_logEcentres, effArea_nuderivs, effArea_nubarderivs,
     & effArea_AngResderivs, effArea_nusigma, effArea_nubarsigma,
     & effArea_AngRessigma, BGangdist_derivs, BGangdist_sigma,
     & BGangdist_phi, BGangdist_prob, BGnchandist_nchan, phi_max_rad,
     & phi_max_deg, log10mwimp, theoryErr, EAErr, BGpvalPoissonian,
     & BGangdist_norm, effArea_AngRes, exp_time, Eshare, thetashare, 
     & theta_Snu, theta_Snubar, annrate, 
     & BGangdist_conenorm, hist_LogE, hist_logEcentres, hist_nchan,
     & hist_prob, hist_derivs, hist_sigma, edisp_prob, edisp_derivs,
     & edisp_sigma, effArea_nubar,effArea_syserr,effArea_staterr,
     & nBinsEA, nBinsBGAng, nBinsBGE, nEvents, nEvents_in_file,
     & nHistograms, nnchan_total, nchan_min, nchan_max,
     & nchan_hist2BGoffset, FullSkyBG, ptypeshare, nchanshare,
     & nchansaved, pvalBGPoisComputed, sysErrDist_logNorm,
     & analysis, nAnalyses, analysis_name_array, likelihood_version
      save /nulike_comm/

      ! These parameters will be initialized in a block data routine
      ! (nulike_FlagBlock, see init.f).  The external line below
      ! ensures the block data routine is loaded at compile time.
      logical nulike_init_called     
      common /nulike_flags/ nulike_init_called
      save /nulike_flags/
      external nulike_FlagBlock

!*************************** end of nulike.h *****************************
