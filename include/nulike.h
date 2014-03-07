!*         -*- mode: fortran -*-
!*######################################################################*
!*                       i n c l u d e     f i l e                      *
!*######################################################################*

!************************************************************************
!***                           nulike.h                               ***
!***         this piece of code is needed as a separate file          ***
!----------------------------------------------------------------------c
!  author: Pat Scott (patscott@physics.mcgill.ca), March 26 2011

      integer lun, ridiculousNumberOfChannels
      integer max_nBinsEA,max_nBinsBGAng,max_nBinsBGE,max_nEvents
      integer max_nHistograms,max_nnchan,nchan_maxallowed
      integer angular, nchannels, events
      character (len=15) hstring(3), nulike_version

      parameter(lun = 20)
      parameter(ridiculousNumberOfChannels = 1000)
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

      integer nBinsEA,nBinsEAError,nBinsBGAng,nBinsBGE,nEvents
      integer nHistograms
      integer nchan_min, nchan_max, nnchan_total

      real*8 effArea_logE(2,max_nBinsEA)
      real*8 effArea_logEcentres(max_nBinsEA)
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

      real*8 EAlogE_inEAErrBins(2,max_nBinsEA)
      real*8 EAErr_inEAErrBins(max_nBinsEA)
      integer nEvents_inEAErrBins(max_nBinsEA)
      integer maxEAErrIndex

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
      integer bestGuessBin(max_nnchan)
      real*8  relProb(max_nnchan,max_nBinsEA)
      integer nchan_hist2BGoffset
      real*8  edisp_prob(max_nHistograms)
      real*8  edisp_derivs(max_nHistograms)
      real*8  edisp_sigma(max_nHistograms)

      integer events_nchan(max_nBinsEA,max_nEvents)
      real*8  events_cosphi(max_nBinsEA,max_nEvents)
      real*8  events_cosphiErr(max_nBinsEA,max_nEvents)

      real*8  theta_BG(max_nBinsEA)
      real*8  theta_S(max_nBinsEA)
      real*8  theta_Snu(max_nBinsEA)
      real*8  theta_Snubar(max_nBinsEA)
      real*8  phi_max_rad, phi_max_deg, exp_time, theta_S_total
      real*8  Eshare, thetashare, log10mwimp, theoryErr, EAErr_max
      real*8  BGpvalPoissonian, annrate, BGangdist_conenorm
      integer FullSkyBG, ptypeshare, nchanshare, nchansaved
      logical pvalBGPoisComputed, sysErrDist_logNorm

      common /nulike_comm/ events_nchan,events_cosphi,events_cosphiErr,
     & effArea_logE,EAlogE_inEAErrBins,effArea_nu,relProb,theta_BG,
     & theta_S,EAErr_inEAErrBins,nEvents_inEAErrBins, BGnchandist_prob,
     & effArea_logEcentres, effArea_nuderivs, effArea_nubarderivs,
     & effArea_AngResderivs, effArea_nusigma, effArea_nubarsigma,
     & effArea_AngRessigma, BGangdist_derivs, BGangdist_sigma,
     & BGangdist_phi, BGangdist_prob, BGnchandist_nchan, phi_max_rad,
     & phi_max_deg, log10mwimp, theoryErr, EAErr_max, BGpvalPoissonian,
     & BGangdist_norm, effArea_AngRes, exp_time, Eshare, thetashare, 
     & theta_Snu, theta_Snubar, theta_S_total, annrate, 
     & BGangdist_conenorm, hist_LogE, hist_logEcentres, hist_nchan,
     & hist_prob, hist_derivs, hist_sigma, edisp_prob, edisp_derivs,
     & edisp_sigma, effArea_nubar,effArea_syserr,effArea_staterr,
     & bestGuessBin, nBinsEA, nBinsEAError, nBinsBGAng,nBinsBGE,nEvents,
     & nHistograms, nnchan_total, nchan_min, nchan_max, maxEAErrIndex,
     & nchan_hist2BGoffset, FullSkyBG, ptypeshare, nchanshare,
     & nchansaved, pvalBGPoisComputed, sysErrDist_logNorm, nulike_version
      save /nulike_comm/

      ! These parameters will be initialized in a block data routine
      ! (nulike_FlagBlock, see init.f).  The external line below
      ! ensures the block data routine is loaded at compile time.
      logical nulike_init_called     
      common /nulike_flags/ nulike_init_called
      save /nulike_flags/
      external nulike_FlagBlock

!*************************** end of nulike.h *****************************
