!*         -*- mode: fortran -*-
!*######################################################################*
!*                       i n c l u d e     f i l e                      *
!*######################################################################*

!************************************************************************
!***                           nulike.h                               ***
!***         this piece of code is needed as a separate file          ***
!----------------------------------------------------------------------c
!  author: Pat Scott (p.scott@imperial.ac.uk), March 26 2011, June 3, 6 2014

      include 'nucommon.h'

      integer nchan_maxallowed
      parameter(nchan_maxallowed = 100) !Maximum value of nchan

      character (len=nulike_clen) analysis_name_array(max_analyses)
      integer likelihood_version(max_analyses)
      logical sysErrDist_logNorm(max_analyses)
      logical no_bias(max_analyses)
      real*8 phi_max_deg(max_analyses)
      real*8 exp_time(max_analyses), theoryErr(max_analyses)
      real*8 theta_BG(max_analyses)

      integer nSensBins(max_analyses), nBiasBins(max_analyses)
      integer nBinsBGAng(max_analyses), nBinsBGE(max_analyses)
      integer nEvents(max_analyses)
      integer nHistograms(max_analyses), nnchan_total(max_analyses)
      real*8  ee_min(max_analyses), ee_max(max_analyses)

      real*8 sens_logE(2,max_nSensBins+1,max_analyses)
      real*8 sens_logEcentres(max_nSensBins+1,max_analyses)
      real*8 sens_nu(max_nSensBins+1,max_analyses)
      real*8 sens_nubar(max_nSensBins+1,max_analyses)
      real*8 sens_syserr(max_nSensBins+1,max_analyses)
      real*8 sens_staterr(max_nSensBins+1,max_analyses)
      real*8 sens_AngRes(max_nSensBins+1,max_analyses)
      real*8 sens_nuderivs(max_nSensBins+1,max_analyses)
      real*8 sens_nubarderivs(max_nSensBins+1,max_analyses)
      real*8 sens_AngResderivs(max_nSensBins+1,max_analyses)
      real*8 sens_nusigma(max_nSensBins+1,max_analyses)
      real*8 sens_nubarsigma(max_nSensBins+1,max_analyses)
      real*8 sens_AngRessigma(max_nSensBins+1,max_analyses)

      real*8 bias_logE(2,max_nBiasBins+1,max_analyses)
      real*8 bias_logEcentres(max_nBiasBins+1,max_analyses)
      real*8 bias_nu(max_nBiasBins+1,max_analyses)
      real*8 bias_nubar(max_nBiasBins+1,max_analyses)
      real*8 bias_nuderivs(max_nBiasBins+1,max_analyses)
      real*8 bias_nubarderivs(max_nBiasBins+1,max_analyses)
      real*8 bias_nusigma(max_nBiasBins+1,max_analyses)
      real*8 bias_nubarsigma(max_nBiasBins+1,max_analyses)

      real*8  BGangdist_phi(max_nBinsBGAng,max_analyses)
      real*8  BGangdist_prob(max_nBinsBGAng,max_analyses)
      real*8  BGangdist_derivs(max_nBinsBGAng,max_analyses)
      real*8  BGangdist_sigma(max_nBinsBGAng,max_analyses)
      real*8  BGangdist_norm(max_analyses), BGangdist_conenorm(max_analyses)
      real*8  BGeedist_ee(max_nBinsBGE,max_analyses)
      real*8  BGeedist_prob(max_nBinsBGE,max_analyses)
      real*8  BGeedist_derivs(max_nBinsBGE,max_analyses)
      real*8  BGeedist_sigma(max_nBinsBGE,max_analyses)
      integer FullSkyBG(max_analyses)

      real*8  hist_logE(2,max_nHistograms+2,max_analyses)
      real*8  hist_logEcentres(max_nHistograms+2,max_analyses)
      integer hist_nchan(max_nHistograms+2,max_ncols,max_analyses)
      real*8  hist_prob(max_nHistograms+2,max_ncols,max_analyses)
      real*8  hist_derivs(max_nHistograms+2,max_ncols,max_analyses)
      real*8  hist_sigma(max_nHistograms+2,max_ncols,max_analyses)
      integer nchan_hist2BGoffset(max_analyses)

      real*8  events_nchan(max_nEvents,max_analyses)
      real*8  events_cosphi(max_nEvents,max_analyses)
      real*8  events_cosphiErr(max_nEvents,max_analyses)

      logical pvalBGPoisComputed(max_analyses)
      real*8  BGpvalPoissonian(max_analyses)

      integer nPrecompE(max_nPrecompE)
      integer start_index(max_analyses), start_index_noL(max_analyses)
      real*8 precomp_log10E(max_nPrecompE,max_analyses)
      real*8 precomp_weights(max_nPrecompE,max_nEvents,2,max_analyses)
      real*8 precomp_derivs(max_nPrecompE,max_nEvents,2,max_analyses)
      real*8 precomp_sigma(max_nPrecompE,max_nEvents,2,max_analyses)
      real*8 precompEA_weights(max_nPrecompE,2,max_analyses)
      real*8 precompEA_derivs(max_nPrecompE,2,max_analyses)
      real*8 precompEA_sigma(max_nPrecompE,2,max_analyses)
      real*8 precompEAnoL_weights(max_nPrecompE,2,max_analyses)
      real*8 precompEAnoL_derivs(max_nPrecompE,2,max_analyses)
      real*8 precompEAnoL_sigma(max_nPrecompE,2,max_analyses)

      real*8  thetashare, annrateshare, nchanshare

      type(c_ptr) context_shared

      abstract interface
        real(c_double) function nuyield_signature(log10E,ptype,context)
          use iso_c_binding, only: c_ptr, c_double, c_int
          implicit none
          real(c_double), intent(in) :: log10E
          integer(c_int), intent(in) :: ptype
          type(c_ptr), intent(inout) :: context
        end function nuyield_signature
      end interface
      type nuyield_ptr_container
        SEQUENCE
        procedure (nuyield_signature), pointer, nopass :: f
      end type nuyield_ptr_container
      type (nuyield_ptr_container) :: nuyield_ptr

      common /nulike_comm/ context_shared, nuyield_ptr,
     & events_nchan, events_cosphi,
     & events_cosphiErr, sens_logE,sens_nu,theta_BG, BGeedist_prob,
     & sens_logEcentres, sens_nuderivs, sens_nubarderivs,
     & sens_AngResderivs, sens_nusigma, sens_nubarsigma,
     & sens_AngRessigma, BGangdist_derivs, BGangdist_sigma,
     & BGangdist_phi, BGangdist_prob, BGeedist_ee,
     & phi_max_deg, theoryErr, BGpvalPoissonian, BGangdist_norm,
     & sens_AngRes, exp_time, thetashare, annrateshare,
     & BGangdist_conenorm, hist_LogE, hist_logEcentres, hist_nchan,
     & hist_prob, hist_derivs, hist_sigma, sens_nubar,sens_syserr,
     & sens_staterr, precomp_log10E, precomp_weights, precomp_derivs,
     & precomp_sigma, precompEAnoL_weights, precompEAnoL_derivs,
     & precompEAnoL_sigma, precompEA_weights, precompEA_derivs,
     & precompEA_sigma, bias_logE, bias_logEcentres, bias_nu,
     & bias_nubar, bias_nuderivs, bias_nubarderivs, bias_nusigma,
     & bias_nubarsigma, ee_min, ee_max, nchanshare,
     & BGeedist_derivs, BGeedist_sigma,
     & nBinsBGE, nBinsBGAng, nEvents, start_index, start_index_noL,
     & nSensBins, nBiasBins, nPrecompE, nHistograms, nnchan_total,
     & nchan_hist2BGoffset, FullSkyBG,
     & pvalBGPoisComputed, sysErrDist_logNorm, no_bias,
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
