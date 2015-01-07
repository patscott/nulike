!*         -*- mode: fortran -*-
!************************************************************************
!***                           nucommon.h                             ***
!----------------------------------------------------------------------c
!  author: Pat Scott (patscott@physics.mcgill.ca), June 12, 2014

      include "nuconst.h"

      integer analysis, ptypeshare, eventnumshare
      real*8  Eshare, min_detectable_logE
      real*8  events_ee(max_nEvents,max_analyses)

      common /nu_comm/ Eshare, min_detectable_logE, events_ee, analysis, ptypeshare, eventnumshare
      save /nu_comm/

      ! This parameter will be initialized in a block data routine
      ! (see flagblocks.f).  The external line below
      ! ensures the block data routine is loaded at compile time.
      logical credits_rolled
      common /nulike_credit_flag/ credits_rolled
      save /nulike_credit_flag/
      external nulike_credit_flagblock

      interface
        subroutine nulike_bounds(analysis_name, mwimp, annrate, 
     &   nuyield, Nsignal_predicted, NBG_expected, Ntotal_observed, 
     &   lnlike, pvalue, liketype, pvalFromRef, referenceLike, dof,
     &   context)
          use iso_c_binding, only: c_ptr
          implicit none
          include "nuconst.h"
          integer Ntotal_observed, liketype
          real*8 Nsignal_predicted, NBG_expected, lnlike, pvalue, referenceLike, dof, nuyield, mwimp, annrate
          logical pvalFromRef
          character (len=nulike_clen) analysis_name
          type(c_ptr), value :: context
          external nuyield
        end subroutine
      end interface

!*************************** end of nulike.h *****************************
