!*         -*- mode: fortran -*-
!************************************************************************
!***                           nucommon.h                             ***
!----------------------------------------------------------------------c
!  author: Pat Scott (p.scott@imperial.ac.uk), June 12, 2014

      include "nuconst.h"

      integer analysis, ptypeshare
      integer eventnumshare(max_threads)
      real*8  Eshare, min_detectable_logE
      real*8  events_ee(max_nEvents,max_analyses)
      real*8  EAErr(max_analyses)

      common /nu_comm/ Eshare, min_detectable_logE, events_ee, EAErr,
     & analysis, ptypeshare, eventnumshare
      save /nu_comm/

      ! This parameter will be initialized in a block data routine
      ! (see flagblocks.f).  The external line below
      ! ensures the block data routine is loaded at compile time.
      logical credits_rolled
      common /nulike_credit_flag/ credits_rolled
      save /nulike_credit_flag/
      external nulike_credit_flagblock

!*************************** end of nucommon.h **************************
