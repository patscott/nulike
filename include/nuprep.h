!*         -*- mode: fortran -*-
!*######################################################################*
!*                       i n c l u d e     f i l e                      *
!*######################################################################*

!************************************************************************
!***                           nuprep.h                               ***
!***         this piece of code is needed as a separate file          ***
!----------------------------------------------------------------------c
!  author: Pat Scott (patscott@physics.mcgill.ca), June 12, 2014


      real*8 events_ee(max_nEvents)
      !real*8 hist_enrgyest(max_nHistograms,max_ncols)

      common /nu_prep_comm/ events_ee
      save /nu_prep_comm/


