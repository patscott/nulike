!*         -*- mode: fortran -*-
!************************************************************************
!***                           nuconst.h                              ***
!----------------------------------------------------------------------c
!  author: Pat Scott (p.scott@imperial.ac.uk), June 12, 2014

      integer max_nPrecompE,max_nSensBins,max_nBinsBGAng,max_nBinsBGE
      integer max_nEvents,max_nHistograms,max_ncols,max_analyses,nulike_clen
      parameter(nulike_clen     = 300)  !Length of fixed-length strings in nulike
      parameter(max_nPrecompE   = 1000) !Maximum number partial likelihoods per event
      parameter(max_nSensBins   = 50)   !Maximum number of bins in the effective area or volume
      parameter(max_nBinsBGAng  = 180)  !Maximum number of angular bins in the background
      parameter(max_nBinsBGE    = 100)  !Maximum number of spectral bins in the background
      parameter(max_nEvents     = 2000) !Maximum number of events for any analysis
      parameter(max_nHistograms = 20)   !Max number of energy windows to give the energy dispersion in
      parameter(max_ncols       = 80)   !Max number of columns (entries) in a single energy-window histogram
      parameter(max_analyses    = 15)   !Maximum number of analyses that can be loaded up

      integer angular, enrgyest, events, lun, lun2
      parameter(angular = 1, enrgyest = 2, events = 3, lun = 20, lun2 = 21)

      integer Simplex, HyperQuad
      parameter (Simplex = 1, HyperQuad = 2)

      real*8 pi
      parameter (pi=3.141592653589793238d0)

      character (len=15) hstring(3)
      parameter(hstring = (/'####--Angular--',
     &                      '####--Nchan--  ',
     &                      '####--Nevents--'/))

!*************************** end of nuconst.h *****************************
