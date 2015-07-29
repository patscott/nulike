!*         -*- mode: fortran -*-
!************************************************************************
!***                           nuconst.h                              ***
!----------------------------------------------------------------------c
!  author: Pat Scott (p.scott@imperial.ac.uk), June 12, 2014

      integer max_nPrecompE,max_nSensBins,max_nBinsBGAng,max_nBinsBGE,max_nBiasBins
      integer max_nEvents,max_nHistograms,max_ncols,max_analyses,max_threads,nulike_clen
      parameter(nulike_clen     = 300)  !Length of fixed-length strings in nulike
      parameter(max_nPrecompE   = 1000) !Maximum number partial likelihoods per event
      parameter(max_nSensBins   = 50)   !Maximum number of bins in the effective area or volume
      parameter(max_nBiasBins   = 200)  !Maximum number of bins in the biased effective area 
      parameter(max_nBinsBGAng  = 180)  !Maximum number of angular bins in the background
      parameter(max_nBinsBGE    = 100)  !Maximum number of spectral bins in the background
      parameter(max_nEvents     = 5000) !Maximum number of events for any analysis
      parameter(max_nHistograms = 20)   !Max number of energy windows to give the energy dispersion in
      parameter(max_ncols       = 100)  !Max number of columns (entries) in a single energy-window histogram
      parameter(max_analyses    = 5)    !Maximum number of analyses that can be loaded up simultaneously
      parameter(max_threads     = 256)  !Maximum number of threads that nulike can imagine being run with

      integer angular, enrgyest, events, lun, lun2, lun3
      parameter(angular = 1, enrgyest = 2, events = 3, lun = 20, lun2 = 21, lun3 = 22)

      integer Simplex, HyperQuad
      parameter (Simplex = 1, HyperQuad = 2)

      real*8 pi
      parameter (pi=3.141592653589793238d0)

      real*8 effZero, logZero
      parameter(effZero = 1.d-300)
      parameter(logZero = log10(effZero))

      character (len=15) hstring(3)
      parameter(hstring = (/'####--Angular--',
     &                      '####--Nchan--  ',
     &                      '####--Nevents--'/))

!*************************** end of nuconst.h *****************************
