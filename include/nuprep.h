!*         -*- mode: fortran -*-
!*######################################################################*
!*                       i n c l u d e     f i l e                      *
!*######################################################################*

!************************************************************************
!***                           nuprep.h                               ***
!***         this piece of code is needed as a separate file          ***
!----------------------------------------------------------------------c
!  author: Pat Scott (p.scott@imperial.ac.uk), June 12, 2014


      real*8 m_water !Mass of a single water molecule, in kg
      parameter (m_water = 18.01528d0 * 1.66053892d-27)

      real*8 eps_partials !Integration accuracy for partial likelihoods
      parameter (eps_partials = 1.d-3)

      real*8 m_p, m_n, lepmass(3) !Proton, neutron and lepton masses (GeV)
      parameter (m_p = 0.938272d0)
      parameter (m_n = 0.939566d0)
      parameter (lepmass = (/0.510999d-3,0.105658d0,1.77682d0/))

      real*8 phi_max, numdens_n, numdens_p
      real*8 hist_ee_flip(max_ncols,max_nHistograms+2)
      real*8 hist_prob_flip(max_ncols,max_nHistograms+2)
      real*8 hist_derivs_flip(max_ncols,max_nHistograms+2)
      real*8 hist_sigma_flip(max_ncols,max_nHistograms+2)
      real*8 hist_single_ee_prob(max_nHistograms+2)
      real*8 hist_single_ee_derivs(max_nHistograms+2)
      real*8 hist_single_ee_sigma(max_nHistograms+2)
      real*8 hist_logEnergies(max_nHistograms+2)

      integer nhist, leptypeshare

      abstract interface
        real*8 function func (E,x,y,nu,targ,interaction)
          real*8 E, x, y
          integer nu
          character*1 targ
          character*2 interaction
        end function func
      end interface
      type dsdxdy_ptr_container
        SEQUENCE
        procedure (func), pointer, nopass :: f
      end type dsdxdy_ptr_container
      type (dsdxdy_ptr_container) :: dsdxdy_ptr

      common /nu_prep_comm/ hist_ee_flip, hist_prob_flip,
     & hist_derivs_flip, hist_sigma_flip, hist_single_ee_prob,
     & hist_single_ee_derivs, hist_single_ee_sigma, hist_logEnergies,
     & phi_max, numdens_n, numdens_p, nhist, leptypeshare, dsdxdy_ptr
      save /nu_prep_comm/


