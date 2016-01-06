!*         -*- mode: fortran -*-
!************************************************************************
!***                           nucommon.h                             ***
!----------------------------------------------------------------------c
!  author: Pat Scott (p.scott@imperial.ac.uk), June 12, 2014

      include "nuconst.h"

      interface

        subroutine nulike_bounds(analysis_name_in, mwimp, annrate, 
     &   nuyield, Nsignal_predicted, NBG_expected, Ntotal_observed, 
     &   lnlike, pvalue, liketype, theoryError, speed, pvalFromRef, 
     &   referenceLike, dof, context, threadsafe)  BIND(C)
          use iso_c_binding, only: c_ptr, c_char, c_double, c_int, c_bool
          implicit none
          include "nuconst.h"
          integer(c_int), intent(inout) :: Ntotal_observed
          integer(c_int), intent(in) :: liketype, speed
          real(c_double), intent(inout) :: Nsignal_predicted, NBG_expected, lnlike, pvalue
          real(c_double), intent(in) :: referenceLike, dof, mwimp, annrate, theoryError
          logical(c_bool), intent(in) :: pvalFromRef, threadsafe
          character(kind=c_char), dimension(nulike_clen), intent(inout):: analysis_name_in
          type(c_ptr), intent(inout) :: context

          interface
            real(c_double) function nuyield(log10E,ptype,context) bind(c)
              use iso_c_binding, only: c_ptr, c_double, c_int
              implicit none
              real(c_double), intent(in) :: log10E
              integer(c_int), intent(in) :: ptype
              type(c_ptr), intent(inout) :: context
            end function
          end interface

        end subroutine

      end interface

!*************************** end of nulike.h *****************************

