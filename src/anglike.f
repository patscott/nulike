***********************************************************************
*** nulike_anglike returns the contribution of a single event to the
*** unbinned likelihood, based on the reconstructed angle between
*** its arrival direction in the detector and the direction of the Sun.
*** This routine is used only with the 2012 likelihood.
***
*** input:  cosphi     cos(reconstructed angle from solar position)
***         parabsigma paraboloid (or other) sigma corresponding to
***                    1 sigma angular error when considering a single
***                    dimension of a 2D Gaussian PDF (i.e. 39,3% C.L.,
***                    not 68%). (degrees)
***         f_S        signal fraction; percentage of predicted counts
***                    within analysis cut cone and energy range.
***                    expected to be due to signal rather than back-
***                    ground.
*** output:            ln(Likelihood / degrees^-1)
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: May 6, 2011
*** Modified: Jun 3, 6 2014
***********************************************************************

      double precision function nulike_anglike(cosphi,
     & parabsigma,f_S)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      real*8 cosphi, parabsigma, f_S
      real*8 signalpartiallike,bgpartiallike,nulike_psf,nulike_bgangpdf

      !Return low likelihood automatically if cosphi = 1.0 exactly
      if (cosphi .eq. 1.d0) then
        nulike_anglike = logZero
        return
      endif

      !Calculate the signal part of the likelihood
      signalpartiallike = f_S*nulike_psf(acos(cosphi)*180.d0/pi,
     & 0.d0, parabsigma)

      !Calculate the background part of the likelihood
      bgpartiallike = (1.d0-f_S) * nulike_bgangpdf(cosphi)

      !Put them together
      nulike_anglike = log(signalpartiallike + bgpartiallike)

      end function nulike_anglike
