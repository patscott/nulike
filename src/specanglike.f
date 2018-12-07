***********************************************************************
*** nulike_specanglike returns the contribution of a single event to the
*** unbinned likelihood, based on values of the arrival direction and
*** the energy estimator for the event, as well as the theoretical
*** energy spectrum of incoming neutrinos.
*** This routine is used only with the 2015 likelihood.
***
*** input:  event_number
***         theta_S      total predicted number of signal events
***                       within analysis window (cut cone)
***         f_S          signal fraction; percentage of predicted counts
***                       expected to be due to signal rather than back-
***                       ground.
***         annrate      Annihilation rate (s^-1)
***         logmw        log_10(m_WIMP / GeV)
***         logEmin      log_10(Emin/GeV), where Emin is the lower energy
***                       boundary of the analysis energy range
***
*** output:              ln(Likelihood / chan^-1)
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Jun 8, 2014
***********************************************************************


      double precision function nulike_specanglike(event_number,theta_S,
     & f_S,annrate,logmw,logEmin,nuyield,context,speed)

      use iso_c_binding, only: c_ptr, c_int
      use Precision_Model
      use CUI
      use omp_lib

      implicit none
      include 'nulike_internal.h'

      real*8 theta_S, f_S, annrate, logmw, logEmin, nuyield
      real*8 bgpartial, sigpartial, cosphi, ee, TSINTL
      real*8 nulike_bgangpdf, nulike_bgspec, weight, log10E
      real*8 yvals(max_nPrecompE), yderivs(max_nPrecompE)
      real*8 sigma(max_nPrecompE), working(2*max_nPrecompE-2)
      real*8 SAbsErr, SVertices(1,2), eps, thresh
      real*8 logEmin_true, logEmax_true
      real*8 sigpartial_accurate, measure
      integer event_number, i, ptype, IER, SRgType
      parameter (SRgType = Simplex)
      integer(c_int), intent(in) :: speed
      type(c_ptr) context
      logical, save :: revert_to_accurate_likelihood
      external nuyield

      interface
        function nulike_specangintegrand(NumFun,X) result(Value)
          integer, intent(in) :: NumFun
          real*8, intent(in) :: X(:)
          real*8 :: Value(NumFun)
        end function nulike_specangintegrand
      end interface

      ! Choose integration eps and calculation threshold according to fast likelihood switch
      select case (speed)
      case (0)
        eps = 1d-2
        thresh = 1.d-2
      case (1)
        eps = 0.3d0
        thresh = 1.d-2
      case (2)
        eps = 0.3d0
        thresh = 1.d-1
      case (3)
        thresh = 1.d-1
      case default
        write(*,*) 'Error in nulike_specanglike: unrecognised speed: ', speed
        stop
      end select

      ! Don't bother with energies above or below which there is no estimate of the effective area
      logEmin_true = max(precomp_log10E(1,analysis),logEmin)
      logEmax_true = min(precomp_log10E(nPrecompE(analysis),analysis),logmw)

      ! If this is the first event...
      if (event_number .eq. 1) then
        ! Reset the likelihood accuracy flag
        revert_to_accurate_likelihood = .false.
        ! Set the global context pointers unable to be passed through CUBPACK
        context_shared = context
        nuyield_ptr%f => nuyield
      endif

      ! Retrieve the event info
      cosphi = events_cosphi(event_number,analysis)
      ee = events_ee(event_number,analysis)

      ! Include the background component
      bgpartial = nulike_bgangpdf(cosphi) * nulike_bgspec(ee,2015)

      ! Include the signal component
      if (logEmin_true .gt. logEmax_true) then

        sigpartial = 0.d0

      else

        ! If the signal fraction is big enough to impact the likelihood, go ahead.
        ! Also go ahead if this is the first event and we need to use it to judge
        ! whether to interpolate or integrate in order to get the signal partial likelihood.
        if (f_S/(1.d0-f_S) .gt. thresh * bgpartial .or.
     &      (event_number == 1 .and. speed .ne. 0 .and. speed .ne. 3)) then

          ! Do the likelihood calculation using a fast interpolation
          if (speed .gt. 0 .and. .not. revert_to_accurate_likelihood) then

            yvals = 0.d0
            do i = 1, nPrecompE(analysis)
              log10E = precomp_log10E(i,analysis)
              do ptype = 1,2
                weight = precomp_weights(i,event_number,ptype,analysis)
                if (weight - logZero .gt. epsilon(logZero)) then
                  yvals(i) = yvals(i) + nuyield(log10E,ptype,context) * 10.d0**weight
                endif
              enddo
              yvals(i) = yvals(i) * 10.d0**log10E
            enddo

            !Set up interpolation in spectral-angular likelihood for this event
            call TSPSI(nPrecompE(analysis),precomp_log10E(:,analysis),yvals,
     &       1,0,.false.,.false.,2*nPrecompE(analysis)-2,working,yderivs,
     &       sigma,IER)
            if (IER .lt. 0) then
              write(*,*) 'Error in nulike_specanglike: TSPSI failed with error'
              write(*,*) 'code',IER,' in setting up neutrino effective area.'
              stop
            endif

            !Use the interpolant to calculate the integral
            sigpartial = TSINTL(logEmin_true,logEmax_true,nPrecompE(analysis),
     &       precomp_log10E(:,analysis),yvals,yderivs,sigma,IER)
            if (IER .lt. 0) then
              write(*,*) 'Error in nulike_specanglike: TSPSI failed with error'
              write(*,*) 'code',IER,' in setting up neutrino effective area.'
              stop
            endif

          else

            sigpartial = 0.d0

          endif

          ! Do the likelihood calculation using a proper integration
          if (speed == 0 .or.
     &        revert_to_accurate_likelihood .or.
     &        (event_number == 1 .and. speed .ne. 3)) then

            eventnumshare(omp_get_thread_num()+1) = event_number
            IER = 0
            SVertices(1,:) = (/logEmin_true, logEmax_true/)
            call CUBATR(1,nulike_specangintegrand,SVertices,SRgType,
     &       sigpartial_accurate,SAbsErr,IER,EpsAbs=effZero,MaxPts=5000000,EpsRel=eps,Job=2,Key=2)
            if (IER .ne. 0) then
              write(*,*) 'Error raised by CUBATR in nulike_specanglike: ', IER
              stop
            endif
            call CUBATR()

          else

            sigpartial_accurate = 0.d0

          endif

          ! Revert to the accurate likelihood for this model if speed = 1 or 2 and the fast likelihood
          ! is too far off for the first event.
          if (event_number == 1 .and. speed .ne. 0 .and. speed .ne. 3) then
            measure = f_S * nEvents(analysis) * abs(log(sigpartial) - log(sigpartial_accurate))
            if (measure .gt. 10.d0*thresh) revert_to_accurate_likelihood = .true.
          endif

          if (theta_S .gt. epsilon(theta_S)) then
            if (revert_to_accurate_likelihood .or. speed .eq. 0) sigpartial = sigpartial_accurate
            sigpartial = exp_time(analysis) / theta_S * annrate * dlog(10.d0) * sigpartial
          else
            sigpartial = 0.d0
          endif

        else
        ! If the signal fraction is so small that it can't significantly impact the likelihood,
        ! skip the calculation and just zero for the signal partial likelihood.
        sigpartial = 0.d0

        endif

      endif

      ! Combine background and signal components
      nulike_specanglike = f_S * sigpartial + (1.d0-f_S) * bgpartial

      ! Take ln of total likelihood
      if (nulike_specanglike .le. effZero) then
        nulike_specanglike = logZero
      else
        nulike_specanglike = log(nulike_specanglike)
      endif

      if (nulike_specanglike .ge. 0.d0) then
        write(*,*) 'Error in nulike_specanglike for event ', event_number
        write(*,*) 'Spectral-angular likelihood greater than 1.'
        stop
      endif

      end function nulike_specanglike
