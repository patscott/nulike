***********************************************************************
*** nulike_bounds returns counts, likelihoods and p-values from neutrino
*** telescope searches for dark matter annihilation in the Sun.
***
*** input: analysis_name Name of the nulike analysis to use for calculating
***                      likelihood and/or p-value.
***
***        mwimp         WIMP mass (GeV)
***
***        annrate       Annihilation rate (s^-1)
***
***        nuyield       Name of a function that takes arguments
***                       real*8   log10E  log_10(E_nu / GeV)
***                       integer  ptype   1=nu, 2=nubar
***                       c_ptr    context A void C pointer that can be used for callback
***                      and returns
***                       real*8   differential neutrino flux at the detector
***                                (m^-2 GeV^-1 annihilation^-1)
***
***        liketype      Sets combination of data to use in likelihood
***                      calculations.  Note that liketype = 2 and 3 only
***                      have any meaning when using the 2012 likelihood.
***                        1 => Number of events only
***                        2 => Number of events and event arrival angles
***                        3 => Number of events and energy estimator (nchan
***                             = number of lit DOMs)
***                        4 => Number of events, event arrival angles and
***                             energy estimator
***
***        theoryError   theoretical error to incorporate into likelihood
***                      and p-value calculations (given as a fractional
***                      relative error).
***
***        speed         In the case of liketype=4 with the 2015 likelihood,
***                      take some shortcuts on the spectral-angular
***                      part of the likelihood.
***                      0 => no shortcuts (slow)
***                           -full integration of partial likelihood (not log) with 1% accuracy
***                           -neglect signal contribution to spectral-anglular likelihood
***                            when it cannot be > 1% of the likelihood.
***                      1 => fast
***                           -use interpolated likelihood where it does not introduce an
***                            error of more than 0.1 in final spectral-angular log-likelihood
***                           -otherwise use full integration of partial likelihood with 30%
***                            accuracy in likelihood (not log)
***                           -neglect signal contribution to spectral-anglular likelihood
***                            when it cannot be > 1% of the likelihood.
***                      2 => faster
***                           -use interpolated likelihood where it does not introduce an
***                            error of more than 1.0 in final spectral-angular log-likelihood
***                           -otherwise use full integration of partial likelihood with 30%
***                            accuracy in likelihood (not log)
***                           -neglect signal contribution to spectral-anglular likelihood
***                            when it cannot be > 10% of the likelihood.
***                      3 => fastest
***                           -always use interpolated likelihood,
***                           -neglect signal contribution to spectral-anglular likelihood
***                            when it cannot be > 10% of the likelihood.
***
***        pvalFromRef   T => calculate the p-value with reference to a user-
***                           specified likelihood, assuming that the log-
***                           likelihood follows a chi^2 distriubution in
***                           the vicinity of the reference point (i.e.
***                           essentially assume that the reference corresponds
***                           to the best-fit log likelihood, and that Wilks'
***                           theorem is satisfied).
***                      F => calculate the p-value with reference to the
***                           observed background
***
***        referenceLike reference value of the (natural) log likelihood to
***                      use for the calculation of the p-value if
***                      pvalFromRef=T; ignored otherwise.
***
***        dof           degrees of freedom to assume in the calculation of the
***                      p-value in the case that pvalFromRef=T; ignored
***                      otherwise.  This should normally be set to the
***                      number of parameters over which the likelihood has
***                      not been maximised when the solution is perturbed
***                      from the minimum (i.e. the number of free parameters
***                      in a scan minus the number profiled over).
***
***        context       A c_ptr passed in to nuyield when it is called
***
***        threadsafe    Indicate that nulike can treat nuyield as threadsafe.
***
***
*** output: Nsignal_predicted  Predicted number of signal events within
***                      phi_cut (includes solar coronal BG)
***
***         NBG_expected Expected number of (non-solar) background
***                      events insided phi_cut
***
***         Ntotal_observed  Observed number of events inside phi_cut
***
***         lnlike       natural log of chosen likelihood
***
***         pvalue       derived p-value for chosen model
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Mar 20, 2011
*** Updated: Jul 21, 2011
***          Mar 6, 2014
***          Jun 3, 6, 8 2014
***********************************************************************


      subroutine nulike_bounds(analysis_name_in, mwimp, annrate,
     & nuyield, Nsignal_predicted, NBG_expected, Ntotal_observed,
     & lnlike, pvalue, liketype, theoryError, speed, pvalFromRef,
     & referenceLike, dof, context, threadsafe) bind(c)

      use iso_c_binding, only: c_ptr, c_char, c_double, c_int, c_bool
      use omp_lib

      implicit none
      include 'nulike_internal.h'

      integer(c_int), intent(inout) :: Ntotal_observed
      integer(c_int), intent(in) :: liketype, speed
      real(c_double), intent(inout) :: Nsignal_predicted, NBG_expected, lnlike, pvalue
      real(c_double), intent(in) :: referenceLike, dof, mwimp, annrate, theoryError
      logical(c_bool), intent(in) :: pvalFromRef, threadsafe
      character(kind=c_char), dimension(nulike_clen), intent(inout) :: analysis_name_in
      type(c_ptr), intent(inout) :: context

      integer j
      integer counted1, counted2, countrate, nulike_amap
      real*8 nulike_pval, theta_S, DGAMIC, DGAMMA
      real*8 nLikelihood, angularLikelihood, spectralLikelihood, logmw
      real*8 theta_tot, f_S, nulike_anglike, nulike_speclike, nulike_nlike
      real*8 deltalnlike, specAngLikelihood, nulike_signal, nulike_specanglike
      logical nulike_speclike_reset, doProfiling
      character(len=nulike_clen) analysis_name
      !Hidden option for doing speed profiling
      parameter (doProfiling = .false.)

      interface
        real(c_double) function nuyield(log10E,ptype,context) bind(c)
          use iso_c_binding, only: c_ptr, c_double, c_int
          implicit none
          real(c_double), intent(in) :: log10E
          integer(c_int), intent(in) :: ptype
          type(c_ptr), intent(inout) :: context
        end function
      end interface


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 1. Initialisation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      analysis_name = analysis_name_in(1)
      do j = 1, nulike_clen
        analysis_name = analysis_name(1:j-1)//analysis_name_in(j)
      enddo

      if (doProfiling) call system_clock(counted1,countrate)

      !Check validity of user-provided liketype
      if (liketype .lt. 1 .or. liketype .gt. 4) then
        write(*,*) 'Invalid liketype provided to nulike_bounds.'
        write(*,*) 'Quitting...'
        stop
      endif

      !If nulike_init has not yet been called, quit.
      if (.not. nulike_init_called) then
        write(*,*) "Please call nulike_init before nulike_bounds."
        write(*,*) 'Quitting...'
        stop
      endif

      !Look up the analysis requested by the user.
      analysis = nulike_amap(analysis_name)
      if (analysis .eq. 0) then
        write(*,*) "Analysis '"//trim(analysis_name)//"' requested of nulike_bounds"
        write(*,*) 'is not one of the ones that has already been loaded.'
        write(*,*) 'Quitting...'
        stop
      endif

      !Make sure the user has not tried to use the 2015 like with only angular or only spectral likelihood.
      if (likelihood_version(analysis) .eq. 2015 .and.
     & (liketype .eq. 2 .or. liketype .eq. 3) ) then
        write(*,*) "Analysis '"//trim(analysis_name)//"' requested of nulike_bounds"
        write(*,*) 'uses the 2015 likelihood, which is incompatible with liketype = ',liketype
        write(*,*) 'Quitting...'
        stop
      endif

      !Take log of WIMP mass
      logmw = dlog10(mwimp)

      if (doProfiling) then
        call system_clock(counted2,countrate)
        write(*,*) 'Elapsed time initialising (s): ',
     &   real(counted2 - counted1)/real(countrate)
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 2. Signal calculation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Calculate signal counts and spectrum.
      theta_S = nulike_signal(nuyield, context, annrate, logmw, likelihood_version(analysis))
      !Calculate the total predicted number of events
      theta_tot = theta_BG(analysis) + theta_S
      !Calculate the signal fraction.
      f_S = theta_S / theta_tot

      if (doProfiling) then
        call system_clock(counted1,countrate)
        write(*,*) 'Elapsed time on signal calc (s): ',
     &   real(counted1 - counted2)/real(countrate)
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 3. Number likelihood
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Calculate the number likelihood
      nLikelihood = nulike_nlike(nEvents(analysis),
     & theta_tot,theta_S,EAErr(analysis),theoryError)

      if (doProfiling) then
        call system_clock(counted1,countrate)
        write(*,*) 'Elapsed time on number likelihood calc (s): ',
     &   real(counted1 - counted2)/real(countrate)
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 4. Angular and/or spectral likelihood
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Switch according to likelihood version.
      select case (likelihood_version(analysis))

      ! 2012 likelihood, as per arXiv:1207.0810
      case (2012)
        angularLikelihood = 0.d0
        spectralLikelihood = 0.d0
        if (liketype .ne. 1) then
          !Reset the saved spectral likelihoods next time nulike_speclike is run
          nulike_speclike_reset = .true.
          !Step through the individual events
          do j = 1, nEvents(analysis)
            !Add in angular likelihood for this event
            if (liketype .eq. 2 .or. liketype .eq. 4)
     &       angularLikelihood = angularLikelihood + nulike_anglike(
     &       events_cosphi(j,analysis),events_cosphiErr(j,analysis),f_S)
            !Add in spectral likelihood for this event
            if (liketype .eq. 3 .or. liketype .eq. 4)
     &       spectralLikelihood = spectralLikelihood + nulike_speclike(
     &       events_nchan(j,analysis),theta_S,f_S,annrate,logmw,
     &       nulike_speclike_reset,
     &       sens_logE(1,1,analysis),
     &       sens_logE(2,nSensBins(analysis),analysis),
     &       nuyield, context)
          enddo
        endif
        specAngLikelihood = angularLikelihood + spectralLikelihood

      !2015 likelihood, as per arXiv:1601.00653
      case (2015)
        specAngLikelihood = 0.d0
        if (liketype .eq. 4) then
          if (threadsafe .and. omp_get_max_threads() .gt. max_threads) then
            write(*,*) "Nulike actually hasn't been configured to deal with"
            write(*,*) "such a powerhouse machine.  Please increase "
            write(*,*) "max_threads in nuconst.h."
            stop
          endif
          specAngLikelihood = nulike_specanglike(1,theta_S, f_S, annrate,
     &     logmw, sens_logE(1,1,analysis), nuyield, context, speed)
          if (threadsafe) then
!$omp parallel do reduction(+:specAngLikelihood)
            do j = 2, nEvents(analysis)
              specAngLikelihood = specAngLikelihood + nulike_specanglike(j,
     &         theta_S, f_S, annrate, logmw, sens_logE(1,1,analysis),
     &         nuyield, context, speed)
            enddo
          else
            do j = 2, nEvents(analysis)
              specAngLikelihood = specAngLikelihood + nulike_specanglike(j,
     &         theta_S, f_S, annrate, logmw, sens_logE(1,1,analysis),
     &         nuyield, context, speed)
            enddo
          endif
        endif

      case default
        write(*,*) "Unrecognised likelihood version in nulike_bounds."
        write(*,*) "Quitting..."
        stop

      end select

      if (doProfiling) then
        call system_clock(counted2,countrate)
        write(*,*) 'Elapsed time on ang/spec likelihood calc (s): ',
     &   real(counted2 - counted1)/real(countrate)
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 5. Combined likelihood
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Put together the number, angular and spectral likelihoods
      lnlike = nLikelihood + specAngLikelihood


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 6. p value
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Calculate pvalue
      if (pValFromRef) then !compute p-value with reference to input log likelihood

        deltalnlike = max(0.d0, referenceLike-lnlike)
        pvalue = DGAMIC(dof*0.5d0,deltalnlike)/DGAMMA(dof*0.5d0)

      else                  !compute p-value with reference to background

        !p-value from Poissonian statistics
        if (.not. pvalBGPoisComputed(analysis)) call nulike_bglikeprecomp
        pvalue = nulike_pval(nEvents(analysis), theta_tot, theta_S, theoryError)
        pvalue = pvalue / BGpvalPoissonian(analysis)

      endif

      if (doProfiling) then
        call system_clock(counted1,countrate)
        write(*,*) 'Elapsed time on pval calc (s): ',
     &   real(counted1 - counted2)/real(countrate)
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 7. Data return
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Export various counts
      NBG_expected = theta_BG(analysis)
      Nsignal_predicted = theta_S
      Ntotal_observed = nEvents(analysis)

      end subroutine nulike_bounds

