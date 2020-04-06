***********************************************************************
*** nulike_partials computes model-independent partial angular-spectral
*** likelihoods using data for a given analysis, and writes them to
*** disk.  These files are required for computing the final 2015-type
*** likelihoods.
*** This routine is used only with the 2015 likelihood.
***
*** input:
***   eventfile         path to the file containing IceCube event data
***                      and total exposure time.
***   effvolfile        path to the file containing the IceCube
***                      effective volume and angular resolution.
***   edispfile         path to the file containing distributions of
***                      the energy estimator (e.g. the number of DOMs
***                      hit in the IceCube detector) for neutrinos of
***                      different energies.
***   partialfolder     path to the folder in which to save the files
***                      with the results of the partial likelihood
***                      calculation.  If the folder does not exist, it
***                      is created,
***   nEnergies         number of neutrino energies to tabluate the
***                      partial likelihoods for.
***   logE_min,logE_max limits of the tabulation energies
***   phi_cut           cutoff angle; likelihoods and p-values will be
***          based only on events with reconstructed
***          directions within this angle of the solar centre.
***          [degrees]
***   dsdxdy            differential cross-section function, with signature
***                     Input:  real*8      E      Neutrino energy in GeV
***                             real*8      x      Bjorken-x
***                             real*8      y      Bjorken-y
***                             integer     nu     Neutrino type
***                                                 1 = nu_e
***                                                 2 = nu_e-bar
***                                                 3 = nu_mu
***                                                 4 = nu_mu-bar
***                                                 5 = nu_tau
***                                                 6 = nu_tau-bar
***                             character*1 targ   Target
***                                                 'p' = proton
***                                                 'n' = neutron
***                             character*2 int    Interaction type
***                                                 'CC'=charged current
***                                                 'NC'=neutral current
***                     Output: real*8 Differential cross section in cm^2
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Jun 15 2014
***********************************************************************


      subroutine nulike_partials(eventfile, effvolfile, edispfile,
     & partialfolder, nEnergies, logE_min, logE_max, phi_cut, dsdxdy)

      use iso_c_binding, only: c_ptr
      use Precision_Model
      use CUI

      implicit none
      include 'nucommon.h'
      include 'nuprep.h'

      character (len=*) eventfile, effvolfile, edispfile, partialfolder
      character (len=100) instring
      character (len=6) eventstring, evnmshrfmt
      character (len=9) shapenames(2)
      real*8 phi_cut, ee_min, ee_max, exp_time, density, log10E
      real*8 logE_min, logE_max, working(2*max_nHistograms+2)
      real*8 eff_areas(nEnergies,2), partial_likes(nEnergies,2), dsdxdy
      real*8 nEvents2, nEnergies2, phi_cut2, logE_min2, logE_max2
      real*8 SAbsErr, SValue, SVertices(2,3)
      logical abserr(2), is_fishy, ultracautious
      logical this_result_needs_correcting(2)
      logical previous_result_needs_correcting(2)
      integer shapes(2), job_indices(4), i, j, k, l, m, eventnum
      integer like, ncols(max_nHistograms+2), nEvents_in_file
      integer nEvents, nbins_effvol, nEnergies, IER, SRgType
      interface
        function nulike_partials_handoff(NumFun,X) result(Value)
        integer, intent(in) :: NumFun
        double precision, intent(in) :: X(:)
        double precision :: Value(NumFun)
        end function nulike_partials_handoff
      end interface
      external dsdxdy, is_fishy
      parameter (ultracautious = .false.)

      !Roll credits.
      call nulike_credits

      !Only compile one set of partial likelihoods at a time.
      analysis = 1

      !Set the internal cut angle
      phi_max = phi_cut

      !Set the lepton/neutino generation.  This must be elevated to a datafile input to support e and/or tau neutrinos.
      leptypeshare = 2

      !Set the differential cross-section function pointer
      dsdxdy_ptr%f => dsdxdy

      !Set the limits of integration in x and y. The x=1 contribution is tiny and causes issues for CTEQ6 DIS PDFs.
      SVertices(:,1) = (/0.9999999999d0, 0.d0/)
      SVertices(:,2) = (/0.9999999999d0, 1.d0/)
      SVertices(:,3) = (/0.d0,           0.d0/)

      !Set the integration adaptive shape type
      SRgType = Simplex

      !Open event file, determine the total number of events, likelihood version
      call nulike_preparse_eventfile(eventfile, nEvents_in_file, exp_time, like)
      !Die if the event file is not meant for 2015-type likelihoods.
      if (like .ne. 2015) stop 'Error: event file is not for 2015 likelihood.'
      !Set nEvents to zero to indicate to eventinit that it has not been determined elsewhere.
      nEvents = 0
      !Read in the actual details of all events.
      call nulike_eventinit(eventfile, nEvents_in_file, nEvents, dcos(phi_cut*pi/180.d0), 2015)

      !Open neutrino effective volume file and determine number of bins
      call nulike_preparse_effarea_or_volume(effvolfile, nbins_effvol, density, 2015)
      !Read in the actual effective volume and PSF data.
      call nulike_sensinit(effvolfile,nbins_effvol,2015)

      !Open file of energy estimators, determine how many histograms and how many bins in each histogram.
      call nulike_preparse_energy_dispersion(edispfile, nhist, ncols, ee_min, ee_max, 2015)
      !Read in the actual energy estimator response histograms and rearrange them into energy dispersion estimators
      call nulike_edispinit(edispfile, nhist, ncols, ee_min, 2015)

      !Calculate neutron and proton number per m^3 in detector, scaled up by 10^5 for later unit conversion.
      numdens_n = 8.d5 * density / m_water
      numdens_p = 10.d5 * density / m_water

      !Try to open the files for saving the unbiased effective areas
      open(lun, file=partialfolder//'/unbiased_effective_area.dat', form='unformatted',
     & action='WRITE', status='NEW', err=20, recl=nEnergies*2*8)
      open(lun2, file=partialfolder//'/unbiased_effective_area_noL.dat', form='unformatted',
     & action='WRITE', status='NEW', err=20, recl=nEnergies*2*8)

      !That worked, so now we can at least make some sort of claim about what is to come...
      write(*,*)
      write(*,*) 'Computing partial likelihoods and unbiased '
      write(*,*) 'effective areas, and saving to '
      write(*,*) trim(partialfolder)//'.'
      write(*,*) 'Note that existing files will be retained; '
      write(*,*) 'Please manually delete any partial likelihood'
      write(*,*) 'files that need to be regenerated.'
      write(*,*)

      !Save auxiliary info about the partial likelihoods
      open(lun3, file=partialfolder//'/partlikes.aux', action='WRITE')
      write(lun3, fmt=*) '#This file provides auxiliary information about partial likelihoods.'
      write(lun3, fmt=*) '#It was generated automatically by nulike_partials.'
      write(lun3, fmt=*) '# #events #energies, phi_cut, log10E_min, log10E_max, fractional_err_EA'
      write(lun3, fmt='(A1,2I8,4E16.5)') ' ', nEvents, nEnergies, phi_cut, logE_min, logE_max, EAErr(analysis)
      close(lun3)

      !Compute and save the effective area
      write(*,*) 'Computing unbiased effective areas.'
      !Step through each energy
      do i = 1, nEnergies
        log10E = logE_min + dble(i-1)/dble(nEnergies-1)*(logE_max - logE_min)
        Eshare = 10.d0**log10E
        write(*,*) '    Computing unbiased effective areas for E = ',Eshare,' GeV'
        !If the neutrino already has less energy than the lowest-E lepton that can be detected,
        !we know the effective volume will always be zero and the efffective area is zero.
        if (log10E .lt. min_detectable_logE) then
          partial_likes(i,:) = 0.d0
          eff_areas(i,:) = 0.d0
        else !Otherwise, we might see some leptons, so iterate over CP eigenstates
          eventnumshare(1) = 0
          do ptypeshare = 1, 2
            IER = 0
            call CUBATR(2,nulike_partials_handoff,SVertices,SRgType,
     &       SValue,SAbsErr,IFAIL=IER,EpsAbs=1.d-150,EpsRel=eps_partials,MaxPts=25000000,Job=11)
            if (IER .ne. 0) then
              write(*,*) 'Error raised by CUBATR in nulike_partials: ', IER
              stop
            endif
            call CUBATR()
            partial_likes(i,ptypeshare) = SValue
          enddo
          eventnumshare(1) = -1
          do ptypeshare = 1, 2
            IER = 0
            call CUBATR(2,nulike_partials_handoff,SVertices,SRgType,
     &       SValue,SAbsErr,IFAIL=IER,EpsAbs=1.d-150,EpsRel=eps_partials,MaxPts=25000000,Job=11)
            if (IER .ne. 0) then
              write(*,*) 'Error raised by CUBATR in nulike_partials: ', IER
              stop
            endif
            call CUBATR()
            eff_areas(i,ptypeshare) = SValue
          enddo
        endif
        write(*,*) '      Effective area (m^2), nu:                    ',partial_likes(i,1)
        write(*,*) '      Effective area (m^2), nubar:                 ',partial_likes(i,2)
        write(*,*) '      Effective area (m^2), nu, no angular cut:    ',eff_areas(i,1)
        write(*,*) '      Effective area (m^2), nubar, no angular cut: ',eff_areas(i,2)
      enddo
      !Save effective areas for this event.
      write(lun) partial_likes
      write(lun2) eff_areas
      close(lun)
      close(lun2)

      !Check that the already-saved auxiliary info matches
20    open(lun, file=partialfolder//'/partlikes.aux', err=60, action='READ')
      !Skip header
      instring = '#'
      do while (instring(1:1) .eq. '#' .or. instring(2:2) .eq. '#')
        read(lun, fmt='(A100)') instring
      enddo
      read(instring, fmt=*) nEvents2, nEnergies2, phi_cut2, logE_min2, logE_max2
      if (nEvents .ne. nEvents2 .or. nEnergies .ne. nEnergies2 .or. phi_cut .ne. phi_cut2
     & .or. logE_min .ne. logE_min2 .or. logE_max .ne. logE_max2) then
        write(*,*) 'Current nEvents, nEnergies, phi_cut, logE_min, logE_max:'
        write(*,*) nEvents, nEnergies, phi_cut, logE_min, logE_max
        write(*,*) 'Saved nEvents, nEnergies, phi_cut, logE_min, logE_max:'
        write(*,*) nEvents2, nEnergies2, phi_cut2, logE_min2, logE_max2
        write(*,*) 'Check your parameters, or delete '//partialfolder//' to start over.'
        stop
      endif

      !Step through the events and compute partial likelihoods for each one
      write(*,*) 'Computing partial likelihoods for ',nEvents,' events.'
      do eventnum = 1, nEvents
        eventnumshare(1) = eventnum

        write(eventstring,fmt=evnmshrfmt(eventnumshare(1))) eventnumshare(1)
        open(lun, file=partialfolder//'/partlike_event'//trim(eventstring)//'.dat', form='unformatted',
     &   action='WRITE', status='NEW', err=40, recl=nEnergies*2*8)

        write(*,*) '  Computing partial likelihood for event ',eventnumshare(1)

        !Arrange the energy dispersion function for this event.
        do i = 1, nhist
          !If the measured value of the energy estimator is outside the range of this histogram, set the prob to zero.
          if (events_ee(eventnumshare(1),analysis) .lt. hist_ee_flip(1,i) .or.
     &        events_ee(eventnumshare(1),analysis) .gt. hist_ee_flip(ncols(i),i) ) then
            hist_single_ee_prob(i) = 0.d0
          else !If the measured value of the energy estimator is in the range of this histogram, set the prob by interpolating.
            call TSVAL1(ncols(i),hist_ee_flip(:,i),hist_prob_flip(:,i),hist_derivs_flip(:,i),
     &       hist_sigma_flip(:,i),0,1,events_ee(eventnumshare(1),analysis),hist_single_ee_prob(i),IER)
            if (IER .lt. 0) then
              write(*,*) 'TSVAL1 error from energy dispersion'
              write(*,*) 'in nulike_partials, code: ', IER, 'i=',i
              stop
            endif
          endif
        enddo
        !Initialise the interpolator in energy of the energy dispersion for this event.
        call TSPSI(nhist,hist_logEnergies,hist_single_ee_prob,2,0,.false.,.false.,
     &   2*nhist-2,working,hist_single_ee_derivs,hist_single_ee_sigma,IER)
        if (IER .lt. 0) then
          write(*,*) 'TSPSI error from energy dispersion'
          write(*,*) 'in nulike_partials, code: ', IER
          stop
        endif

        !Step through each energy
        do i = 1, nEnergies

          log10E = logE_min + dble(i-1)/dble(nEnergies-1)*(logE_max - logE_min)
          Eshare = 10.d0**log10E
          write(*,*) '    Computing partial likelihoods for E = ',Eshare,' GeV'
          !If the neutrino already has less energy than the lowest-E lepton that can be detected,
          !we know the effective volume will always be zero and the partial likelihoods are zero.

          if (log10E .lt. min_detectable_logE) then

            partial_likes(i,:) = 0.d0
            previous_result_needs_correcting(:) = .false.

          else !Otherwise, we might see some leptons, so iterate over CP eigenstates

            do ptypeshare = 1, 2

              this_result_needs_correcting(ptypeshare) = .false.

              IER = 0
              call CUBATR(2,nulike_partials_handoff,SVertices,Simplex,
     &         SValue,SAbsErr,IFAIL=IER,EpsAbs=effZero,EpsRel=eps_partials,MaxPts=2100000000,Job=11)
              if (IER .ne. 0) then
                write(*,*) 'Error raised by CUBATR in nulike_partials: ', IER
                stop
              endif
              call CUBATR()
              partial_likes(i,ptypeshare) = SValue

              if (i .gt. 1) then
                !Possibly repeat the calculation, depending on the comparison to the previous energy bin.
                if (is_fishy(partial_likes(i,ptypeshare), partial_likes(i-1,ptypeshare))) then
                  write(*,*) '      Result looks fishy.  Systematically retrying with all possible '
                  write(*,*) '      integrator settings until it looks right.'
                  !Not very efficient to set these here, but it makes things clearer when you need to come back and mess with them.
                  shapes = (/Simplex, HyperQuad/)
                  shapenames = (/"Simplex  ", "HyperQuad"/)
                  abserr = (/.true., .false./)
                  job_indices = (/12,11,2,1/)
                  do j = 1,2
                    do k = 1,2
                      do l = 1,4
                        do m = 0,5
                          call try_integration(partial_likes(i,ptypeshare),partial_likes(i-1,ptypeshare),m,
     &                                         job_indices(l),shapes(j),shapenames(j),abserr(k),SVertices)
                        enddo
                      enddo
                    enddo
                  enddo
                  if (is_fishy(partial_likes(i,ptypeshare), partial_likes(i-1,ptypeshare))) then
                    write(*,*) '        All integration options exhausted.  Largest result: ',
     &               partial_likes(i,ptypeshare)
                    if (ultracautious .or. previous_result_needs_correcting(ptypeshare)) then
                      !Bail if result still looks suspicious and it can't be fixed.
                      if (previous_result_needs_correcting(ptypeshare)) then
                        write(*,*) '        Previous result was bogus too.'
                      else
                        write(*,*) '          This still does not look trustworthy.'
                      endif
                      write(*,*) '            Sorry, you will need to try adjusting the integrator in'
                      write(*,*) '            partials.f yourself, or decreasing the tolerance'
                      write(*,*) '            eps_partials in include/nuprep.h.'
                      stop 'Quitting.'
                    else
                      partial_likes(i,ptypeshare) = partial_likes(i-1,ptypeshare)
                      if (i .eq. nEnergies) then
                        write(*,*) '        Adopting result for previous energy: ', partial_likes(i,ptypeshare)
                      else
                        write(*,*) '        This value will be filled in by linear interpolation between the'
                        write(*,*) '        logs of the previous and next results.'
                        this_result_needs_correcting(ptypeshare) = .true.
                      endif
                    endif
                  endif
                endif
              endif

            enddo

          endif

          if (any(previous_result_needs_correcting)) then
            partial_likes(i-1,1) = 10.d0**(0.5d0*(log10(partial_likes(i-2,1)) + log10(partial_likes(i,1))))
            partial_likes(i-1,2) = 10.d0**(0.5d0*(log10(partial_likes(i-2,2)) + log10(partial_likes(i,2))))
            write(*,*) '      Deferred partial likelihood, nu:    ',partial_likes(i-1,1)
            write(*,*) '      Deferred partial likelihood, nubar: ',partial_likes(i-1,2)
          endif

          if (any(this_result_needs_correcting)) then
            write(*,*) '      Partial likelihood, nu:    deferred'
            write(*,*) '      Partial likelihood, nubar: deferred'
          else
            write(*,*) '      Partial likelihood, nu:    ',partial_likes(i,1)
            write(*,*) '      Partial likelihood, nubar: ',partial_likes(i,2)
          endif

          previous_result_needs_correcting = this_result_needs_correcting

        enddo

        !Save partial likelihoods for this event.
        write(lun) partial_likes
        close(lun)

40    enddo

      write(*,*)
      write(*,*) 'Done.'
      write(*,*)

      return

      !Try to work out if auxiliary file is just missing or if the whole directory is not there.
60    open(lun, file=partialfolder//'/.direxists', action='WRITE', status='REPLACE', err=80)
      write(*,*) 'Missing file:'
      write(*,*) partialfolder//'/partlikes.aux.'
      stop ' Nulike_partials cannot continue.'
80    write(*,*) 'Folder '//partialfolder//' is missing.'
      stop ' Please create it and try again.'

      end subroutine nulike_partials


      function nulike_partials_handoff(NumFun,X) result(Value)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'
      include 'nuprep.h'

      integer, intent(in) :: NumFun
      double precision, intent(in) :: X(:)
      double precision :: Value(NumFun)
      double precision :: nulike_partintegrand
      external nulike_partintegrand

      if (associated(dsdxdy_ptr%f)) then
        Value(1) = nulike_partintegrand(X(1), X(2), dsdxdy_ptr%f, eventnumshare(1),
     &   Eshare, ptypeshare, leptypeshare)
      else
        stop('Cannot call nulike_partials_handoff without setting dsdxdy_ptr%f')
      endif

      end function nulike_partials_handoff


      !Determine whether a new partial likelihood looks fishy or not.
      logical function is_fishy(now, before)
      implicit none
      double precision, intent(in) :: now, before
      is_fishy = (now .lt. 1.d-40 .and. now .lt. 1.d-5*before)
      end function is_fishy


      !Try a partial likelihood integration
      subroutine try_integration(now, before, key, job, intshape, intshapename, abserr, SVertices)

      use CUI

      implicit none
      include 'nucommon.h'
      include 'nuprep.h'

      double precision, intent(inout) :: now
      double precision, intent(in) :: before, SVertices(2,3)
      integer, intent(in) :: key, job, intshape
      character(len=9), intent(in) :: intshapename
      logical, intent(in) :: abserr
      logical :: is_fishy
      integer :: IER
      double precision:: SValue, SAbsErr
      interface
        function nulike_partials_handoff(NumFun,X) result(Value)
        integer, intent(in) :: NumFun
        double precision, intent(in) :: X(:)
        double precision :: Value(NumFun)
        end function nulike_partials_handoff
      end interface
      external is_fishy

      if (is_fishy(now, before)) then
        write(*,'(A,I2,A,I2,A,L1)') "         Retrying with "//intshapename//": Key=",
     &   key,": Job=",job,": Absolute error=",abserr
        IER = 1
        if (abserr) then
          call CUBATR(2,nulike_partials_handoff,SVertices,intshape,SValue,
     &     SAbsErr,IFAIL=IER,EpsRel=eps_partials,EpsAbs=effZero,
     &     MaxPts=10000000,Key=key,Job=job)
        else
          call CUBATR(2,nulike_partials_handoff,SVertices,intshape,SValue,
     &     SAbsErr,IFAIL=IER,EpsRel=eps_partials,MaxPts=10000000,
     &     Key=key,Job=job)
        endif
        if (IER .eq. 0) now = max(now, SValue)
        if (IER .eq. 1) then
          write(*,*) '           Failed to converge.  Result:',SValue,' Relative error: ',SAbsErr/SValue
          if (SAbsErr*1.d-1/SValue .le. eps_partials) then
            write(*,*) '           Relative error is within a factor of 10 of request; accepting as possibly correct.'
            now = max(now, SValue)
          endif
        endif
        call CUBATR()
      endif

      end subroutine try_integration

