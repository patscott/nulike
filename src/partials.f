***********************************************************************
*** nulike_partials computes model-independent partial angular-spectral
*** likelihoods using data for a given analysis, and writes them to 
*** disk.  These files are required for computing the final 2014-type 
*** likelihoods.
*** This routine is used only with the 2014 likelihood.
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
***   phi_cut	        cutoff angle; likelihoods and p-values will be 
***			 based only on events with reconstructed 
***			 directions within this angle of the solar centre.
***			 [degrees]
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
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Jun 15 2014
***********************************************************************


      subroutine nulike_partials(eventfile, effvolfile, edispfile, 
     & partialfolder, nEnergies, logE_min, logE_max, phi_cut, dsdxdy) 

      implicit none
      include 'nucommon.h'
      include 'nuprep.h'

      character (len=*) eventfile, effvolfile, edispfile, partialfolder
      character (len=100) instring
      character (len=6) eventstring, evnmshrfmt
      real*8 phi_cut, ee_min, ee_max, exp_time, density, logE_min 
      real*8 logE_max, dsdxdy, working(2*max_nHistograms-2)
      real*8 nulike_simpson, nulike_partintegrand1
      real*8 partial_likes(nEnergies,2)
      real*8 nEvents2, nEnergies2, phi_cut2, logE_min2, logE_max2
      integer like, ncols(max_nHistograms), nEvents_in_file 
      integer nEvents, nbins_effvol, nEnergies, i, IER
      external dsdxdy, nulike_partintegrand1

      !Roll credits.
      call nulike_credits

      !Only compile one set of partial likelihoods at a time.
      analysis = 1

      !Set the internal cut angle
      phi_max = phi_cut

      !Set the lepton/neutino generation.  This must be elevated to a datafile input to support e and/or tau neutrinos.
      leptypeshare = 2

      !Open event file, determine the total number of events, likelihood version
      call nulike_preparse_eventfile(eventfile, nEvents_in_file, exp_time, like)
      !Die if the event file is not meant for 2014-type likelihoods.
      if (like .ne. 2014) stop 'Error: event file is not for 2014 likelihood.'
      !Set nEvents to zero to indicate to eventinit that it has not been determined elsewhere.
      nEvents = 0
      !Read in the actual details of all events.
      call nulike_eventinit(eventfile, nEvents_in_file, nEvents, dcos(phi_cut*pi/180.d0), 2014)

      !Open neutrino effective volume file and determine number of bins
      call nulike_preparse_effarea_or_volume(effvolfile, nbins_effvol, density, 2014)
      !Read in the actual effective volume and PSF data.
      call nulike_sensinit(effvolfile,nbins_effvol)

      !Open file of energy estimators, determine how many histograms and how many bins in each histogram.
      call nulike_preparse_energy_dispersion(edispfile, nhist, ncols, ee_min, ee_max, 2014)
      !Read in the actual energy estimator response histograms and rearrange them into energy dispersion estimators
      call nulike_edispinit(edispfile, nhist, ncols, ee_min, 2014)
    
      !Calculate neutron and proton number per m^3 in detector, scaled up by 10^5 for later unit conversion.
      numdens_n = 8.d5 * density / m_water
      numdens_p = 10.d5 * density / m_water

      !Try to open the file for saving the effective area
      open(lun, file=partialfolder//'/effective_area.dat', form='unformatted', 
     & action='WRITE', status='NEW', err=20, recl=nEnergies*2*8)

      !That worked, so now we can at least make some sort of claim about what is to come...
      write(*,*)
      write(*,*) 'Computing partial likelihoods and effective '
      write(*,*) 'area, and saving to '
      write(*,*) trim(partialfolder)//'.'
      write(*,*) 'Note that existing files will be retained; '
      write(*,*) 'Please manually delete any partial likelihood'
      write(*,*) 'files that need to be regenerated.'
      write(*,*)

      !Save auxiliary info about the partial likelihoods
      open(lun2, file=partialfolder//'/partlikes.aux', action='WRITE')
      write(lun2, fmt=*) '#This file provides auxiliary information about partial likelihoods.'
      write(lun2, fmt=*) '#It was generated automatically by nulike_partials.'
      write(lun2, fmt=*) '# #events #energies, phi_cut, log10E_min, log10E_max'
      write(lun2, fmt='(A1,2I8,3E16.5)') ' ', nEvents, nEnergies, phi_cut, logE_min, logE_max
      close(lun2)

      !Compute and save the effective area
      write(*,*) 'Computing effective areas.'
      eventnumshare = 0
      !Step through each energy
      do i = 1, nEnergies
        log10Eshare = logE_min + dble(i-1)/dble(nEnergies-1)*(logE_max - logE_min)
        Eshare = 10.d0**log10Eshare
        write(*,*) '    Computing effective areas for E = ',Eshare,' GeV'
        !If the neutrino already has less energy than the lowest-E lepton that can be detected,
        !we know the effective volume will always be zero and the efffective area is zero.
        if (log10Eshare .lt. min_detectable_logE) then
          partial_likes(i,:) = 0.d0
        else !Otherwise, we might see some leptons, so iterate over CP eigenstates
          do ptypeshare = 1, 2
            partial_likes(i,ptypeshare) = nulike_simpson(nulike_partintegrand1,
     &       dsdxdy,0.d0,0.9999999999d0,eps_partials) !x=1 contribution is tiny and causes issues for CTEQ6 DIS PDFs
          enddo
        endif
        write(*,*) '      Effective area (m^2), nu:    ',partial_likes(i,1)
        write(*,*) '      Effective area (m^2), nubar: ',partial_likes(i,2)
      enddo
      !Save effective area for this event.
      write(lun) partial_likes
      close(lun) 

      !Check that the already-saved auxiliary info matches 
20    open(lun, file=partialfolder//'/partlikes.aux', err=60, action='READ')
      !Skip header
      instring = '#'
      do while (instring(1:1) .eq. '#' .or. instring(2:2) .eq. '#')
        read(lun, fmt='(A100)'), instring
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
      do eventnumshare = 1, nEvents

        write(eventstring,fmt=evnmshrfmt(eventnumshare)) eventnumshare
        open(lun, file=partialfolder//'/partlike_event'//trim(eventstring)//'.dat', form='unformatted',
     &   action='WRITE', status='NEW', err=40, recl=nEnergies*2*8)

        write(*,*) '  Computing partial likelihood for event ',eventnumshare

        !Arrange the energy dispersion function for this event.
        do i = 1, nhist
          !If the measured value of the energy estimator is outside the range of this histogram, set the prob to zero.
          if (events_ee(eventnumshare,analysis) .lt. hist_ee_flip(1,i) .or. 
     &        events_ee(eventnumshare,analysis) .gt. hist_ee_flip(ncols(i),i) ) then
            hist_single_ee_prob(i) = 0.d0
          else !If the measured value of the energy estimator is in the range of this histogram, set the prob by interpolating.
            call TSVAL1(ncols(i),hist_ee_flip(:,i),hist_prob_flip(:,i),hist_derivs_flip(:,i),
     &       hist_sigma_flip(:,i),0,1,events_ee(eventnumshare,analysis),hist_single_ee_prob(i),IER)
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
          log10Eshare = logE_min + dble(i-1)/dble(nEnergies-1)*(logE_max - logE_min)
          Eshare = 10.d0**log10Eshare
          write(*,*) '    Computing partial likelihoods for E = ',Eshare,' GeV'
          !If the neutrino already has less energy than the lowest-E lepton that can be detected,
          !we know the effective volume will always be zero and the partial likelihoods are zero.
          if (log10Eshare .lt. min_detectable_logE) then
            partial_likes(i,:) = 0.d0
          else !Otherwise, we might see some leptons, so iterate over CP eigenstates
            do ptypeshare = 1, 2
              partial_likes(i,ptypeshare) = nulike_simpson(nulike_partintegrand1,
     &         dsdxdy,0.d0,0.9999999999d0,eps_partials) !x=1 contribution is tiny and causes issues for CTEQ6 DIS PDFs
            enddo
          endif
          write(*,*) '      Partial likelihood, nu:    ',partial_likes(i,1)
          write(*,*) '      Partial likelihood, nubar: ',partial_likes(i,2)
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
