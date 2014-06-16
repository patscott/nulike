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
***   partialfile       path to the file to be created or overwritten
***                      with the results of the partial likelihood
***                      calculation.
***   nEnergies         number of neutrino energies to tabluate the
***                      partial likelihoods for.
***   logE_min,logE_max limits of the tabulation energies    
***   phi_cut	        cutoff angle; likelihoods and p-values will be 
***			 based only on events with reconstructed 
***			 directions within this angle of the solar centre.
***			 [degrees]
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Jun 15 2014
***********************************************************************


      subroutine nulike_partials(eventfile, effvolfile, edispfile, 
     & partialfile, nEnergies, logE_min, logE_max, phi_cut) 

      implicit none
      include 'nucommon.h'
      include 'nuprep.h'

      character (len=*) eventfile, effvolfile, edispfile, partialfile
      real*8 phi_cut, ee_min, ee_max, exp_time, density, logE_min 
      real*8 logE_max, working(2*max_nHistograms-2)
      real*8 partial_likes(max_nPrecompE,max_nEvents,2)
      integer like, ncols(max_nHistograms), nEvents_in_file 
      integer nEvents, nbins_effvol, nEnergies, i, IER

      !Roll credits.
      call nulike_credits
      write(*,*)
      write(*,*) 'Computing partial likelihoods and saving to '
      write(*,*) trim(partialfile)//'.'

      !Only compile one set of partial likelihoods at a time.
      nAnalyses = 1
      analysis = 1

      !Set the internal cut angle
      phi_max = phi_cut

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
    
      !Step through the events and compute partial likelihoods for each one.
      do eventnumshare = 1, nEvents

        write(*,*) '  Computing partial likelihood for event ',eventnumshare

        !Arrange the energy dispersion function for this event.
        do i = 1, nhist
          !If the measured value of the energy estimator is outside the range of this histogram, set the prob to zero.
          if (events_ee(eventnumshare) .lt. hist_ee_flip(1,i) .or. 
     &        events_ee(eventnumshare) .gt. hist_ee_flip(ncols(i),i) ) then
            hist_single_ee_prob(i) = 0.d0
          else !If the measured value of the energy estimator is in the range of this histogram, set the prob by interpolating.
            call TSVAL1(ncols(i),hist_ee_flip(:,i),hist_prob_flip(:,i),hist_derivs_flip(:,i),
     &       hist_sigma_flip(:,i),0,1,events_ee(eventnumshare),hist_single_ee_prob(i),IER)
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
          Eshare = logE_min + dble(i-1)/dble(nEnergies-1)*(logE_max - logE_min)
          !Iterate over CP eigenstates
          do ptypeshare = 1, 2
            !partial_likes(i,eventnumshare,ptype) = nulike_simpson()
          enddo
        enddo

      enddo


      !Save auxiliary info about the partial likelihoods
      open(lun, file=partialfile//'.aux', action='WRITE')
      write(lun, fmt='(2I8,3E16.5)') nEvents, nEnergies, phi_cut, logE_min, logE_max
      close(lun)

      !Save partial likelihoods all in one hit.
      !open(lun, file=partialfile, form='unformatted', action='WRITE', recl=nEvents*nEnergies*2*8)
      !write(lun, rec=1) partial_likes(1:nEnergies,1:nEvents,:)
      !close(lun)
      
      write(*,*)
      write(*,*) 'Done.'
      write(*,*)

      end subroutine nulike_partials
