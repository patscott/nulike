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
***   phi_cut	        cutoff angle; likelihoods and p-values will be 
***			 based only on events with reconstructed 
***			 directions within this angle of the solar centre.
***			 [degrees]
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Jun 15 2014
***********************************************************************


      subroutine nulike_partials(eventfile, effvolfile, edispfile, 
     & partialfile, phi_cut) 

      implicit none
      include 'nucommon.h'
      include 'nuprep.h'

      character (len=*) eventfile, effvolfile, edispfile, partialfile
      real*8 phi_cut, ee_min, ee_max, exp_time, density
      integer like, nhist, ncols(max_nHistograms), nEvents_in_file, nEvents, nbins_effvol

      !Roll credits.
      call nulike_credits

      !Only compile one set of partial likelihoods at a time.
      nAnalyses = 1
      analysis = 1

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

      
      end subroutine nulike_partials
