*****************************************************************************
*** nulike_init initialises neutrino telescope data from user-supplied
*** files.  The only things actually done here rather than in subroutines
*** are to determine the amount of data in each file, and whether nulike
*** has been correctly configured to store that much data.
***
*** input:
***   analysis_name     a name by which to refer to this particular analysis.
***   eventfile         path to the file containing IceCube event data
***                      and total exposure time.
***   BGfile            path to the file containing the distribution
***                      of arrival directions and number of hit DOMs
***                      for the observed background, as well as the
***                      total number of background events.
***   effareafile       If the eventfile indicates that the likelihood is
***                      2012 type: path to the file containing the IceCube
***                                 effective area and (nu) angular resolution.
***                      2015 type: path to the file containing the IceCube
***                                 effective area, or just 'no-bias' if the
***                                 impacts of analysis-dependent biases are to
***                                 be ignored.
***   file4             If the eventfile indicates that the likelihood is
***                      2012 type: path to the file containing distributions
***                                 of the number of DOMs in the IceCube
***                                 detector triggered by neutrinos of
***                                 different energies.
***                      2015 type: path to the folder containing the partial
***                                 angular-spectral likelihoods, the
***                                 unbiased effective area, the cut angle
***                                 and all other parameters that they have
***                                 been computed with.
***   phi_cut           If the eventfile indicates that the likelihood is
***                      2012 type: cutoff angle; likelihoods and p-values
***                                 will be based only on events with
***                                 reconstructed directions within this
***                                 angle of the solar centre. [degrees]
***                      2015 type: ignored
***   uselogNorm        if false, assume a Gaussian distribution for the
***                      PDF of systematic errors (from the effective area/vol
***                      and theory errors).  If true, use a log-normal
***                      distribution instead.
***   BGLikePrecompute  If true, nulike_init precomputes the Possonian p-value
***                      for the background estimate to save time. This is
***                      later used to calculate the modified frequentist
***                      p-value for each model when nulike_bounds is called
***                      with pvalFromRef = F.
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Mar 20, 2011
*** Modified: Mar, Jun 2014
*****************************************************************************

      subroutine nulike_init(analysis_name, eventfile, BGfile,
     & effareafile, file4, phi_cut, uselogNorm, BGLikePrecompute)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      character (len=nulike_clen) analysis_name, eventfile, BGfile, effareafile, file4
      integer nnchan(max_nHistograms+2), nAnalyses
      integer BGfirst, BGsecond, nulike_amap, nbins, nhist, nEvents_available
      real*8 phi_cut, cosphimax, dummy
      logical BGLikePrecompute, uselogNorm
      external nulike_amap
      data nAnalyses /0/
      save nAnalyses

      !Roll credits.
      call nulike_credits

      !Set the flag indicating that initialisation has been done.
      if (.not. nulike_init_called) nulike_init_called = .true.

      !Make sure that this analysis is not already loaded.
      analysis = nulike_amap(analysis_name)
      if (analysis .ne. 0) then
        write(*,*) "Analysis '"//analysis_name//"' requested for load"
        write(*,*) 'in nulike_init is already loaded.'
        write(*,*) 'Returning without doing anything.'
        return
      endif

      !Register this analysis.
      nAnalyses = nAnalyses + 1
      analysis = nAnalyses
      analysis_name_array(analysis) = trim(analysis_name)

      !Choose whether to have a Gaussian distribution for the assumed PDF of
      !systematic errors on the effective area/volume or a log-normal distribution
      sysErrDist_logNorm(analysis) = uselogNorm

      !Open event file, determine the total number of events and likelihood version
      call nulike_preparse_eventfile(eventfile, nEvents_available, exp_time(analysis), likelihood_version(analysis))

      !Open background file, determine numbers of bins for angular
      !and nchan distributions, and which comes first
      call nulike_preparse_bgfile(BGfile, nBinsBGAng(analysis), nBinsBGE(analysis), BGfirst, BGsecond)


      !Switch according to likelihood version.
      select case (likelihood_version(analysis))

      !2012 likelihood, as per arXiv:1207.0810 (load the effective area, PSF and energy dispersion.)
      case (2012)

        !Make sure the user isn't confused about bias factors
        if (trim(effareafile) .eq. 'no-bias') stop 'nulike_init: \"no-bias\" only allowed with 2015 likelihood.'

        !Set maximum opening angle from solar centre to consider
        phi_max_deg(analysis) = phi_cut
        cosphimax = dcos(phi_cut*pi/180.d0)

        !Open neutrino effective area file and determine number of bins
        call nulike_preparse_effarea_or_volume(effareafile, nbins, dummy, 2012)
        !Read in the actual effective area and PSF data.
        call nulike_sensinit(effareafile, nbins, 2012)

        !Open file of nchan response histograms (energy dispersions), determine how many histograms
        !and how many bins in each histogram.
        call nulike_preparse_energy_dispersion(file4, nhist, nnchan,
     &   ee_min(analysis), ee_max(analysis), 2012)
        nnchan_total(analysis) = nint(ee_max(analysis) - ee_max(analysis)) + 1

        !Read in the actual background data
        call nulike_bginit(BGfile, nBinsBGAng(analysis), nBinsBGE(analysis), BGfirst, BGsecond, 2012)

        !Read in the actual nchan response histograms and rearrange them into energy dispersion estimators
        call nulike_edispinit(file4, nhist, nnchan, ee_min(analysis), 2012)

      !2015 likelihood, as per arXiv:1601.00653 (load the precalculated unbiased effective area and partial likelihoods.)
      case (2015)

        !Keep track of whether a no-bias estimate is being done or not.
        no_bias(analysis) = (trim(effareafile) .eq. 'no-bias')

        !Read in partial likelihoods, return number of events, cut angle and no. of energies tabulated over
        call nulike_specanginit(file4,nEvents(analysis),phi_max_deg(analysis),nPrecompE(analysis))
        cosphimax = dcos(phi_max_deg(analysis)*pi/180.d0)

        !Read in the simulated effective area and calculate the bias factors
        call nulike_biasinit(effareafile)

        !Read in the actual background data
        call nulike_bginit(BGfile, nBinsBGAng(analysis), nBinsBGE(analysis), BGfirst, BGsecond, 2015)

      case default
        write(*,*) "Unrecognised likelihood version in nulike_init: ", likelihood_version(analysis)
        write(*,*) "Quitting..."
        stop

      end select

      !Read in the actual details of all events.
      call nulike_eventinit(eventfile, nEvents_available, nEvents(analysis), cosphimax, likelihood_version(analysis))

      !Calculate the expected background count.
      call nulike_bgpredinit(cosphimax)

      !Precompute the background p-value (confidence level) for the Poissonian likelihood if requested.
      !This is used for calculation of the final p-value for each model if nulike_bounds is called with pvalFromRef = F.
      pvalBGPoisComputed(analysis) = .false.
      if (BGLikePrecompute) call nulike_bglikeprecomp

      write(*,*) "Initialisation of nulike analysis '"//trim(analysis_name)//"' complete."

      end subroutine nulike_init
