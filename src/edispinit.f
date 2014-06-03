***********************************************************************
*** nulike_edispinit initialises the IceCube energy dispersion function.
*** This entails reading in the nchan response distributions for neutrinos
*** with energies in certain bands.
***
*** input:  filename     name of file containig energy dispersion data
***         nbins_Ein    number of histograms corresponding to responses
***                       to neutrinos with energies in certain bands
***         nbins_nchan  array of size nbins_Ein, indicating how many
***                       values of nchan are included in each histogram
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: April 8, 2011
*** Modified: Jun 3, 2014
***********************************************************************

      subroutine nulike_edispinit(filename, nbins_Ein, 
     & nbins_nchan)

      implicit none
      include 'nulike.h'

      character (len=*) filename
      character (len=20) instring
      integer nbins_Ein, nbins_nchan(nbins_Ein), dummyint, i, j, k, l
      integer hist_nchan_temp(max_nHistograms, max_nnchan)
      integer Eindex(2,max_nBinsEA), IER
      real*8  hist_prob_temp(max_nHistograms, max_nnchan)
      real*8  localdiff, normFactor
      real*8  energydiff(nHistograms)
      real*8  bestRelProb, tempRelProb
      real*8  working(2*nHistograms-2)

      !Read in nchan response distribution for each incoming neutrino energy band
      open(lun,file=filename, ACTION='READ')

      !Skip over header
      instring = '#'
      do while (instring .eq. '#')
        read(lun, fmt='(A1)'), instring
      enddo

      !Read actual data
      do i = 1, nbins_Ein
        read(lun, *) instring, hist_logE(1,i), hist_logE(2,i)
        hist_logEcentres(i) = 0.5d0*(hist_logE(1,i)+hist_logE(2,i))
        do j = 1, nbins_nchan(i)
          read(lun, *), instring, dummyint,
     &     hist_nchan_temp(i,j), hist_prob_temp(i,j)
        enddo
        read(lun,*) instring
        if (i .ne. nbins_Ein) read(lun,*) instring
      enddo

      close(lun)

      !Arrange histograms so they all cover the same range in nchan
      hist_prob = 0.d0
      do k = 1, nnchan_total
        do i = 1, nbins_Ein
          hist_nchan(i,k) = k - 1 + nchan_min
          do j = 1, nbins_nchan(i)
            if (hist_nchan_temp(i,j) .eq. hist_nchan(i,k)) then
              hist_prob(i,k) = hist_prob_temp(i,j)
            endif
          enddo
        enddo
      enddo

      !Work out where the indexing of nchan values in energy dispersion lines
      !up with indexing of nchan values in observed background spectrum.
      nchan_hist2BGoffset = -1
      do k = 1, nnchan_total
        if (hist_nchan(1,k) .eq. BGnchandist_nchan(1)) then 
          nchan_hist2BGoffset = k-1
        endif
      enddo


      !Work out the most likely energy bin (in terms of the 'superbins' over which the 
      !systematic error on the efective area is deemed to vary) that events with
      !each nchan are expected to have originated from.  This calculation effectively
      !estimates the integral of P(E|nchan) over each such bin, i.e. the probability 
      !of an event originating from that bin given an observed value of nchan.  To 
      !do this requires an assumption about the source shape; here we just assume 
      !a power law spectral model and allow the user to set the index 1.2d0.  Note
      !that this assumption *only* enters in classification of events (and background)
      !into different effective area systematic error bins - it does not enter into 
      !the unbinned likelihood calculation within each bin.  These probabilities are
      !also saved for estimating what fraction of observed background events should
      !be expected in each superbin.

      !Set up the bounds that say which histograms are to be included in full in the
      !count for which bins in effective area. 
      i = 1
      do k = 1,2
        do while (hist_logE(k,i) .lt. EAlogE_inEAErrBins(k,1))
          i = i + 1
        enddo
        Eindex(k,1) = i - k + 1
      enddo

      !Set up the spectral weighting factors for each histogram
      do k = 1, nHistograms
        energyDiff(k) = dexp(hist_logE(2,k)*(1.d0-1.2d0)) - 
     &                  dexp(hist_logE(1,k)*(1.d0-1.2d0))
      enddo

      !For each nchan, work out the relative probability of each 
      !effective area error bin and identify the bin in which the
      !probability is highest.
      do i = 1, nnchan_total

        BestGuessBin(i) = 1
        bestRelProb = 0.d0
        normFactor = 0.d0

        !Add on the integral over the first partial histogram
        l = Eindex(1,1)-1
        if (l .eq. 0) then
          tempRelProb = 0.d0
        else
          localdiff = dexp(hist_logE(2,l)*(1.d0-1.2d0)) -
     &                dexp(EAlogE_inEAErrBins(1,1)*(1.d0-1.2d0))
          tempRelProb = hist_prob(l,i) * localdiff
        endif

        !Add on the integrals over all the internal full histograms
        if (Eindex(1,1) .le. Eindex(2,1)) then
          do k = Eindex(1,1), Eindex(2,1)
            tempRelProb = tempRelProb + hist_prob(k,i) * energyDiff(k)
          enddo
        endif

        !Add on the integral over the final partial histogram
        l = Eindex(2,1)+1
        if (l .le. nHistograms) then
          localdiff = EAlogE_inEAErrBins(2,1)**(1.d0-1.2d0) -
     &                hist_logE(1,l)**(1.d0-1.2d0)
          tempRelProb = tempRelProb + hist_prob(l,i) * localdiff
        endif

        !Save relative probabilities, update normalisation factor
        !(abs is just because the neglected constant factor in the absolute 
        !probability goes -ve when tempRelProb also goes -ve)
        relProb(i,1) = abs(tempRelProb)
        normFactor = normFactor + relProb(i,1)

        !If relative probability of current bin is the best so far, save it
        if (abs(tempRelProb) .gt. bestRelProb) then
          bestRelProb = relProb(i,1)
          bestGuessBin(i) = 1
        endif

        !Make relative probabilities unitary
        if (normFactor .gt. 0.d0) then
          relProb(i,1) = relProb(i,1) / normFactor
        else !For nchans with no events in modelled histograms, make 
             !each bin just as likely as the next.
          relProb(i,1) = 1.d0
        endif
        
      enddo


      !Set up interpolation in energy histograms for use as energy dispersion estimator

      do i = 1, nnchan_total

        do j = 1, nHistograms
          edisp_prob(j) = hist_prob(j,i)
        enddo

        call TSPSI(nHistograms,hist_logEcentres,edisp_prob,
     &   2,0,.false.,.false.,2*nHistograms-2,working,edisp_derivs,
     &   edisp_sigma,IER)
        if (IER .lt. 0) then
          write(*,*) 'Error in nulike_edispinit: TSPSI failed with error'
          write(*,*) 'code',IER, ' at i=',i
          stop
        endif

        do j = 1, nHistograms
          hist_derivs(j,i) = edisp_derivs(j)
          hist_sigma(j,i) = edisp_sigma(j)
        enddo

      enddo 
      !Indicate which nchan data are currently loaded for
      nchansaved = nnchan_total + nchan_min - 1


      end subroutine nulike_edispinit

