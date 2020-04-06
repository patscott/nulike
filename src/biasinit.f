***********************************************************************
*** nulike_biasinit initialises the energy-dependent analysis bias
*** factor, using the precomputed unbiased effective area and the
*** input effective area file.
***
*** This routine is used only with the 2015 likelihood.
***
*** input:  filename  name effective area file, or 'no-bias' if a zero-
***         bias calculation is to be done.
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Jul 27 2015
***********************************************************************

      subroutine nulike_biasinit(filename)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      character (len=*) filename
      character (len=1) instring
      integer nbins, i, ptype, IER
      real*8 temp1(1), temp2(1), EA_unbiased(2), totalSystematic(max_nBiasBins + 1)
      real*8 biased_sens_syserr(max_nBiasBins + 1,max_analyses)
      real*8 biased_sens_staterr(max_nBiasBins + 1,max_analyses)
      real*8 working(2*max_nBiasBins - 2), dummy

      !If a no-bias analysis is requested, just abort.
      if (trim(filename) .eq. 'no-bias') return

      !Open neutrino effective area file and determine number of bins
      call nulike_preparse_effarea_or_volume(filename, nbins, dummy, 2012)

      !Save bin number for other routines.  Extra 'bin' is at low-E edge of first bin.
      nBiasBins(analysis) = nbins + 1

      !Open effective area/volume file for reading
      open(lun,file=filename, ACTION='READ')

      !Skip header
      instring = '#'
      do while (instring .eq. '#')
        read(lun, fmt='(A1)') instring
      enddo
      !Skip density block if it exists
      if (instring .eq. 'r') then
        read(lun, fmt='(A1)') instring
        read(lun, fmt='(A1)') instring
      endif

      !Read in effective neutrino and anti-neutrino areas, uncertainties on eff areas and angular
      !resolutions (if present), and use them to calculate the bias factors at each energy.
      do i = 2, nbins + 1

        !Read from the file
        read(lun, *) instring, bias_logE(1,i,analysis), bias_logE(2,i,analysis)
        read(lun, *) instring, bias_nu(i,analysis), bias_nubar(i,analysis)
        read(lun, *) instring, biased_sens_syserr(i,analysis), biased_sens_staterr(i,analysis)
        read(lun, *) instring, dummy
        if (instring(1:1) .ne. '#') read(lun, fmt='(A1)') instring
        if (i .ne. nbins + 1) read(lun, fmt='(A1)') instring

        !Convert to log(GeV) scale
        bias_logE(1,i,analysis) = dlog10(bias_logE(1,i,analysis))
        bias_logE(2,i,analysis) = dlog10(bias_logE(2,i,analysis))

        !Find log-weighted bin centres
        bias_logEcentres(i,analysis) = 0.5d0*(bias_logE(1,i,analysis)+bias_logE(2,i,analysis))

        !Call interpolator to get unbiased, non angularly-corrected effective areas for this energy
        temp1(1) = bias_logEcentres(i,analysis)
        do ptype = 1, 2
          call TSVAL1(nPrecompE(analysis)-start_index_noL(analysis)+1,
     &     precomp_log10E(start_index_noL(analysis):,analysis),
     &     precompEAnoL_weights(start_index_noL(analysis):,ptype,analysis),
     &     precompEAnoL_derivs(:,ptype,analysis),
     &     precompEAnoL_sigma(:,ptype,analysis),
     &     0,1,temp1,temp2,IER)
           EA_unbiased(ptype) = temp2(1)
        enddo

        !Rescale the biased effective area by the unbiased effective area to get the bias factor
        bias_nu(i,analysis) = merge(10.d0**(log10(bias_nu(i,analysis)) - EA_unbiased(1)), 0.d0, EA_unbiased(1) .gt. logZero)
        bias_nubar(i,analysis) = merge(10.d0**(log10(bias_nubar(i,analysis)) - EA_unbiased(2)), 0.d0, EA_unbiased(2) .gt. logZero)

      enddo

      close(lun)

      if (IER .lt. 0) then
        write(*,*) 'TSVAL1 error from unbiased no-L effective area'
        write(*,*) 'in nulike_biasinit, code:', IER
        stop
      endif

      !Make the bias factor go to zero at the lower edge of the first bin.
      bias_logE(1,1,analysis) = bias_logE(1,2,analysis)
      bias_logE(2,1,analysis) = bias_logE(1,2,analysis)
      bias_logEcentres(1,analysis) = bias_logE(1,2,analysis)
      bias_nu(1,analysis) = 0.d0
      bias_nubar(i,analysis) = 0.d0
      !Set the first entry to have the same errors as the first bin.
      biased_sens_syserr(1,analysis) = biased_sens_syserr(2,analysis)
      biased_sens_staterr(1,analysis) = biased_sens_staterr(2,analysis)
      !Set the entries above the uppermost bin to have the same errors as the edge bins.
      if (nbins .lt. max_nBiasBins) then
        biased_sens_syserr(nbins+2:,analysis) = biased_sens_syserr(nbins+1,analysis)
        biased_sens_staterr(nbins+2:,analysis) = biased_sens_staterr(nbins+1,analysis)
      endif

      !Calculate the percentage total systematic error from the effective
      !area estimation in each bin, as the quadrature sum of the true
      !systematic and the Monte Carlo statistical error.  (Both contributions
      !have the character of a pure systematic when considered in the
      !context of using a fixed effective area/volume to do an analysis of real events.)
      totalSystematic = dsqrt(biased_sens_syserr(:,analysis)**2 +
     & biased_sens_staterr(:,analysis)**2) * 0.01d0
      !Overwrite the value read in from the partial likelihood auxiliarly file, based on the
      !effective volume (i.e. prefer the error on the effective area error if using bias factors).
      EAErr(analysis) = maxval(totalSystematic)

      !Set up interpolation in the neutrino bias factor
      call TSPSI(nbins+1,bias_logEcentres(:,analysis),bias_nu(:,analysis),
     & 2,0,.false.,.false.,2*nbins,working,bias_nuderivs(:,analysis),
     & bias_nusigma(:,analysis),IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_biasinit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up nu/lepton effective area/vol.'
        stop
      endif

      !Set up interpolation in anti-neutrino bias factor
      call TSPSI(nbins+1,bias_logEcentres(:,analysis),bias_nubar(:,analysis),
     & 2,0,.false.,.false.,2*nbins,working,bias_nubarderivs(:,analysis),
     & bias_nubarsigma(:,analysis),IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_biasinit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up anti-nu eff area/vol.'
        stop
      endif

      end subroutine nulike_biasinit
