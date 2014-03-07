***********************************************************************
*** nulike_aeinit initialises the IceCube effective area and its errors,
*** as well as the angular resolution, and calculates superbin
*** boundaries if required.
***
*** input: filename       name of effective area file
***        nbins          number of energy bins in file
***        doSuperBinning indicate whether events will be sorted into
***                        superbins or just analysed as a single block
***        superBinWid    the maximum change in the error (as a fraction
***                        of the effective area) allowed within a
***                        single superbin.  The resultant superbinning
***                        will thus be such that the error on the
***                        effective area will be constant within each
***                        superbin to within 100%*superBinWid of the
***                        effective area in that bin.     
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: April 8, 2011
***********************************************************************

      subroutine nulike_eainit(filename,nbins,doSuperBinning,superBinWid)

      implicit none
      include 'nulike.h'

      character (len=*) filename
      character (len=1) instring
      integer nbins, i, prevStartIndex, IER
      real*8 totalSystematic(max_nBinsEA), superBinWid
      real*8 runningAverage, prevWidth, prevStartlogE, firstSystematic
      real*8 firstWidth, totalWidth, newWidth, diff, firstdiff, enddiff
      real*8 prevContribution, newContribution, working(2*nbins-2)
      logical doSuperBinning, clusterdebug
      parameter(clusterdebug = .false.)

      !Open effective area file for reading
      open(lun,file=filename, ACTION='READ')

      !Skip header
      instring = '#'
      do while (instring .eq. '#')
        read(lun, fmt='(A1)'), instring
      enddo

      !Read in effective neutrino and anti-neutrino areas, uncertainties on eff area and angular resolutions
      do i = 1,nbins
        read(lun, *) instring, 
     &   effArea_logE(1,i), effArea_logE(2,i)
        !Convert to log(GeV) scale 
        effArea_logE(1,i) = dlog10(effArea_logE(1,i))
        effArea_logE(2,i) = dlog10(effArea_logE(2,i))
        !Find log-weighted bin centres 
        effArea_logEcentres(i) = 0.5d0*(effArea_logE(1,i)
     &   +effArea_logE(2,i))
        read(lun, *) instring, 
     &   effArea_nu(i), effArea_nubar(i)
        read(lun, *) instring, 
     &   effArea_syserr(i), effArea_staterr(i)
        read(lun, *) instring, 
     &   effArea_AngRes(i)
        read(lun, fmt='(A1)') instring
        if (i .ne. nbins) read(lun, fmt='(A1)'), instring
      enddo

      close(lun)

      !Calculate the percentage total systematic error from the effective 
      !area estimation in each bin, as the quadrature sum of the true
      !systematic and the Monte Carlo statistical error.  (Both contributions
      !have the character of a pure systematic when considered in the 
      !context of using a fixed effective area to do an analysis of real events.)
      totalSystematic = dsqrt(effArea_syserr*effArea_syserr + 
     & effArea_staterr*effArea_staterr) * 0.01d0

      !A very rough hueuristic 1D clustering algorithm,
      !meant to identify an appropriate rebinning such that
      !each new bin has constant error to within some specified
      !width.
      nBinsEAError = 0
      runningAverage = 0.d0
      prevStartIndex = 1
      prevWidth = 0.d0
      prevStartlogE = effArea_logE(1,1)
      firstSystematic = totalSystematic(1)
      firstWidth = effArea_logE(2,1) - prevStartlogE

      do i = 1,nbins

        EAErr_max = max(totalSystematic(i),EAErr_max)

        if (doSuperBinning) then

          diff = abs(totalSystematic(i)-firstSystematic)
          if (i .eq. prevStartIndex + 1) firstdiff = diff

          if (diff .gt. superBinWid) then
            !Errors have changed more than superBinWidth over sliding window; 
            !time to make a new bin.
            nBinsEAError = nBinsEAError + 1
            enddiff = abs(totalSystematic(i)-totalSystematic(i-1))
            EAlogE_inEAErrBins(1,nBinsEAError) = prevStartlogE
            if (enddiff.ge.firstdiff .or. i.eq.prevStartIndex+1) then
              !The difference at the end of the window is greater than
              !at the start; finish the last bin at the previous step
              !and start the new bin at this step.
              if (clusterdebug) write(*,*) 'End drop at bin',i
              prevStartIndex = i
              prevStartlogE = effArea_logE(1,i)
              prevWidth = 0.d0      
              EAlogE_inEAErrBins(2,nBinsEAError) = prevStartlogE
              EAErr_inEAErrBins(nBinsEAError) = runningAverage
              runningAverage = 0.d0
              firstWidth = effArea_logE(2,i) - effArea_logE(1,i)
              firstSystematic = totalSystematic(i)
            else
              !The difference at the start of the window is greater than
              !at the end; chop the first step off the current bin and
              !make it a dedicated bin; continue adding steps to the current
              !bin.
              if (clusterdebug) write(*,*) 'Back Drop at bin',i
              EAErr_inEAErrBins(nBinsEAError) = 
     &                             totalSystematic(prevStartIndex)
              totalWidth = effArea_logE(2,i-1) - prevStartlogE
              runningAverage = (runningAverage*totalWidth -
     &         firstSystematic*firstWidth)/(totalWidth-firstWidth)
              prevwidth = prevwidth - firstwidth
              prevStartIndex = prevStartIndex+1
              firstdiff = abs(totalSystematic(prevStartIndex+1)-
     &                       totalSystematic(prevStartIndex))
              firstWidth = effArea_logE(2,prevStartIndex)
     &                      - effArea_logE(1,prevStartIndex)
              firstSystematic = totalSystematic(prevStartIndex)
              prevStartlogE = effArea_logE(1,prevStartIndex)
              EAlogE_inEAErrBins(2,nBinsEAError) = prevStartlogE
            endif
          endif

          totalWidth = effArea_logE(2,i) - prevStartlogE
          newWidth = effArea_logE(2,i) - effArea_logE(1,i)
          prevContribution = runningAverage * prevWidth/totalWidth
          newContribution = totalSystematic(i) * newWidth/totalWidth

          runningAverage = prevContribution + newContribution
          prevWidth = prevWidth + newWidth

        endif

      enddo

      nBinsEAError = nBinsEAError + 1
      EAlogE_inEAErrBins(1,nBinsEAError) = prevStartlogE
      EAlogE_inEAErrBins(2,nBinsEAError) = effArea_logE(2,nbins)
      if (doSuperBinning) then
        EAErr_inEAErrBins(nBinsEAError) = runningAverage
      else
        EAErr_inEAErrBins(nBinsEAError) = EAErr_max
      endif
      
      !Find which bin gives the most conservative estimate of the
      !error.
      maxEAErrIndex = 1
      do i = 1, nBinsEAError
        if (EAErr_inEAErrBins(i) .gt. 
     &   EAErr_inEAErrBins(maxEAErrIndex)) then
          maxEAErrIndex = i
        endif
      enddo
      

      !Now need to init the interpolators in effective area and angular resolution.

      !Set up interpolation in neutrino effective area
      call TSPSI(nbins,effArea_logEcentres,effArea_nu,
     & 2,0,.false.,.false.,2*nbins-2,working,effArea_nuderivs,
     & effArea_nusigma,IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_eainit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up neutrino eff area.'
        stop
      endif

      !Set up interpolation in anti-neutrino effective area
      call TSPSI(nbins,effArea_logEcentres,effArea_nubar,
     & 2,0,.false.,.false.,2*nbins-2,working,effArea_nubarderivs,
     & effArea_nubarsigma,IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_eainit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up anti-nu eff area.'
        stop
      endif

      !Set up interpolation in angular resolution
      call TSPSI(nbins,effArea_logEcentres,effArea_AngRes,
     & 2,0,.false.,.false.,2*nbins-2,working,effArea_AngResderivs,
     & effArea_AngRessigma,IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_eainit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up angular resolution.'
        stop
      endif

      end subroutine nulike_eainit
