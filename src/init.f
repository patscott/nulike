***********************************************************************
*** nulike_init initialises IceCube data from included or user-supplied
*** files.  The only things actually done here rather than in subroutines
*** are to determine the amount of data in each file, and whether DarkSUSY
*** (being a stubborn adherent to strict F77 and its static allocation
*** rules) has been correctly configured to store that much data.
***
*** input: eventfile    path to the file containing IceCube event data
***                      and total exposure time.
***   nchandistfile     path to the file containing distributions of
***                      the number of DOMs in the IceCube detector  
***                      triggered by neutrinos of different energies.
***   BGfile            path to the file containing the distribution
***                      of arrival directions and number of hit DOMs
***                      for the observed background, as well as the
***                      total number of background events.
***   effareafile       path to the file containing the IceCube
***                      effective area (or volume) and angular resolution.
***   phi_cut	        cutoff angle; likelihoods and p-values will be 
***			 based only on events with reconstructed 
***			 directions within this angle of the solar centre.
***			 [degrees]
***   theoryError       theoretical error to incorporate into likelihood
***                      and p-value calculations (given as a fractional
***                      relative error).
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
*** Hidden options (hard-coded to be turned off):
***   doSuperBinning    T => Perform the likelihood calculations
***                          in multiple 'superbins', where unbinned
***                          likelihoods are calculated within a small
***                          number of energy superbins, so that the error
***                          on the effective area is allowed to very
***                          with energy.
***                     F => Perform the entire calculation using a
***                          constant systematic error on the effective
***                          area, given by the largest value of the 
***                          error at any energy.
***   superbinWid       Sets the target range over which the systematic
***                      error on the effective area is allowed to vary
***                      over the width of each superbin. A percentage 
***                      of the effective area, i.e. if = 0.01, the super-
***                      bins will be chosen such that the fractional error
***                      on the effective area varies by 0.01 over the
***                      width of the superbins.  Ignored if doSuperBinning
***                      = false.
***   pLawIndx          Spectral index of model spectrum for Bayesian
***                      unfolding, used for sorting of events into 
***                      superbins. Ignored if doSuperBinning = false.
*** Note that superbinning only applies to the likelihood calculation; the p 
*** value is always calculated over the entire sample, without any superbinning. 
***
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Mar 20, 2011
*** Modified: Mar 5, 2014
***********************************************************************


      subroutine nulike_init(eventfile, nchandistfile, BGfile, 
     & effareafile, phi_cut, theoryError, uselogNorm, BGLikePrecompute)

      implicit none

      include 'nulike.h'

      character (len=*) eventfile, nchandistfile, BGfile, 
     & effareafile
      character (len=30) instring, instring2
      integer IFAIL, i, nBinsRunning, nnchan(max_nHistograms), nchan
      integer BGfirst, BGsecond, BGcurrent, altind(2)
      real*8 superBinWid, pLawIndx, phi_cut
      real*8 theoryError
      logical doSuperBinning, BGLikePrecompute, uselogNorm

      !Hard-coded superbinning options
      parameter (doSuperBinning = .false.)
      parameter (superBinWid = 0.01d0)
      parameter (pLawIndx = 3.4d0)

      if (.not. nulike_init_called) nulike_init_called = .true.

      !Choose whether to have a Gaussian distribution for the assumed PDF of 
      !systematic errors on the effective area/volume or a log-normal distribution
      sysErrDist_logNorm = uselogNorm

      !Set maximum opening angle from solar centre to consider
      phi_max_deg = phi_cut
      phi_max_rad = phi_max_deg*pi/180.d0

      !Set percentage theoretical error
      theoryErr = theoryError

    
      !Open neutrino effective area/volume file, determine number of bins
      
      open(lun,file=effareafile,IOSTAT=IFAIL, ACTION='READ')
      if (IFAIL .ne. 0) then
        write(*,*) 'Error opening effective area/volume file. ',effareafile,'.'
        write(*,*) 'Quitting...'
        stop
      endif

      instring = '#'
      do while (instring .eq. '#')
        read(lun, fmt='(A1,I1)'), instring        
      enddo

      if (instring .ne. 'B') then
        write(*,*) 'Bad format in effective area/volume file ',effareafile,'.'
        write(*,*) 'First non-comment line begins with: ',instring
        write(*,*) 'Quitting...'
        stop
      endif 
      
      nBinsEA = 0
      do
        do i = 1,5
          read(lun, fmt=*, IOSTAT=IFAIL, END=10), instring
        enddo
        read(lun, fmt=*, IOSTAT=IFAIL, END=10) instring, nBinsEA
        if (IFAIL .ne. 0) then
         write(*,*) 'Bad format in effective area/volume file ',effareafile,'.'
         write(*,*) 'Quitting...'
         stop
        endif
      enddo

10    close(lun)
      
      nBinsEA = nBinsEA + 1

      if (nBinsEA .gt. max_nBinsEA) then
        write(*,*) 'Effective area/volume file contains more bins than'
        write(*,*) 'DarkSUSY has been configured to handle.'
        write(*,*) 'Increase max_nEffAreaBins in nulike.h and' 
        write(*,*) 'recompile.'
        stop
      endif


      !Open background file, determine numbers of bins for angular 
      !and nchan distributions, and which comes first

      open(lun,file=BGfile,IOSTAT=IFAIL, ACTION='READ')
      if (IFAIL .ne. 0) then
        write(*,*) 'Error opening IC background file. ',BGfile,'.'
        write(*,*) 'Quitting...'
        stop
      endif 

      instring = '#'
      do while (instring .ne. hstring(1) .and.
     & instring .ne. hstring(2) .and. 
     & instring .ne. hstring(3))
        read(lun, fmt=*), instring
      enddo
      
      do i = 1,3
        if (instring .eq. hstring(i)) BGfirst = i
      enddo
      BGcurrent = BGfirst

      do
        read(lun, fmt='(A20)', IOSTAT=IFAIL, END=20) instring2
       
        altind = mod(BGcurrent+[0,1],3) + 1

        do i = 1,2
          if (instring2 .eq. hstring(altind(i))) then
            if (BGcurrent .eq. BGfirst) BGsecond = altind(i)
            if (BGcurrent .eq. angular) nBinsBGAng = nBinsRunning + 1
            if (BGcurrent .eq. nchannels) nBinsBGE = nBinsRunning + 1
            BGcurrent = altind(i)
            read(lun, fmt=*, IOSTAT=IFAIL) instring2
          endif
        enddo
      
        read(instring2, fmt=*, IOSTAT=IFAIL) instring, nBinsRunning

        read(lun, fmt=*, IOSTAT=IFAIL, END=20), instring
        if (BGcurrent .ne. events)
     &   read(lun, fmt=*, IOSTAT=IFAIL, END=20), instring
   
        if (IFAIL .ne. 0) then
         write(*,*) 'Bad format in IC background file ',BGfile,'.'
         write(*,*) 'Quitting...'
         stop
        endif

      enddo

20    close(lun)
      
      if (BGfirst.ne.angular .and. BGsecond.ne.angular) then
        nBinsBGAng = nBinsBGAng + 1
      else if (BGfirst.ne.nchannels .and. BGsecond.ne.nchannels) then
        nBinsBGE = nBinsRunning + 1
      endif

      if (nBinsBGAng .gt. max_nBinsBGAng .or.
     &    nBinsBGE .gt. max_nBinsBGE) then
        write(*,*) 'Background file contains more data than'
        write(*,*) 'DarkSUSY has been configured to handle.'
        write(*,*) 'Increase max_nBinsBGE or max_nBinsBGAng'
        write(*,*) 'in nulike.h and recompile.'
        stop
      endif


      !Open event file, determine the total number of events

      open(lun,file=eventfile,IOSTAT=IFAIL, ACTION='READ')
      if (IFAIL .ne. 0) then
        write(*,*) 'Error opening IC event file. ',eventfile,'.'
        write(*,*) 'Quitting...'
        stop
      endif 

      instring = '#'
      do while (instring .ne. '###--Likelihood--')
        read(lun, fmt=*) instring
      enddo

      read(lun, fmt=*) instring, nulike_version

      if (instring(1:3) .ne. '[v]') then
        write(*,*) 'Bad format in IC event file ',eventfile,'.'
        write(*,*) 'First non-comment line begins with: ',instring
        write(*,*) 'Quitting...'
        stop
      endif 

      write(*,*) 
      write(*,*) '**********************************************************'
      write(*,*) '*                   nulike 1.0                           *'
      write(*,*) '*              Pat Scott, Chris Savage                   *'
      write(*,*) '*                 arXiv:1207.0810                        *'       
      write(*,*) '* Neutrino telescope likelihood version: ', nulike_version,' *' 
      write(*,*) '**********************************************************'

      do while (instring .ne. '###--Exposure--')
        read(lun, fmt=*) instring
      enddo

      read(lun, fmt=*) instring, instring2

      if (instring(1:3) .ne. '[t]') then
        write(*,*) 'Bad format in IC event file ',eventfile,'.'
        write(*,*) 'Second non-comment line begins with: ',instring
        write(*,*) 'Quitting...'
        stop
      endif 

      read(instring2, fmt=*) exp_time

      do
        read(lun, fmt=*, IOSTAT=IFAIL, END=30) instring
        read(lun, fmt=*, IOSTAT=IFAIL, END=30) instring, nEvents
        if (IFAIL .ne. 0) then
         write(*,*) 'Bad format in IC event file.',eventfile,'.'
         write(*,*) 'Quitting...'
         stop
        endif
        do i = 1,3
          read(lun, fmt=*, IOSTAT=IFAIL, END=30), instring
        enddo
      enddo

30    close(lun)
      
      if (nEvents .gt. max_nEvents) then
        write(*,*) 'IC event file contains more bins than'
        write(*,*) 'DarkSUSY has been configured to handle.'
        write(*,*) 'Increase max_nEffAreaBins in nulike.h and' 
        write(*,*) 'recompile.'
        stop
      endif

      !Open file of nchan response histograms (energy dispersions), determine how many histograms
      !and how many bins in each histogram.

      open(lun,file=nchandistfile,IOSTAT=IFAIL, ACTION='READ')
      if (IFAIL .ne. 0) then
        write(*,*) 'Error opening IC nchan distribution'
        write(*,*) ' file. ',nchandistfile,'. Quitting...'
        stop
      endif

      instring = '#'
      do while (instring .eq. '#')
        read(lun, fmt='(A1)'), instring
      enddo

      if (instring .ne. 'H') then
        write(*,*) 'Bad format in IC nchan histogram file'
        write(*,*) nchandistfile,'.'
        write(*,*) 'First non-comment line begins with: ',instring
        write(*,*) 'Quitting...'
        stop
      endif 
        
      nchan_min = ridiculousNumberOfChannels
      nchan_max = 0
      nHistograms = 0

      do

        read(lun, fmt='(A20)', IOSTAT=IFAIL, END=40) instring

        if (instring(1:1) .eq. 'h') then
          read(instring, fmt=*, IOSTAT=IFAIL) instring2, 
     &     nnchan(nHistograms+1), nchan
          if (nchan .lt. nchan_min) nchan_min = nchan
          if (nchan .gt. nchan_max) nchan_max = nchan
        else if (instring(1:1) .eq. 'H') then
          read(instring, fmt=*, IOSTAT=IFAIL) instring2, nHistograms
          nnchan(nHistograms) = nnchan(nHistograms) + 1
        endif     

        if (IFAIL .ne. 0) then
        write(*,*) 'Bad format in IC nchan histogram file'
        write(*,*) nchandistfile,'.'
        write(*,*) 'Quitting...'
        stop
          stop
        endif

      enddo

40    close(lun)

      nnchan_total = nchan_max - nchan_min + 1

      if (nnchan_total .gt. max_nnchan) then
        write(*,*) 'IC nchan histogram file contains more'
        write(*,*) 'nchan values than DarkSUSY has'
        write(*,*) 'been configured to handle.  Increase' 
        write(*,*) 'max_nnchan in nulike.h and recompile.'
        stop
      endif

      nHistograms = nHistograms + 1
      nnchan(nHistograms) = nnchan(nHistograms) + 1
      
      if (nHistograms .gt. max_nHistograms) then
        write(*,*) 'IC nchan histogram file contains more histograms'
        write(*,*) 'than DarkSUSY has been configured to handle.'
        write(*,*) 'Increase max_nHistograms in nulike.h and' 
        write(*,*) 'recompile.'
        stop
      endif
      
      !Read in the actual effective area/volume data, and create the effective area/volume superbins 
      call nulike_eainit(effareafile,nBinsEA,doSuperBinning,superBinWid)

      !Read in the actual background data
      call nulike_bginit(BGfile, nBinsBGAng, nBinsBGE, BGfirst, BGsecond)

      !Read in the actual nchan response histograms and rearrange them into energy dispersion estimators
      call nulike_edispinit(nchandistfile, nHistograms, nnchan, pLawIndx)

      !Read in the actual detail of all events, and sort them into the superbins determined from the errors on the effective area/volume
      call nulike_eventinit(eventfile, nEvents)

      !Calculate the expected background counts in each superbin
      call nulike_bgpredinit

      !Precompute the background p-value (confidence level) for the Poissonian likelihood if requested.  This is used for calculation of 
      !the final p-value for each model if nulike_bounds is called with pvalFromRef = F.
      pvalBGPoisComputed = .false.
      if (BGLikePrecompute) call nulike_bglikeprecomp

      write(*,*) 'Initialisation of nulike complete.'

      end subroutine nulike_init


C**********************************************************************
C Initialization of common block elements must appear in a block data
C routine
C**********************************************************************
      block data nulike_FlagBlock
      logical nulike_init_called     
      common /nulike_flags/ nulike_init_called
      data nulike_init_called/.false./
      end


