***********************************************************************
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
***   effareafile       path to the file containing the IceCube
***                      effective area (or volume) and angular resolution.
***                      Ignored if eventfile indicates that the analysis
***                      uses the 2014 likelihood.
***   nchandistfile     path to the file containing distributions of
***                      the number of DOMs in the IceCube detector  
***                      triggered by neutrinos of different energies.
***                      Ignored if eventfile indicates that the analysis
***                      uses the 2014 likelihood.
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
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Mar 20, 2011
*** Modified: Mar 5, Jun 3, 6 2014
***********************************************************************


      subroutine nulike_init(analysis_name, eventfile, BGfile, 
     & effareafile, nchandistfile, phi_cut, theoryError, uselogNorm, 
     & BGLikePrecompute)

      implicit none
      include 'nulike.h'

      character (len=*) analysis_name, eventfile, nchandistfile, 
     & BGfile, effareafile
      character (len=30) instring, instring2
      integer IFAIL, i, nBinsRunning, nnchan(max_nHistograms), nchan
      integer BGfirst, BGsecond, BGcurrent, altind(2), nulike_amap
      integer nBinsBGE
      real*8 phi_cut, theoryError
      logical BGLikePrecompute, uselogNorm
      external nulike_amap

      !Roll credits.
      if (.not. nulike_init_called) then 
        write(*,*) 
        write(*,*) 'I like, you like...'
        write(*,*) '**********************************************************'
        write(*,*) '*                      nulike 1.0                        *'
        write(*,*) '*               Pat Scott, Chris Savage                  *'
        write(*,*) '*         JCAP (2012) 11:057, arXiv:1207.0810)           *'
        write(*,*) '*         JCAP (2014) xx:xxx, arXiv:141y.yyyy)           *'
        write(*,*) '**********************************************************'
        nulike_init_called = .true.      
      endif

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
 
      !Set maximum opening angle from solar centre to consider
      phi_max_deg(analysis) = phi_cut
      phi_max_rad(analysis) = phi_cut*pi/180.d0

      !Set percentage theoretical error
      theoryErr(analysis) = theoryError

    
      !Open event file, determine the total number of events and likelihood version

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

      read(lun, fmt=*) instring, instring2

      if (instring .ne. '[v]') then
        write(*,*) 'Bad format in IC event file ',eventfile,'.'
        write(*,*) 'First non-comment line begins with: ',instring
        write(*,*) 'Quitting...'
        stop
      endif 

      read(instring2, fmt='(I4)') likelihood_version(analysis)

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

      read(instring2, fmt=*) exp_time(analysis)

      do
        read(lun, fmt=*, IOSTAT=IFAIL, END=30) instring
        read(lun, fmt=*, IOSTAT=IFAIL, END=30) instring, nEvents_in_file(analysis)
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
      
      if (nEvents_in_file(analysis) .gt. max_nEvents) then
        write(*,*) 'Event file '//trim(eventfile)//' contains more bins than'
        write(*,*) 'nulike has been configured to handle.'
        write(*,*) 'Increase max_nEffAreaBins in nulike.h and' 
        write(*,*) 'recompile.'
        stop
      endif

      !Read in the actual details of all events.
      call nulike_eventinit(eventfile, nEvents_in_file(analysis))


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
            if (BGcurrent .eq. angular) nBinsBGAng(analysis) = nBinsRunning + 1
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
        nBinsBGAng(analysis) = nBinsBGAng(analysis) + 1
      else if (BGfirst.ne.nchannels .and. BGsecond.ne.nchannels) then
        nBinsBGE = nBinsRunning + 1
      endif

      if (nBinsBGAng(analysis) .gt. max_nBinsBGAng .or.
     &    nBinsBGE .gt. max_nBinsBGE) then
        write(*,*) 'Background file contains more data than'
        write(*,*) 'nulike has been configured to handle.'
        write(*,*) 'Increase max_nBinsBGE or max_nBinsBGAng'
        write(*,*) 'in nulike.h and recompile.'
        stop
      endif


      ! Switch according to likelihood version.
      select case (likelihood_version(analysis))

      ! 2012 likelihood, as per arXiv:1207.0810 (load the effective area, PSF and energy dispersion.)
      case (2012)

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
      
        nBinsEA(analysis) = 0
        do
          do i = 1,5
            read(lun, fmt=*, IOSTAT=IFAIL, END=10), instring
          enddo
          read(lun, fmt=*, IOSTAT=IFAIL, END=10) instring, nBinsEA(analysis)
          if (IFAIL .ne. 0) then
           write(*,*) 'Bad format in effective area/volume file ',effareafile,'.'
           write(*,*) 'Quitting...'
           stop
          endif
        enddo

10      close(lun)
      
        nBinsEA(analysis) = nBinsEA(analysis) + 1

        if (nBinsEA(analysis) .gt. max_nBinsEA) then
          write(*,*) 'Effective area/volume file contains more bins than'
          write(*,*) 'nulike has been configured to handle.'
          write(*,*) 'Increase max_nEffAreaBins in nulike.h and' 
          write(*,*) 'recompile.'
          stop
        endif

        !Read in the actual effective area and PSF data. 
        call nulike_eainit(effareafile,nBinsEA(analysis))


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
        
        nchan_min(analysis) = dble(ridiculousNumberOfChannels)
        nchan_max(analysis) = 0.d0
        nHistograms(analysis) = 0

        do

          read(lun, fmt='(A20)', IOSTAT=IFAIL, END=40) instring

          if (instring(1:1) .eq. 'h') then
            read(instring, fmt=*, IOSTAT=IFAIL) instring2, 
     &       nnchan(nHistograms(analysis)+1), nchan
            if (nchan .lt. nint(nchan_min(analysis))) nchan_min(analysis) = dble(nchan)
            if (nchan .gt. nint(nchan_max(analysis))) nchan_max(analysis) = dble(nchan)
          else if (instring(1:1) .eq. 'H') then
            read(instring, fmt=*, IOSTAT=IFAIL) instring2, nHistograms(analysis)
            nnchan(nHistograms(analysis)) = nnchan(nHistograms(analysis)) + 1
          endif     

          if (IFAIL .ne. 0) then
            write(*,*) 'Bad format in energy dispersion histogram file'
            write(*,*) nchandistfile,'.'
            write(*,*) 'Quitting...'
            stop
          endif

        enddo

40      close(lun)

        nnchan_total(analysis) = nint(nchan_max(analysis) - nchan_min(analysis)) + 1

        if (nnchan_total(analysis) .gt. max_nnchan) then
          write(*,*) 'IC nchan histogram file contains more'
          write(*,*) 'nchan values than nulike has'
          write(*,*) 'been configured to handle.  Increase' 
          write(*,*) 'max_nnchan in nulike.h and recompile.'
          stop
        endif

        nHistograms(analysis) = nHistograms(analysis) + 1
        nnchan(nHistograms(analysis)) = nnchan(nHistograms(analysis)) + 1
      
        if (nHistograms(analysis) .gt. max_nHistograms) then
          write(*,*) 'IC nchan histogram file contains more histograms'
          write(*,*) 'than nulike has been configured to handle.'
          write(*,*) 'Increase max_nHistograms in nulike.h and' 
          write(*,*) 'recompile.'
          stop
        endif
      
        !Read in the actual background data
        call nulike_bginit(BGfile, nBinsBGAng(analysis), nBinsBGE, BGfirst, BGsecond, likelihood_version(analysis))

        !Read in the actual nchan response histograms and rearrange them into energy dispersion estimators
        call nulike_edispinit(nchandistfile, nHistograms(analysis), nnchan)


      !2014 likelihood, as per arXiv:141x.xxxx (load the precalculated effective area and partial likelihoods.)
      case (2014)!FIXME

        !Read in the actual background data
        call nulike_bginit(BGfile, nBinsBGAng(analysis), nBinsBGE, BGfirst, BGsecond, likelihood_version(analysis))

        !Set nchan_min and nchan_max
        !do it from vector of nchan abcissae for bg spectrum

        !Read in the partial likelihoods
        !call nulike_bgspecanginit()

      case default
        write(*,*) "Unrecognised likelihood version in nulike_init."
        write(*,*) "Quitting..."
        stop

      end select



      !Calculate the expected background count.
      call nulike_bgpredinit

      !Precompute the background p-value (confidence level) for the Poissonian likelihood if requested.  
      !This is used for calculation of the final p-value for each model if nulike_bounds is called with pvalFromRef = F.
      pvalBGPoisComputed(analysis) = .false.
      if (BGLikePrecompute) call nulike_bglikeprecomp

      write(*,*) "Initialisation of nulike analysis '"//trim(analysis_name)//"' complete."

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


