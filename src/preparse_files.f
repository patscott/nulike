***********************************************************************
*** These routines pre-parse various data files and return basic info
*** about them.
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: June 15 2014
***********************************************************************


      !Open event file, determine the total number of events and likelihood version
      subroutine nulike_preparse_eventfile(eventfile, nevents, exp_time, like)

      implicit none
      include 'nucommon.h'
      character (len=*) eventfile
      character (len=30) instring, instring2
      integer nevents, like, IFAIL, i
      real*8 exp_time

      open(lun, file=eventfile, IOSTAT=IFAIL, ACTION='READ', STATUS='OLD')
      if (IFAIL .ne. 0) then
        write(*,*) 'Error opening IC event file: ',trim(eventfile),'.'
        write(*,*) 'Quitting...'
        stop
      endif

      instring = '#'
      do while (instring .ne. '###--Likelihood--')
        read(lun, fmt=*) instring
      enddo

      read(lun, fmt=*) instring, instring2

      if (instring .ne. '[v]') then
        write(*,*) 'Bad format in neutrino telescope event file: ',trim(eventfile),'.'
        write(*,*) 'First non-comment line begins with: ',instring
        write(*,*) 'Quitting...'
        stop
      endif

      read(instring2, fmt='(I4)') like

      do while (instring .ne. '###--Exposure--')
        read(lun, fmt=*) instring
      enddo

      read(lun, fmt=*) instring, instring2

      if (instring(1:3) .ne. '[t]') then
        write(*,*) 'Bad format in neutrino telescope event file: ',trim(eventfile),'.'
        write(*,*) 'Second non-comment line begins with: ',instring
        write(*,*) 'Quitting...'
        stop
      endif

      read(instring2, fmt=*) exp_time

      nevents = 0
      do
        do i = 1,5
          read(lun, fmt=*, IOSTAT=IFAIL, END=30) instring
        enddo
        if (IFAIL .ne. 0) then
         write(*,*) 'Bad format in neutrino telescope event file: ',trim(eventfile),'.'
         write(*,*) 'Quitting...'
         stop
        endif
        nevents = nevents + 1
      enddo

30    close(lun)

      end subroutine nulike_preparse_eventfile


      !Open background file, determine numbers of bins for angular
      !and energy-estimator distributions, and which comes first
      subroutine nulike_preparse_bgfile(BGfile, nbins_ang, nbins_E, first, second)

      implicit none
      include 'nucommon.h'
      character (len=*) BGfile
      character (len=30) instring, instring2
      integer nbins_ang, nbins_E, nbins_run
      integer first, second, current, altind(2), IFAIL, i

      open(lun,file=BGfile,IOSTAT=IFAIL, ACTION='READ', STATUS='OLD')
      if (IFAIL .ne. 0) then
        write(*,*) 'Error opening neutrino telescope background file: ',trim(BGfile),'.'
        write(*,*) 'Quitting...'
        stop
      endif

      instring = '#'
      do while (instring .ne. hstring(1) .and.
     & instring .ne. hstring(2) .and.
     & instring .ne. hstring(3))
        read(lun, fmt=*) instring
      enddo

      do i = 1,3
        if (instring .eq. hstring(i)) first = i
      enddo
      current = first

      nbins_run = 0
      do
        read(lun, fmt='(A20)', IOSTAT=IFAIL, END=20) instring2

        altind = mod(current+[0,1],3) + 1

        do i = 1,2
          if (instring2 .eq. hstring(altind(i))) then
            if (current .eq. first) second = altind(i)
            if (current .eq. angular) nbins_ang = nbins_run
            if (current .eq. enrgyest) nbins_E = nbins_run
            nbins_run = 0
            current = altind(i)
            read(lun, fmt=*, IOSTAT=IFAIL) instring2
          endif
        enddo
        read(lun, fmt=*, IOSTAT=IFAIL, END=20) instring
        if (current .ne. events) read(lun, fmt=*, IOSTAT=IFAIL, END=20) instring

        if (IFAIL .ne. 0) then
         write(*,*) 'Bad format in neutrino telescope background file: ',BGfile,'.'
         write(*,*) 'Quitting...'
         stop
        endif

        nbins_run = nbins_run + 1

      enddo

20    close(lun)

      if (current .eq. angular) nbins_ang = nbins_run
      if (current .eq. enrgyest) nbins_E = nbins_run

      if (nbins_ang .gt. max_nBinsBGAng .or. nbins_E .gt. max_nBinsBGE) then
        write(*,*) 'Background file contains more data than'
        write(*,*) 'nulike has been configured to handle.'
        write(*,*) 'Increase max_nBinsBGE or max_nBinsBGAng'
        write(*,*) 'in nuconst.h and recompile.'
        stop
      endif


      end subroutine nulike_preparse_bgfile


      !Open neutrino effective area or volume file and determine number of bins
      subroutine nulike_preparse_effarea_or_volume(fname, nbins, rho, like)

      implicit none
      include 'nucommon.h'
      character (len=*) fname
      character (len=30) instring, instring2
      integer nbins, like, IFAIL, i
      real*8 rho

      open(lun, file=fname, IOSTAT=IFAIL, ACTION='READ', STATUS='OLD')
      if (IFAIL .ne. 0) then
        write(*,*) 'Error opening effective area/volume file: ',trim(fname),'.'
        write(*,*) 'Quitting...'
        stop
      endif

      instring = '#'
      do while (instring(1:1) .eq. '#' .and. instring .ne. '###--Density--')
        read(lun, fmt=*) instring
      enddo

      if (instring .eq. '###--Density--') then
        read(lun, fmt=*) instring, instring2
        if (instring .ne. 'rho') then
          write(*,*) 'Bad format in effective area/volume file: ',trim(fname),'.'
          write(*,*) 'First line in density section begins with: ',instring
          write(*,*) 'Quitting...'
          stop
        endif
        read(instring2, fmt=*) rho
        read(lun, fmt=*) instring
        read(lun, fmt=*) instring
      else
        if (like .eq. 2015) stop 'Density required in effective volume file for 2015 likelihood.'
        rho = 0.d0
      endif

      if (instring(1:1) .ne. 'B') then
        write(*,*) 'Bad format in effective area/volume file: ',trim(fname),'.'
        write(*,*) 'First line in response section begins with: ',instring
        write(*,*) 'Quitting...'
        stop
      endif

      nbins = 0
      do
        do i = 1,5
          read(lun, fmt=*, IOSTAT=IFAIL, END=10) instring
        enddo
        read(lun, fmt=*, IOSTAT=IFAIL, END=10) instring
        if (IFAIL .ne. 0) then
          write(*,*) 'Bad format in effective area/volume file: ',trim(fname),'.'
          write(*,*) 'Quitting...'
          stop
        endif
        nbins = nbins + 1
      enddo

10    close(lun)

      nbins = nbins + 1

      if (nbins .gt. max_nSensBins) then
        write(*,*) 'Effective area/volume file contains more bins than'
        write(*,*) 'nulike has been configured to handle.'
        write(*,*) 'Increase max_nSensBins in nuconst.h and'
        write(*,*) 'recompile.'
        stop
      endif

      end subroutine nulike_preparse_effarea_or_volume


      !Open file of energy estimator response histograms, determine how
      !many histograms and how many bins in each histogram.
      subroutine nulike_preparse_energy_dispersion(edispfile, nhist, ncol, ee_max, ee_min, like)

      implicit none
      include 'nucommon.h'
      character (len=*) edispfile
      character (len=30) instring, instring2
      integer nhist, ncol(max_nHistograms+2), like, IFAIL, dummy
      real*8 ee, ee_max, ee_min

      open(lun, file=edispfile, IOSTAT=IFAIL, ACTION='READ', STATUS='OLD')
      if (IFAIL .ne. 0) then
        write(*,*) 'Error opening distribution file of'
        write(*,*) 'neutrino telescope energy estimator:'
        write(*,*) trim(edispfile),'. Quitting...'
        stop
      endif

      instring = '#'
      do while (instring .eq. '#')
        read(lun, fmt='(A1)') instring
      enddo

      if (instring .ne. 'H') then
        write(*,*) 'Bad format in distribution file of'
        write(*,*) 'neutrino telescope energy estimator:'
        write(*,*) trim(edispfile),'.'
        write(*,*) 'First non-comment line begins with: ',instring
        write(*,*) 'Quitting...'
        stop
      endif

      nhist = 1
      ncol = 0
      ee_min = huge(ee_min)
      ee_max = -huge(ee_max)

      do

        read(lun, fmt='(A20)', IOSTAT=IFAIL, END=40) instring

        if (instring(1:1) .eq. 'h') then
          read(instring, fmt=*, IOSTAT=IFAIL) instring2, dummy, ee
          if (ee .lt. ee_min) ee_min = ee
          if (ee .gt. ee_max) ee_max = ee
          ncol(nhist) = ncol(nhist) + 1
        else if (instring(1:1) .eq. 'H') then
          read(instring, fmt=*, IOSTAT=IFAIL) instring2, dummy
          nhist = nhist + 1
        endif

        if (IFAIL .ne. 0) then
          write(*,*) 'Bad format in energy dispersion histogram file:'
          write(*,*) trim(edispfile),'.'
          write(*,*) 'Quitting...'
          stop
        endif

      enddo

40    close(lun)

      if (nhist .gt. max_nHistograms) then
        write(*,*) 'Neutrino telescope energy dispersion histogram'
        write(*,*) 'file contains more histograms than nulike has'
        write(*,*) 'been configured to handle. Increase '
        write(*,*) 'max_nHistograms in nuconst.h and recompile.'
        stop
      endif

      if (
     &     (maxval(ncol) .gt. max_ncols)
     &     .or.
     &     (like .eq. 2012 .and. nint(ee_max - ee_min) + 1 .gt. max_ncols)
     &   ) then
        write(*,*) 'Neutrino telescope energy dispersion histogram'
        write(*,*) 'file contains more energy estimator values'
        write(*,*) 'than nulike has been configured to handle.'
        write(*,*) 'Increase max_ncols in nuconst.h and recompile.'
        stop
      endif

      end subroutine nulike_preparse_energy_dispersion
