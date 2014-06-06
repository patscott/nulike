***********************************************************************
*** nulike_eventinit reads in and initialises the IceCube event data.
***
*** input:   filename     name of file containing all event data
***          totalevents  total number of events contained in file
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: April 8, 2011
*** Modified: March 5 2014
***********************************************************************

      subroutine nulike_eventinit(filename, totalevents)

      implicit none
      include 'nulike.h'

      character (len=*) filename
      character (len=20) instring, instring2
      integer totalevents, j, nchan
      real*8 cosphimin, cosphi, cosphierr


      !Define cutoff angle above which events will be ignored
      cosphimin = dcos(phi_max_rad(analysis))

      !Open event file and read in events
      open(lun,file=filename, ACTION='READ')

      !Skip over header lines
      instring = '#'
      do while (instring(1:1) .eq. '#' .or. instring .eq. '[t]' .or. instring .eq. '[v]')
        read(lun, fmt='(A3)'), instring
      enddo

      !Read in events
      nEvents(analysis) = 0
      do j = 1, totalevents
        read(lun, *) instring, nchan
        read(lun, *) instring, cosphi
        read(lun, *) instring, instring2, cosphierr
        read(lun, *) instring
        if (j .ne. nEvents_in_file(analysis)) read(lun, *) instring
        !read in only events that have phi < phi_cut
        if (cosphi .gt. cosphimin) then
          nEvents(analysis) = nEvents(analysis) + 1
          events_nchan(nEvents(analysis)) = nchan
          events_cosphi(nEvents(analysis)) = cosphi
          events_cosphiErr(nEvents(analysis)) = cosphiErr
        endif
      enddo

      close(lun)

      end subroutine nulike_eventinit

