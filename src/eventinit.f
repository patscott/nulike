***********************************************************************
*** nulike_eventinit reads in and initialises the IceCube event data.
***
*** input:   filename     name of file containing all event data
***          totalevents  total number of events contained in file
***          cosphimin    cosine(cutoff angle above which events will 
***                       be ignored)
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: April 8, 2011
*** Modified: March 5, June 15 2014
***********************************************************************

      subroutine nulike_eventinit(filename, totalevents, cutevents, cosphimin, like)

      implicit none
      include 'nulike.h'

      character (len=*) filename
      character (len=20) instring, instring2
      integer totalevents, cutevents, savedevents, j, like
      real*8 cosphimin, cosphi, cosphierr, ee

      !Save the current number of events for later comparison if it is already set.
      if (like .eq. 2014) savedevents = cutevents

      !Open event file and read in events
      open(lun,file=filename, ACTION='READ')

      !Skip over header lines
      instring = '#'
      do while (instring(1:1) .eq. '#' .or. instring .eq. '[t]' .or. instring .eq. '[v]')
        read(lun, fmt='(A3)'), instring
      enddo

      !Read in events
      cutevents = 0
      do j = 1, totalevents
        read(lun, *) instring, ee
        read(lun, *) instring, cosphi
        read(lun, *) instring, instring2, cosphierr
        read(lun, *) instring
        if (j .ne. totalevents) read(lun, *) instring
        !read in only events that have phi < phi_cut
        if (cosphi .gt. cosphimin) then
          cutevents = cutevents + 1
          select case (like)
          case (2012)
            events_nchan(cutevents,analysis) = nint(ee)
          case (2014)
            events_ee(cutevents,analysis) = ee
          case default
            stop 'Unrecognised likelihood type in nulike_eventinit.'
          end select
          if (cosphiErr .eq. 0.d0) stop 'Error in event file: an event with phi_err=0!'
          events_cosphi(cutevents,analysis) = cosphi
          events_cosphiErr(cutevents,analysis) = cosphiErr
        endif
      enddo

      if (like .eq. 2014 .and. savedevents .gt. 0 .and. cutevents .ne. savedevents)
     & stop 'Number of events in partial likelihood file and event file do not match!! Death.'

      close(lun)

      end subroutine nulike_eventinit

