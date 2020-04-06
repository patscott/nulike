***********************************************************************
*** nulike_eventinit reads in and initialises the IceCube event data.
***
*** input:   filename     name of file containing all event data
***          totalevents  total number of events contained in file
***          cosphimin    cosine(cutoff angle above which events will
***                       be ignored)
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: April 8, 2011
*** Modified: March 5, June 15 2014
***********************************************************************

      subroutine nulike_eventinit(filename, totalevents, cutevents, cosphimin, like)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      character (len=*) filename
      character (len=20) instring, instring2
      integer totalevents, cutevents, savedevents, j, like
      real*8 cosphimin, cosphi, cosphierr, ee
      data savedevents /0/

      !Save the current number of events for later comparison if it is already set.
      if (like .eq. 2015) savedevents = cutevents

      !Open event file and read in events
      open(lun,file=filename, ACTION='READ')

      !Skip over header lines
      instring = '#'
      do while (instring(1:1) .eq. '#' .or. instring .eq. '[t]' .or. instring .eq. '[v]')
        read(lun, fmt='(A3)') instring
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
          if (cutevents .gt. max_nEvents) then
            write(*,*) 'Chosen angular cut includes more'
            write(*,*) 'events than nulike has been configured to '
            write(*,*) 'handle.  Increase max_nEvents in nuconst.h and'
            write(*,*) 'recompile.  You may need to reduce max_analyses'
            write(*,*) 'to do this.'
            stop
          endif
          select case (like)
          case (2012)
            events_nchan(cutevents,analysis) = nint(ee)
          case (2015)
            events_ee(cutevents,analysis) = ee
          case default
            stop 'Unrecognised likelihood type in nulike_eventinit.'
          end select
          if (cosphiErr .eq. 0.d0) stop 'Error in event file: an event with phi_err=0!'
          events_cosphi(cutevents,analysis) = cosphi
          events_cosphiErr(cutevents,analysis) = cosphiErr
        endif
      enddo

      if (like .eq. 2015 .and. savedevents .gt. 0 .and. cutevents .ne. savedevents)
     & stop 'Number of events in partial likelihood file and event file do not match!! Death.'

      close(lun)

      end subroutine nulike_eventinit

