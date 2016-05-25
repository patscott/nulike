***********************************************************************
*** nulike_utils contains various simple utility routines.
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Jun 15 2014
***********************************************************************

      subroutine nulike_credits

      implicit none
      include 'nuversion.h'
      include 'nucommon.h'

      if (credits_rolled) return
      write(*,*)
      write(*,*) 'I like, you like...'
      write(*,*) '**********************************************************'
      write(*,*) '*                    nulike '//version//'            *'
      write(*,*) '*               Pat Scott, Chris Savage                  *'
      write(*,*) '*         JCAP (2012) 11:057, arXiv:1207.0810            *'
      write(*,*) '*         JCAP (2016) 04:022, arXiv:1601.00653           *'
      write(*,*) '**********************************************************'
      credits_rolled = .true.

      end subroutine nulike_credits


      character(len=6) function evnmshrfmt(i)

      implicit none
      integer i

      if (i .ge. 100000) then
        evnmshrfmt = '(I6)'
      elseif (i .ge. 10000) then
        evnmshrfmt = '(I5)'
      elseif (i .ge. 1000) then
        evnmshrfmt = '(I4)'
      elseif (i .ge. 100) then
        evnmshrfmt = '(I3)'
      elseif (i .ge. 10) then
        evnmshrfmt = '(I2)'
      else
        evnmshrfmt = '(I1)'
      endif

      end function
