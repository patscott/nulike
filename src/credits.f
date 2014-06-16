***********************************************************************
*** nulike_credits does just what you expect.
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Jun 15 2014
***********************************************************************

      subroutine nulike_credits

      implicit none
      include 'nucommon.h'

      if (credits_rolled) return
      write(*,*) 
      write(*,*) 'I like, you like...'
      write(*,*) '**********************************************************'
      write(*,*) '*                      nulike 1.0                        *'
      write(*,*) '*               Pat Scott, Chris Savage                  *'
      write(*,*) '*         JCAP (2012) 11:057, arXiv:1207.0810)           *'
      write(*,*) '*         JCAP (2014) xx:xxx, arXiv:141y.yyyy)           *'
      write(*,*) '**********************************************************'
      credits_rolled = .true.

      end subroutine nulike_credits
