! This program prepares the partial likelihood data files used by nulike
! for computing the 2014 likelihood.  It requires Fortran 2003 for 
! parsing the command-line arguments, but the rest of nulike does not.
! It requires nusigma for the neutrino-lepton interaction cross-section.
!
! Author: Pat Scott p.scott@imperial.ac.uk
! Date: Jun 12 2014

      program nulike_prep

      implicit none
      include 'nucommon.h'

      character (len=nulike_clen) files(8)
      real*8 phi_cut, lEmin, lEmax, Nudsdxdy
      integer i, gotarg, nE
      external Nudsdxdy

      !Check that the correct number of arguments have been passed on the command line
      if (command_argument_count() .ne. 8) then
        write(*,*)
        write(*,*) 'Usage: nulike_prep f1 f2 f3 f4 nE lEmin lEmax angle'
        write(*,*)
        write(*,*) '   f1  path to the file containing event data '
        write(*,*) '       and total exposure time.'
        write(*,*)
        write(*,*) '   f2  path to the file containing the detector'
        write(*,*) '       effective volume and angular resolution.'
        write(*,*)
        write(*,*) '   f3  path to the file containing distributions of'
        write(*,*) '       the energy estimator (e.g. the number of DOMs'
        write(*,*) '       hit in the IceCube detector) for neutrinos of'
        write(*,*) '       different energies.'
        write(*,*)
        write(*,*) '   f4  path to the folder in which to save the files' 
        write(*,*) '       with the results of the partial likelihood'
        write(*,*) '       calculation.  If the folder does not exist, the'
        write(*,*) '       user must first create it.'
        write(*,*) 
        write(*,*) '  nE   number of neutrino energies to tabluate the'
        write(*,*) '       partial likelihoods for.'
        write(*,*)
        write(*,*) 'lEmin  log10(lower neutrino energy in GeV for tabulation)' 
        write(*,*)
        write(*,*) 'lEmax  log10(upper neutrino energy in GeV for tabulation)'        
        write(*,*)
        write(*,*) 'angle  cutoff angle; likelihoods and p-values will be '
        write(*,*) '       based only on events with reconstructed '
        write(*,*) '       directions within this angle of the solar centre.'
        write(*,*) '       [degrees]'
        write(*,*)
        stop
      endif

      !Retrieve command-line arguments
      do i = 1, 8
        call get_command_argument(i,files(i),status=gotarg)
        if (gotarg .ne. 0) stop 'Error parsing command-line arg in nulike_prep.'
      enddo
      read(files(5),*) nE
      read(files(6),*) lEmin
      read(files(7),*) lEmax
      read(files(8),*) phi_cut

      !Initialise nusigma
      call nusetup

      !Compute the partial likelihoods and output them in partialfile
      call nulike_partials(trim(files(1)),trim(files(2)),trim(files(3)),
     & trim(files(4)),nE,lEmin,lEmax,phi_cut, Nudsdxdy)

      end program nulike_prep
