! This program prepares the partial likelihood data files used by nulike
! for computing the 2014 likelihood.  It requires Fortran 2003 for 
! parsing the command-line arguments, but the rest of nulike does not.
!
! Author: Pat Scott patscott@phsyics.mcgill.ca
! Date: Jun 12 2014

      program nulike_prep

      implicit none
      include 'nucommon.h'
      include 'nuprep.h'

      character (len=300) files(5)
      real*8 phi_cut
      integer i, gotarg

      !Check that the correct number of arguments have been passed on the command line
      if (command_argument_count() .ne. 5) then
        write(*,*)
        write(*,*) 'Usage: nulike_prep f1 f2 f3 f4 angle'
        write(*,*)
        write(*,*) '   f1  path to the file containing IceCube event data '
        write(*,*) '       and total exposure time.'
        write(*,*)
        write(*,*) '   f2  path to the file containing the IceCube'
        write(*,*) '       effective volume and angular resolution.'
        write(*,*)
        write(*,*) '   f3  path to the file containing distributions of'
        write(*,*) '       the energy estimator (e.g. the number of DOMs'
        write(*,*) '       hit in the IceCube detector) for neutrinos of'
        write(*,*) '       different energies.'
        write(*,*)
        write(*,*) '   f4  path to the file to be created or overwritten'
        write(*,*) '       with the results of the partial likelihood'
        write(*,*) '       calculation.'
        write(*,*)
        write(*,*) 'angle  cutoff angle; likelihoods and p-values will be '
        write(*,*) '       based only on events with reconstructed '
        write(*,*) '       directions within this angle of the solar centre.'
        write(*,*) '       [degrees]'
        write(*,*)
        stop
      endif

      !Retrieve command-line arguments
      do i = 1, 5
        call get_command_argument(i,files(i),status=gotarg)
        if (gotarg .ne. 0) stop 'Error parsing command-line arg in nulike_prep.'
      enddo
      read(files(5),*) phi_cut

      !Compute the partial likelihoods and output them in partialfile
      call nulike_partials(trim(files(1)),trim(files(2)),trim(files(3)),
     & trim(files(4)),phi_cut)

      end program nulike_prep
