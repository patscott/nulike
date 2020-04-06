***********************************************************************
*** nulike_specanginit initialises the precomputed unbiased effective
*** areas and partial angular-spectral likelihoods by reading them in
*** from the specified directory.
*** This routine is used only with the 2015 likelihood.
***
*** input:  dirname      name of directory containig partial likelihoods.
***         nevents      number of events for which partial likelihoods
***                       are to be read in.
***         nenergies    number of energies for which partial likelihoods
***                       are to be read in for each event.
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Jun 15 2014
***********************************************************************

      subroutine nulike_specanginit(dirname, n_events, phi_cut, n_energies)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      character (len=*) dirname
      character (len=6) eventstring, evnmshrfmt
      character (len=100) instring
      integer n_events, n_energies, i, IER
      real*8 phi_cut, logE_min, logE_max
      real*8 working(2*max_nPrecompE-2)

      interface
        function nulike_read_weights(local_lun, filename, n_energies)
          implicit none
          real*8 nulike_read_weights(n_energies,2)
          character (len=*) filename
          integer local_lun, n_energies
        end function nulike_read_weights
      end interface

      !Read in the auxiliary info.
      open(lun, file=trim(dirname)//'/partlikes.aux', action='READ')
      !Skip header
      instring = '#'
      do while (instring(1:1) .eq. '#' .or. instring(2:2) .eq. '#')
        read(lun, fmt='(A100)') instring
      enddo
      read(instring, fmt=*) n_events, n_energies, phi_cut, logE_min, logE_max, EAErr(analysis)
      close(lun)

      !Use the tabulation bounds in energy to reconstruct the energies
      forall(i=1:n_energies) precomp_log10E(i,analysis) =
     & logE_min + dble(i-1)/dble(n_energies-1) * (logE_max - logE_min)

      !Open and read in from the unbiased effective area files.
      precompEA_weights(1:n_energies,:,analysis) =
     & nulike_read_weights(lun, trim(dirname)//'/unbiased_effective_area.dat', n_energies)
      precompEAnoL_weights(1:n_energies,:,analysis) =
     & nulike_read_weights(lun, trim(dirname)//'/unbiased_effective_area_noL.dat', n_energies)

      !Determine the index at which the unbiased effective area actually starts
      i = 1
      do
        if (precompEA_weights(i,1,analysis) .gt. 0.99d0*logZero .or.
     &      precompEA_weights(i,2,analysis) .gt. 0.99d0*logZero) then
          start_index(analysis) = i
          exit
        endif
        i = i + 1
      enddo

      !Set up interpolation in unbiased neutrino effective area
      call TSPSI(n_energies-start_index(analysis)+1,
     & precomp_log10E(start_index(analysis):,analysis),
     & precompEA_weights(start_index(analysis):,1,analysis),
     & 2,0,.false.,.false.,2*n_energies-2,working,
     & precompEA_derivs(:,1,analysis),
     & precompEA_sigma(:,1,analysis),IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_specanginit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up unbiased neutrino effective area.'
        stop
      endif

      !Set up interpolation in unbiased anti-neutrino effective area
      call TSPSI(n_energies-start_index(analysis)+1,
     & precomp_log10E(start_index(analysis):,analysis),
     & precompEA_weights(start_index(analysis):,2,analysis),
     & 2,0,.false.,.false.,2*n_energies-2,working,
     & precompEA_derivs(:,2,analysis),
     & precompEA_sigma(:,2,analysis),IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_specanginit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up unbiased nubar effective area.'
        stop
      endif

      !Determine the index at which the unbiased effective area without L actually starts
      i = 1
      do
        if (precompEAnoL_weights(i,1,analysis) .gt. 0.99d0*logZero .or.
     &      precompEAnoL_weights(i,2,analysis) .gt. 0.99d0*logZero) then
          start_index_noL(analysis) = i
          exit
        endif
        i = i + 1
      enddo

      !Set up interpolation in unbiased neutrino effective area without L
      call TSPSI(n_energies-start_index_noL(analysis)+1,
     & precomp_log10E(start_index_noL(analysis):,analysis),
     & precompEAnoL_weights(start_index_noL(analysis):,1,analysis),
     & 2,0,.false.,.false.,2*n_energies-2,working,
     & precompEAnoL_derivs(:,1,analysis),
     & precompEAnoL_sigma(:,1,analysis),IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_specanginit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up unbiased neutrino effective '
        write(*,*) 'area without angular correction.'
        stop
      endif

      !Set up interpolation in unbiased anti-neutrino effective area without L
      call TSPSI(n_energies-start_index_noL(analysis)+1,
     & precomp_log10E(start_index_noL(analysis):,analysis),
     & precompEAnoL_weights(start_index_noL(analysis):,2,analysis),
     & 2,0,.false.,.false.,2*n_energies-2,working,
     & precompEAnoL_derivs(:,2,analysis),
     & precompEAnoL_sigma(:,2,analysis),IER)
      if (IER .lt. 0) then
        write(*,*) 'Error in nulike_specanginit: TSPSI failed with error'
        write(*,*) 'code',IER,' in setting up unbiased nubar effective '
        write(*,*) 'area without angular correction.'
        stop
      endif

      !Loop over events
      do i = 1, n_events

        !Open and read in from the event file
        write(eventstring,fmt=evnmshrfmt(i)) i
        precomp_weights(1:n_energies,i,:,analysis) =
     &   nulike_read_weights(lun, trim(dirname)//'/partlike_event'//trim(eventstring)//'.dat', n_energies)

        !Prime the interpolator for the neutrino partial likelihood of this event
        call TSPSI(n_energies,precomp_log10E(:,analysis),precomp_weights(:,i,1,analysis),
     &   2,0,.false.,.false.,2*n_energies-2,working,precomp_derivs(:,i,1,analysis),
     &   precomp_sigma(:,i,1,analysis),IER)
        if (IER .lt. 0) then
          write(*,*) 'Error in nulike_specanginit: TSPSI failed with error'
          write(*,*) 'code',IER,' in setting up neutrino partial likelihood.'
          stop
        endif

        !Prime the interpolator for the anti-neutrino partial likelihood of this event
        call TSPSI(n_energies,precomp_log10E(:,analysis),precomp_weights(:,i,2,analysis),
     &   2,0,.false.,.false.,2*n_energies-2,working,precomp_derivs(:,i,2,analysis),
     &   precomp_sigma(:,i,2,analysis),IER)
        if (IER .lt. 0) then
          write(*,*) 'Error in nulike_specanginit: TSPSI failed with error'
          write(*,*) 'code',IER,' in setting up nubar partial likelihood.'
          stop
        endif

      enddo

      end subroutine nulike_specanginit


      !Does the reading-in of the weights from a binary file created by nulike_partials
      function nulike_read_weights(local_lun, filename, n_energies)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'
      real*8 nulike_read_weights(n_energies,2)
      character (len=*) filename
      integer local_lun, n_energies

      open(local_lun, file=filename, form='unformatted',
     & action='READ', status='OLD', recl=n_energies*2*8)
      read(local_lun) nulike_read_weights
      close(local_lun)
      where(nulike_read_weights .le. effZero) nulike_read_weights = effZero
      nulike_read_weights = log10(nulike_read_weights)

      end function nulike_read_weights
