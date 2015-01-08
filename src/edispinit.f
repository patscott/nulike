***********************************************************************
*** nulike_edispinit initialises the telescope energy dispersion function.
*** This entails reading in the ee response distributions for neutrinos
*** with energies in certain bands.
***
*** input:  filename     name of file containig energy dispersion data
***         nhgms        number of histograms corresponding to responses
***                       to neutrinos with energies in certain bands
***         nbins_ee     array of size nhgms, indicating how many
***                       values of ee are included in each histogram
***         min_ee       smallest value of the energy estimator in the
***                       dispersion function data.
***         like         likelihood version
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: April 8, 2011
*** Modified: Jun 3, 7, 15 2014
***********************************************************************

      subroutine nulike_edispinit(filename, nhgms, nbins_ee,
     & min_ee, like)

      implicit none
      include 'nulike_internal.h'
      include 'nuprep.h'

      character (len=*) filename
      character (len=20) instring
      integer nhgms, nbins_ee(nhgms), dummyint, i, j, k
      integer like, IER
      real*8  hist_ee_temp(max_nHistograms, max_ncols)
      real*8  hist_prob_temp(max_nHistograms, max_ncols)
      real*8  working(2*nHistograms(analysis)-2), min_ee
      real*8  working2(2*max_ncols-2)

      !Read in ee response distribution for each incoming neutrino energy band
      open(lun,file=filename, ACTION='READ')

      !Skip over header
      instring = '#'
      do while (instring .eq. '#')
        read(lun, fmt='(A1)'), instring
      enddo

      !Read actual data
      do i = 1, nhgms
        read(lun, *) instring, hist_logE(1,i,analysis), hist_logE(2,i,analysis)
        hist_logEcentres(i,analysis) = 0.5d0*(hist_logE(1,i,analysis)+
     &   hist_logE(2,i,analysis))
        do j = 1, nbins_ee(i)
          read(lun, *), instring, dummyint,
     &     hist_ee_temp(i,j), hist_prob_temp(i,j)
        enddo
        read(lun,*) instring
        if (i .ne. nhgms) read(lun,*) instring
      enddo

      close(lun)


      !Switch according to likelihood version.
      select case (like)

      !2012 likelihood, as per arXiv:1207.0810
      case (2012)

        !Arrange histograms so they all cover the same range in nchan
        hist_prob(:,:,analysis) = 0.d0
        do k = 1, nnchan_total(analysis)
          do i = 1, nhgms
            hist_nchan(i,k,analysis) = k - 1 + nint(min_ee)
            do j = 1, nbins_ee(i)
              if (hist_ee_temp(i,j) .eq. hist_nchan(i,k,analysis)) then
                hist_prob(i,k,analysis) = hist_prob_temp(i,j)
              endif
            enddo
          enddo
        enddo

        !Work out where the indexing of nchan values in energy dispersion lines
        !up with indexing of nchan values in observed background spectrum.
        nchan_hist2BGoffset(analysis) = -1
        do k = 1, nnchan_total(analysis)
          if (hist_nchan(1,k,analysis) .eq. BGeedist_ee(1,analysis)) then 
            nchan_hist2BGoffset(analysis) = k-1
          endif
        enddo

        !Set up interpolation in each nchan across energy histograms for use as energy dispersion estimator
        do i = 1, nnchan_total(analysis)
 
          call TSPSI(nHistograms(analysis),hist_logEcentres(:,analysis),hist_prob(:,i,analysis),
     &     2,0,.false.,.false.,2*nHistograms(analysis)-2,working,hist_derivs(:,i,analysis),
     &     hist_sigma(:,i,analysis),IER)
          if (IER .lt. 0) then
            write(*,*) 'Error in nulike_edispinit: TSPSI failed with error'
            write(*,*) 'code ',IER, ' at i=',i,' (like = ',like,').'
            stop
          endif

        enddo 

      !2014 likelihood, as per arXiv:141x.xxxx 
      case (2014)

        hist_ee_flip = transpose(hist_ee_temp)
        hist_prob_flip = transpose(hist_prob_temp)
        hist_logEnergies = hist_logEcentres(:,analysis)


        !Set up interpolation within each energy histogram, for later seeding of event-specific energy dispersion estimator
        do i = 1, nhgms

          call TSPSI(nbins_ee(i),hist_ee_flip(:,i),hist_prob_flip(:,i),2,0,.false.,.false.,
     &     2*nbins_ee(i)-2,working2,hist_derivs_flip(:,i),hist_sigma_flip(:,i),IER)
          if (IER .lt. 0) then
            write(*,*) 'Error in nulike_edispinit: TSPSI failed with error'
            write(*,*) 'code ',IER, ' at i=',i,' (like = ',like,').'
            stop
          endif

        enddo

      case default
        write(*,*) "Unrecognised likelihood version in nulike_edispinit."
        write(*,*) "Quitting..."
        stop

      end select


      end subroutine nulike_edispinit

