***********************************************************************
*** nulike_edispinit initialises the IceCube energy dispersion function.
*** This entails reading in the nchan response distributions for neutrinos
*** with energies in certain bands.
***
*** input:  filename     name of file containig energy dispersion data
***         nbins_Ein    number of histograms corresponding to responses
***                       to neutrinos with energies in certain bands
***         nbins_nchan  array of size nbins_Ein, indicating how many
***                       values of nchan are included in each histogram
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: April 8, 2011
*** Modified: Jun 3, 7 2014
***********************************************************************

      subroutine nulike_edispinit(filename, nbins_Ein, 
     & nbins_nchan)

      implicit none
      include 'nulike.h'

      character (len=*) filename
      character (len=20) instring
      integer nbins_Ein, nbins_nchan(nbins_Ein), dummyint, i, j, k
      integer hist_nchan_temp(max_nHistograms, max_nnchan)
      integer IER
      real*8  hist_prob_temp(max_nHistograms, max_nnchan)
      real*8  working(2*nHistograms(analysis)-2)

      !Read in nchan response distribution for each incoming neutrino energy band
      open(lun,file=filename, ACTION='READ')

      !Skip over header
      instring = '#'
      do while (instring .eq. '#')
        read(lun, fmt='(A1)'), instring
      enddo

      !Read actual data
      do i = 1, nbins_Ein
        read(lun, *) instring, hist_logE(1,i,analysis), hist_logE(2,i,analysis)
        hist_logEcentres(i,analysis) = 0.5d0*(hist_logE(1,i,analysis)+
     &   hist_logE(2,i,analysis))
        do j = 1, nbins_nchan(i)
          read(lun, *), instring, dummyint,
     &     hist_nchan_temp(i,j), hist_prob_temp(i,j)
        enddo
        read(lun,*) instring
        if (i .ne. nbins_Ein) read(lun,*) instring
      enddo

      close(lun)

      !Arrange histograms so they all cover the same range in nchan
      hist_prob(:,:,analysis) = 0.d0
      do k = 1, nnchan_total(analysis)
        do i = 1, nbins_Ein
          hist_nchan(i,k,analysis) = k - 1 + nchan_min(analysis)
          do j = 1, nbins_nchan(i)
            if (hist_nchan_temp(i,j) .eq. hist_nchan(i,k,analysis)) then
              hist_prob(i,k,analysis) = hist_prob_temp(i,j)
            endif
          enddo
        enddo
      enddo

      !Work out where the indexing of nchan values in energy dispersion lines
      !up with indexing of nchan values in observed background spectrum.
      nchan_hist2BGoffset(analysis) = -1
      do k = 1, nnchan_total(analysis)
        if (hist_nchan(1,k,analysis) .eq. BGnchandist_nchan(1,analysis)) then 
          nchan_hist2BGoffset(analysis) = k-1
        endif
      enddo

      !Set up interpolation in energy histograms for use as energy dispersion estimator
      do i = 1, nnchan_total(analysis)

        call TSPSI(nHistograms(analysis),hist_logEcentres(:,analysis),hist_prob(:,i,analysis),
     &   2,0,.false.,.false.,2*nHistograms(analysis)-2,working,hist_derivs(:,i,analysis),
     &   hist_sigma(:,i,analysis),IER)
        if (IER .lt. 0) then
          write(*,*) 'Error in nulike_edispinit: TSPSI failed with error'
          write(*,*) 'code',IER, ' at i=',i
          stop
        endif

      enddo 


      end subroutine nulike_edispinit

