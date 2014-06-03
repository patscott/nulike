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
*** Modified: Jun 3, 2014
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
      real*8  working(2*nHistograms-2)

      !Read in nchan response distribution for each incoming neutrino energy band
      open(lun,file=filename, ACTION='READ')

      !Skip over header
      instring = '#'
      do while (instring .eq. '#')
        read(lun, fmt='(A1)'), instring
      enddo

      !Read actual data
      do i = 1, nbins_Ein
        read(lun, *) instring, hist_logE(1,i), hist_logE(2,i)
        hist_logEcentres(i) = 0.5d0*(hist_logE(1,i)+hist_logE(2,i))
        do j = 1, nbins_nchan(i)
          read(lun, *), instring, dummyint,
     &     hist_nchan_temp(i,j), hist_prob_temp(i,j)
        enddo
        read(lun,*) instring
        if (i .ne. nbins_Ein) read(lun,*) instring
      enddo

      close(lun)

      !Arrange histograms so they all cover the same range in nchan
      hist_prob = 0.d0
      do k = 1, nnchan_total
        do i = 1, nbins_Ein
          hist_nchan(i,k) = k - 1 + nchan_min
          do j = 1, nbins_nchan(i)
            if (hist_nchan_temp(i,j) .eq. hist_nchan(i,k)) then
              hist_prob(i,k) = hist_prob_temp(i,j)
            endif
          enddo
        enddo
      enddo

      !Work out where the indexing of nchan values in energy dispersion lines
      !up with indexing of nchan values in observed background spectrum.
      nchan_hist2BGoffset = -1
      do k = 1, nnchan_total
        if (hist_nchan(1,k) .eq. BGnchandist_nchan(1)) then 
          nchan_hist2BGoffset = k-1
        endif
      enddo

      !Set up interpolation in energy histograms for use as energy dispersion estimator
      do i = 1, nnchan_total

        do j = 1, nHistograms
          edisp_prob(j) = hist_prob(j,i)
        enddo

        call TSPSI(nHistograms,hist_logEcentres,edisp_prob,
     &   2,0,.false.,.false.,2*nHistograms-2,working,edisp_derivs,
     &   edisp_sigma,IER)
        if (IER .lt. 0) then
          write(*,*) 'Error in nulike_edispinit: TSPSI failed with error'
          write(*,*) 'code',IER, ' at i=',i
          stop
        endif

        do j = 1, nHistograms
          hist_derivs(j,i) = edisp_derivs(j)
          hist_sigma(j,i) = edisp_sigma(j)
        enddo

      enddo 

      !Indicate which nchan data are currently loaded for
      nchansaved = nnchan_total + nchan_min - 1


      end subroutine nulike_edispinit

