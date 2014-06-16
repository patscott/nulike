***********************************************************************
*** nulike_specanginit initialises the precomputed effective area and 
*** partial angular-spectral likelihoods by reading them in from the 
*** specified file.
*** This routine is used only with the 2014 likelihood.
***
*** input:  filename     name of file containig partial likelihoods.
***         nevents      number of events for which partial likelihoods 
***                       are to be read in.
***         nenergies    number of energies for which partial likelihoods
***                       are to be read in for each event.
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Jun 15 2014
***********************************************************************

      subroutine nulike_specanginit(filename, n_events, phi_cut, nenergies)

      implicit none
      include 'nulike.h'

      character (len=*) filename
      integer n_events, nenergies
      real*8 phi_cut


      end subroutine nulike_specanginit
