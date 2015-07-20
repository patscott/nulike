***********************************************************************
*** nulike_tabulated_weight provides access to the the precomputed 
*** weights to be used in
*** a) the signal component of the combined angular-spectral likelihood
***    in nulike_specanglike
*** b) the total predicted number of signal events in nulike_signal.
*** This routine is used only with the 2015 likelihood.
***
*** Input:    log10E       log10(neutrino energy/GeV)
***           ptype        1=nu, 2=nubar
***           eventnum     the unique index number of the event
***
*** Output:   weight       eventnum = 0 => effective area (m^2)
***                        eventnum > 0 => weighting for event
***                         eventnum (m^2 chan^-1 degrees^-1 )
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Jun 8, 2014
***********************************************************************


      real*8 function nulike_tabulated_weight(log10E,ptype,eventnum)

      use iso_c_binding, only: c_ptr

      implicit none
      include 'nulike_internal.h'

      real*8 log10E, nulike_tabulated_weight_a(1)
      integer ptype, eventnum, IER

      if (eventnum .eq. 0) then 

        !Call interpolator to get effective area for this energy
        call TSVAL1(nPrecompE(analysis)-start_index(analysis)+1,
     &   precomp_log10E(start_index(analysis):,analysis),
     &   precompEA_weights(start_index(analysis):,ptype,analysis),
     &   precompEA_derivs(:,ptype,analysis),
     &   precompEA_sigma(:,ptype,analysis),
     &   0,1,log10E,nulike_tabulated_weight_a,IER)
        nulike_tabulated_weight = nulike_tabulated_weight_a(1)

      else

        !If this event is just impossible for a signal to produce (e.g. it has ee where the probability
        !of getting this ee from signal is zero), then just return zero.
        if (all(precomp_weights(:,eventnum,ptype,analysis) - logZero .le. epsilon(logZero))) then
          nulike_tabulated_weight = 0.d0
          return
        endif
          
        !Call interpolator for this event to get weight for this energy
        call TSVAL1(nPrecompE(analysis),precomp_log10E(:,analysis),
     &   precomp_weights(:,eventnum,ptype,analysis),
     &   precomp_derivs(:,eventnum,ptype,analysis),
     &   precomp_sigma(:,eventnum,ptype,analysis),
     &   0,1,log10E,nulike_tabulated_weight_a,IER)
        nulike_tabulated_weight = nulike_tabulated_weight_a(1)

      endif

      if (IER .lt. 0) then
        write(*,*) 'TSVAL1 error from weight or effective area'
        write(*,*) 'in nulike_tabulated_weight, code:', IER
        stop
      endif

      if (nulike_tabulated_weight - logZero .lt. epsilon(logZero)) then
        nulike_tabulated_weight = 0.d0
      else
        nulike_tabulated_weight = 10.d0**nulike_tabulated_weight
      endif

      end function nulike_tabulated_weight
