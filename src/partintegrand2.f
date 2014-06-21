***********************************************************************
*** nulike_partialintegrand2 provides the integrand for the inner double
*** integral required to be computed to obtain partial likelihoods.
*** This routine is used only with the 2014 likelihood.
***
*** Input:              y              Bjorken y
***                     dsdxdy         cross-section function
***                                    (see partials.f for details)
***
*** Hidden Input:       eventnumshare  the unique index number of this event
***                     xshare         Bjorken x
***                     Eshare         neutrino energy (GeV)
***                     log10Eshare    log_10(neutrino energy / GeV)
***                     ptypeshare     1=nu, 2=nubar
***                     leptypeshare   1=e, 2=mu, 3=tau          
***
*** Output:             integrand      m^2 degrees^-1 [energy estimator]^-1
***
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Jun 17, 2014
***********************************************************************


      real*8 function nulike_partintegrand2(y,dsdxdy)

      implicit none
      include 'nulike.h'
      include 'nuprep.h'

      real*8 y, Elep, log10Elep, mlep2, phi_obs, phi_err, errlep, plep, philep_p, philep_n
      real*8 spnc, spcc, snnc, sncc, pcont, ncont, edisp, effvol, angloss, dsdxdy 
      real*8 nulike_edisp, nulike_psf, nulike_sens, nulike_angres, cosang
      integer ptype
      external dsdxdy

      Elep = Eshare * (1.d0 - y)                                               ! GeV
      if (Elep .le. 0.d0) then                                                 ! Abort if the lepton energy is zero.
        nulike_partintegrand2 = 0.d0                                                
        return
      endif
      log10Elep = dlog10(Elep)                                       

      !Effective volume
      effvol = nulike_sens(log10Elep, ptypeshare)                              ! km^3
      if (effvol .le. 0.d0) then                                               ! Abort if effective volume is
        nulike_partintegrand2 = 0.d0                                           ! zero.      
        return
      endif

      !Energy dispersion
      if (eventnumshare .ne. 0) then
        edisp = nulike_edisp(log10Elep, events_ee(eventnumshare, analysis), 2014)! [ee]^-1
        if (edisp .le. epsilon(effvol)) then                                     ! Abort if energy dispersion is
          nulike_partintegrand2 = 0.d0                                           ! zero.     
          return
        endif
      endif

      !Angular loss factor
      errlep = nulike_angres(log10Elep)                                        ! degrees
      angloss = 1.d0-exp(-0.5d0*phi_max*phi_max/(errlep*errlep))               ! dimensionless

      !Cross-sections and densities 
      ptype = ptypeshare + 2*(leptypeshare-1)                                  ! Convert to nusigma internal neutrino number.
      spnc = dsdxdy(Eshare,xshare,y,ptype,'p','NC')                            ! cm^2
      spcc = dsdxdy(Eshare,xshare,y,ptype,'p','CC')                            ! cm^2
      snnc = dsdxdy(Eshare,xshare,y,ptype,'n','NC')                            ! cm^2
      sncc = dsdxdy(Eshare,xshare,y,ptype,'n','CC')                            ! cm^2
      pcont = numdens_p * (spnc + spcc)                                        ! 1e-5 m^-3 cm^2 
      ncont = numdens_n * (snnc + sncc)                                        ! 1e-5 m^-3 cm^2

      !Point-spread functions
      if (eventnumshare .ne. 0) then

        mlep2 = lepmass(leptypeshare)*lepmass(leptypeshare)                    ! GeV^2
        plep = dsqrt(Elep*Elep - mlep2)                                        !GeV
        phi_obs = acos(events_cosphi(eventnumshare,analysis)) * 180.d0/pi      ! --> degrees
        phi_err = events_cosphiErr(eventnumshare,analysis)                     ! Already in degrees.

        if (pcont .gt. 0.d0) then
          cosang = (Elep - m_p*xshare*y - mlep2/Eshare) / plep
          if (cosang .gt. 1.d0) cosang = 1.d0 ! Fix floating-point errors (>1 is kinematically disallowed)
          philep_p = acos(cosang) * 180.d0/pi ! --> degrees
          pcont = pcont * nulike_psf(phi_obs, philep_p, phi_err) ! 1e-5 m^-3 cm^2 deg^-1
        else 
          pcont = 0.d0
        endif

        if (ncont .gt. 0.d0) then
          cosang = (Elep - m_n*xshare*y - mlep2/Eshare) / plep
          if (cosang .gt. 1.d0) cosang = 1.d0 ! Fix floating-point errors (>1 is kinematically disallowed)
          philep_n = acos(cosang) * 180.d0/pi ! --> degrees
          ncont = ncont * nulike_psf(phi_obs, philep_n, phi_err) ! 1e-5 m^-3 cm^2 deg^-1
        else 
          ncont = 0.d0
        endif
 
      endif
      
      !Total
      nulike_partintegrand2 = effvol * angloss * (pcont + ncont)  ! km^3 1e-5 m^-3 cm^2 deg^-1 -->m^2 deg^-1
      if (eventnumshare .ne. 0) nulike_partintegrand2 = nulike_partintegrand2*edisp  ! -->m^2 deg^-1 [ee]^-1


      !Debug
      if (nulike_partintegrand2 .ne. nulike_partintegrand2) then
        write(*,*) 'Printing debug info:'
        write(*,*) Eshare,xshare,y
        write(*,*) Elep, effvol, angloss
        write(*,*) ptype,eventnumshare,leptypeshare
        write(*,*) spnc, spcc, snnc, sncc
        write(*,*) pcont, ncont
        if (eventnumshare .ne. 0) then
          write(*,*) plep, mlep2
          write(*,*) edisp, (Elep - m_p*xshare*y - mlep2/Eshare) / plep
          write(*,*) phi_obs, phi_max, phi_err
          write(*,*) philep_n, philep_p
          write(*,*) nulike_psf(phi_obs, philep_n, phi_err), nulike_psf(phi_obs, philep_p, phi_err)
        endif
        stop 'NaN found in nulike_partintegrand2!'
      endif

      end function nulike_partintegrand2
