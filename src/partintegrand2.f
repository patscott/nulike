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
      real*8 nulike_edisp, nulike_psf, nulike_sens, nulike_angres
      integer ptype
      external dsdxdy

      Elep = Eshare * (1.d0 - y)                                               ! GeV
      if (Elep .le. epsilon(Elep)) then                                        ! Abort if the lepton energy is zero.
        nulike_partintegrand2 = 0.d0                                                
        return
      endif
      log10Elep = dlog10(Elep)                                       

      !Effective volume
      effvol = nulike_sens(log10Elep, ptypeshare)                              ! km^3
      if (effvol .le. epsilon(effvol)) then                                    ! Abort if effective volume is
        nulike_partintegrand2 = 0.d0                                                
        return
      endif

      !Energy dispersion
      edisp = nulike_edisp(log10Elep, events_ee(eventnumshare), 2014)          ! [ee]^-1
      if (edisp .le. epsilon(effvol)) then                                     ! Abort if effective volume is
        nulike_partintegrand2 = 0.d0                                                
        return
      endif

      !Angular loss factor
      errlep = nulike_angres(log10Elep)                                        ! degrees
      angloss = 1.d0-exp(-0.5d0*phi_max*phi_max/(errlep*errlep))               ! dimensionless

      !Cross-sections, densities and and point-spread functions
      ptype = ptypeshare + 2*(leptypeshare-1)                                  ! Convert to nusigma internal neutrino number.
      spnc = dsdxdy(Eshare,xshare,y,ptype,'p','NC')  ! cm^2
      spcc = dsdxdy(Eshare,xshare,y,ptype,'p','CC')  ! cm^2
      snnc = dsdxdy(Eshare,xshare,y,ptype,'n','NC')  ! cm^2
      sncc = dsdxdy(Eshare,xshare,y,ptype,'n','CC')  ! cm^2

      mlep2 = lepmass(leptypeshare)*lepmass(leptypeshare)                      ! GeV^2
      plep = dsqrt(Elep*Elep - mlep2)                                          ! GeV
      phi_obs = acos(events_cosphi(eventnumshare,analysis)) * 180.d0/pi        ! --> degrees
      phi_err = events_cosphiErr(eventnumshare,analysis)                       ! Already in degrees.

      if (spnc .gt. 0.d0 .or. spcc .gt. 0.d0) then
        philep_p = acos((Elep - m_p*xshare*y - mlep2/Eshare) / plep) * 180.d0/pi ! --> degrees
        pcont = numdens_p * nulike_psf(phi_obs, philep_p, phi_max, phi_err) * (spnc + spcc) ! 1e-5 m^-3 cm^2 deg^-1
      else 
        pcont = 0.d0
      endif

      if (snnc .gt. 0.d0 .or. sncc .gt. 0.d0) then
        philep_n = acos((Elep - m_n*xshare*y - mlep2/Eshare) / plep) * 180.d0/pi ! --> degrees
        ncont = numdens_n * nulike_psf(phi_obs, philep_n, phi_max, phi_err) * (snnc + sncc) ! 1e-5 m^-3 cm^2 deg^-1
      else 
        ncont = 0.d0
      endif
      
      !Total
      nulike_partintegrand2 = edisp * effvol * angloss * (pcont + ncont)  ! km^3 1e-5 m^-3 cm^2 deg^-1 [ee]^-1-->m^2 deg^-1 [ee]^-1

      if (nulike_partintegrand2 .ne. nulike_partintegrand2) then
        write(*,*) spnc, spcc, snnc, sncc
        write(*,*) Eshare,xshare,y,ptype
        write(*,*) edisp, effvol, angloss
        write(*,*) pcont, ncont
        write(*,*) phi_obs, philep_p, philep_n, phi_max, phi_err
        write(*,*) Elep, lepmass(leptypeshare), plep, Eshare, mlep2, m_p, m_n
        write(*,*) (Elep - m_p*xshare*y - mlep2/Eshare) / plep
        stop 'found NaN'
      endif

      end function nulike_partintegrand2
