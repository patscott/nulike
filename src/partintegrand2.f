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

      use MarcumQ

      implicit none
      include 'nulike_internal.h'
      include 'nuprep.h'

      real*8 y, Elep, log10Elep, mlep2, phi_obs, phi_err, errlep, plep, philep_p, philep_n
      real*8 spnc, spcc, snnc, sncc, pcont, ncont, edisp, effvol, angloss, dsdxdy 
      real*8 nulike_edisp, nulike_offctrpsf, nulike_sens, nulike_angres, cosang
      real*8 twoerrlep2, arg1, arg2, arg3, pupper, plower, q
      integer ptype, ierr1, ierr2, ierr3, ierr4
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

      !Relevant quantities for the angular loss factor
      errlep = nulike_angres(log10Elep)                                        ! degrees
      twoerrlep2 = 2.d0*errlep*errlep                                          ! degrees^2
      arg2 = phi_max*phi_max/twoerrlep2                                        ! dimensionless
      arg3 = 3.24d4/twoerrlep2                                                 ! (180deg)^2/2sigma^2 = dimensionless

      !Cross-sections and densities 
      ptype = ptypeshare + 2*(leptypeshare-1)                                  ! Convert to nusigma internal neutrino number.
      spnc = dsdxdy(Eshare,xshare,y,ptype,'p','NC')                            ! cm^2
      spcc = dsdxdy(Eshare,xshare,y,ptype,'p','CC')                            ! cm^2
      snnc = dsdxdy(Eshare,xshare,y,ptype,'n','NC')                            ! cm^2
      sncc = dsdxdy(Eshare,xshare,y,ptype,'n','CC')                            ! cm^2
      pcont = numdens_p * (spnc + spcc)                                        ! 1e-5 m^-3 cm^2 
      ncont = numdens_n * (snnc + sncc)                                        ! 1e-5 m^-3 cm^2

      !Prime errors 
      ierr1 = 0
      ierr2 = 0
      ierr3 = 0
      ierr4 = 0

      !Point-spread functions
      if (eventnumshare .ne. 0) then

        mlep2 = lepmass(leptypeshare)*lepmass(leptypeshare)                    ! GeV^2
        plep = dsqrt(Elep*Elep - mlep2)                                        ! GeV
        phi_obs = acos(events_cosphi(eventnumshare,analysis)) * 180.d0/pi      ! --> degrees
        phi_err = events_cosphiErr(eventnumshare,analysis)                     ! Already in degrees.

        if (pcont .gt. 0.d0) then
          !Lepton opening angle
          cosang = (Elep - m_p*xshare*y - 0.5d0*mlep2/Eshare) / plep           ! Cosine of lepton scattering angle
          if (cosang .gt. 1.d0) cosang = 1.d0                                  ! Fix floating-pt errs (>1 kinematically disallowed)
          philep_p = acos(cosang) * 180.d0/pi                                  ! --> degrees
          !Angular loss factor
          arg1 = philep_p*philep_p / twoerrlep2                                ! dimensionless
          call marcum(1.d0,arg1,arg2,pupper,q,ierr1)
          if (arg1 .lt. 1.d5 .and. arg2 .lt. 1.d5) then
            call marcum(1.d0,arg1,arg2,pupper,q,ierr1)
          else ! If either arg is huge, the PSF width is tiny compared to either the lepton angle or the size of the cut cone,
               ! so the PSF is either fully in or fully out of the cut cone.
            if (arg1 .gt. arg2) then
              pupper = 0.d0 ! Lepton angle is bigger than cut cone --> psf outside cut cone
            else 
              pupper = 1.d0 ! Lepton angle is less than cut cone --> psf fully inside cut cone.
            endif
          endif
          if (arg3 .lt. 1.d5) then
            call marcum(1.d0,arg1,arg3,plower,q,ierr2)
          else
            plower = 1.d0
          endif
          angloss = pupper/plower                                              ! dimensionless
          !Likelihood
          pcont = pcont * angloss * nulike_offctrpsf(phi_obs, philep_p, phi_err) ! 1e-5 m^-3 cm^2 deg^-1
        else 
          pcont = 0.d0
        endif

        if (ncont .gt. 0.d0) then
          !Lepton opening angle
          cosang = (Elep - m_n*xshare*y - 0.5d0*mlep2/Eshare) / plep           ! Cosine of lepton scattering angle
          if (cosang .gt. 1.d0) cosang = 1.d0                                  ! Fix floating-pt errs (>1 kinematically disallowed)
          philep_n = acos(cosang) * 180.d0/pi                                  ! --> degrees
          !Angular loss factor
          arg1 = philep_n*philep_n / twoerrlep2                                ! dimensionless
          if (arg1 .lt. 1.d5 .and. arg2 .lt. 1.d5) then
            call marcum(1.d0,arg1,arg2,pupper,q,ierr3)
          else ! If either arg is huge, the PSF width is tiny compared to either the lepton angle or the size of the cut cone,
               ! so the PSF is either fully in or fully out of the cut cone.
            if (arg1 .gt. arg2) then
              pupper = 0.d0 ! Lepton angle is bigger than cut cone --> psf outside cut cone
            else 
              pupper = 1.d0 ! Lepton angle is less than cut cone --> psf fully inside cut cone.
            endif
          endif
          if (arg3 .lt. 1.d5) then
            call marcum(1.d0,arg1,arg3,plower,q,ierr4)
          else
            plower = 1.d0 ! If upper cutoff is huge, there is no renormalisation needed.
          endif
          angloss = pupper/plower                                              ! dimensionless
          !Likelihood
          ncont = ncont * angloss * nulike_offctrpsf(phi_obs, philep_n, phi_err) ! 1e-5 m^-3 cm^2 deg^-1
        else 
          ncont = 0.d0
        endif
 
      endif
      
      !Total
      nulike_partintegrand2 = effvol * (pcont + ncont)                         ! km^3 1e-5 m^-3 cm^2 deg^-1 -->m^2 deg^-1
      if (eventnumshare .ne. 0) nulike_partintegrand2 = nulike_partintegrand2*edisp  ! -->m^2 deg^-1 [ee]^-1


      !Debug
      if (ierr1 .gt. 1 .or. ierr2 .gt. 1 .or. ierr3 .gt. 1 .or. ierr4 .gt. 1) then
        write(*,*) 'philep_p, philep_n, arg1, arg2, arg3:',philep_p, philep_n, arg1, arg2, arg3
        write(*,*) 'ierr = ',ierr1,ierr2,ierr3,ierr4
        stop 'Catastrophic error when calling marcum in nulike_partintegrand2!'
      endif
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
          write(*,*) nulike_offctrpsf(phi_obs, philep_n, phi_err), nulike_offctrpsf(phi_obs, philep_p, phi_err)
        endif
        stop 'NaN found in nulike_partintegrand2!'
      endif

      end function nulike_partintegrand2
