***********************************************************************
*** nulike_partialintegrand provides the integrand for the double
*** integral required to be computed to obtain partial likelihoods or
*** effective areas.
*** This routine is used only with the 2015 likelihood.
***
*** Input:              x              Bjorken x
***                     y              Bjorken y
***                     dsdxdy         cross-section function
***                                    (see partials.f for details)
***                     eventnum       the unique index number of this event.
***                                    0  => get integrand for effective area
***                                    -1 => get integrand for effective area,
***                                          without angular loss factor
***                     E              neutrino energy (GeV)
***                     ptype          1=nu, 2=nubar
***                     leptype        1=e, 2=mu, 3=tau          
***
*** Output:             integrand      m^2 degrees^-1 [energy estimator]^-1
***
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Jun 17, 2014
***********************************************************************


      real*8 function nulike_partintegrand(x,y,dsdxdy,eventnum,E,ptype,leptype)

      use iso_c_binding, only: c_ptr
      use MarcumQ

      implicit none
      include 'nulike_internal.h'
      include 'nuprep.h'

      real*8 x, y, E, Elep, log10Elep, mlep2, phi_obs, phi_err, errlep
      real*8 spcc, sncc, pcont, ncont, edisp, effvol, angloss, dsdxdy 
      real*8 nulike_edisp, nulike_offctrpsf, nulike_sens, nulike_angres, cosang
      real*8 twoerrlep2, arg1, arg2, q, plep, philep_p, philep_n
      integer eventnum, ptype, ptype_ns, leptype, ierr1, ierr2, ierr3, ierr4
      external dsdxdy

      Elep = E * (1.d0 - y)                                                    ! GeV
      if (Elep .le. 0.d0) then                                                 ! Abort if the lepton energy is zero.
        nulike_partintegrand = 0.d0                                                
        return
      endif
      log10Elep = dlog10(Elep)                                       

      !Effective volume
      effvol = nulike_sens(log10Elep, 3)                                       ! km^3
      if (effvol .le. 0.d0) then                                               ! Abort if effective volume is
        nulike_partintegrand = 0.d0                                            ! zero.      
        return
      endif

      twoerrlep2 = 0.d0
      arg2 = 0.d0
      if (eventnum .eq. 0) then
        !Relevant quantities for the angular loss factor
        errlep = nulike_angres(log10Elep)                                      ! degrees
        twoerrlep2 = 2.d0*errlep*errlep                                        ! degrees^2
        arg2 = phi_max*phi_max/twoerrlep2                                      ! dimensionless
      elseif (eventnum .gt. 0) then
        !Energy dispersion
        edisp = nulike_edisp(log10Elep, events_ee(eventnum, analysis), 2015)   ! [ee]^-1
        if (edisp .le. epsilon(edisp)) then                                    ! Abort if energy dispersion is
          nulike_partintegrand = 0.d0                                          ! zero.     
          return
        endif
      endif

      !Cross-sections and densities 
      ptype_ns = ptype + 2*(leptype-1)                                         ! Convert to nusigma internal neutrino number.
      spcc = dsdxdy(E,x,y,ptype_ns,'p','CC')                                   ! cm^2
      sncc = dsdxdy(E,x,y,ptype_ns,'n','CC')                                   ! cm^2
      pcont = numdens_p * spcc                                                 ! 1e-5 m^-3 cm^2 
      ncont = numdens_n * sncc                                                 ! 1e-5 m^-3 cm^2

      !Prime errors 
      ierr1 = 0
      ierr2 = 0
      ierr3 = 0
      ierr4 = 0

      !Lepton properties
      mlep2 = lepmass(leptype)*lepmass(leptype)                              ! GeV^2
      plep = dsqrt(Elep*Elep - mlep2)                                        ! GeV

      !Event properties
      if (eventnum .gt. 0) then
        phi_obs = acos(events_cosphi(eventnum,analysis)) * 180.d0/pi         ! --> degrees
        phi_err = events_cosphiErr(eventnum,analysis)                        ! Already in degrees.
      endif

      !Point-spread/angular loss function for proton interactions
      if (pcont .gt. 0.d0) then
        if (eventnum .ge. 0) then
          !Lepton opening angle
          cosang = (Elep - m_p*x*y - 0.5d0*mlep2/E) / plep                     ! Cosine of lepton scattering angle
          if (cosang .gt. 1.d0) cosang = 1.d0                                  ! Fix floating-pt errs (>1 kinematically disallowed)
          philep_p = acos(cosang) * 180.d0/pi                                  ! --> degrees
          if (eventnum .eq. 0) then
            !Angular loss factor
            arg1 = philep_p*philep_p / twoerrlep2                              ! dimensionless
            if (arg1 .lt. 1.d4 .and. arg2 .lt. 1.d4) then
              call marcum(1.d0,arg1,arg2,angloss,q,ierr1)
            else ! If either arg is huge, the PSF width is tiny compared to either the lepton angle or the size of the cut cone,
                 ! so the PSF is either fully in or fully out of the cut cone.
              if (arg1 .gt. arg2) then
                angloss = 0.d0 ! Lepton angle is bigger than cut cone --> psf outside cut cone
              else 
                angloss = 1.d0 ! Lepton angle is less than cut cone --> psf fully inside cut cone.
              endif
            endif          
            pcont = pcont * angloss 
          else
            !PSF and likelihood
            pcont = pcont * nulike_offctrpsf(phi_obs, philep_p, phi_err) ! 1e-5 m^-3 cm^2 deg^-1
          endif
        endif
      else 
        pcont = 0.d0
      endif

      !Point-spread/angular loss function for neutron interactions
      if (ncont .gt. 0.d0) then
        if (eventnum .ge. 0) then
          !Lepton opening angle
          cosang = (Elep - m_n*x*y - 0.5d0*mlep2/E) / plep                     ! Cosine of lepton scattering angle
          if (cosang .gt. 1.d0) cosang = 1.d0                                  ! Fix floating-pt errs (>1 kinematically disallowed)
          philep_n = acos(cosang) * 180.d0/pi                                  ! --> degrees
          if (eventnum .eq. 0) then
            !Angular loss factor
            arg1 = philep_n*philep_n / twoerrlep2                              ! dimensionless
            if (arg1 .lt. 1.d4 .and. arg2 .lt. 1.d4) then
              call marcum(1.d0,arg1,arg2,angloss,q,ierr3)
            else ! If either arg is huge, the PSF width is tiny compared to either the lepton angle or the size of the cut cone,
                 ! so the PSF is either fully in or fully out of the cut cone.
              if (arg1 .gt. arg2) then
                angloss = 0.d0 ! Lepton angle is bigger than cut cone --> psf outside cut cone
              else 
                angloss = 1.d0 ! Lepton angle is less than cut cone --> psf fully inside cut cone.
              endif
            endif
            ncont = ncont * angloss
          else
            !PSF and likelihood
            ncont = ncont * nulike_offctrpsf(phi_obs, philep_n, phi_err)       ! 1e-5 m^-3 cm^2 deg^-1
          endif
        endif
      else 
        ncont = 0.d0
      endif
      
      !Total
      nulike_partintegrand = effvol * (pcont + ncont)                         ! km^3 1e-5 m^-3 cm^2 deg^-1 -->m^2 deg^-1
      if (eventnum .gt. 0) nulike_partintegrand = nulike_partintegrand*edisp  ! -->m^2 deg^-1 [ee]^-1

      !Debug
      if (ierr1 .gt. 1 .or. ierr2 .gt. 1 .or. ierr3 .gt. 1 .or. ierr4 .gt. 1) then
        write(*,*) 'philep_p, philep_n, arg1, arg2:',philep_p, philep_n, arg1, arg2
        write(*,*) 'ierr = ',ierr1,ierr2,ierr3,ierr4
        stop 'Catastrophic error when calling marcum in nulike_partintegrand!'
      endif
      if (nulike_partintegrand .ne. nulike_partintegrand) then
        write(*,*) 'Printing debug info:'
        write(*,*) E,x,y
        write(*,*) Elep, effvol, angloss
        write(*,*) ptype_ns,eventnum,leptype
        write(*,*) spcc, sncc
        write(*,*) pcont, ncont
        if (eventnum .gt. 0) then
          write(*,*) plep, mlep2
          write(*,*) edisp, (Elep - m_p*x*y - mlep2/E) / plep
          write(*,*) phi_obs, phi_max, phi_err
          write(*,*) philep_n, philep_p
          write(*,*) nulike_offctrpsf(phi_obs, philep_n, phi_err), nulike_offctrpsf(phi_obs, philep_p, phi_err)
        endif
        stop 'NaN found in nulike_partintegrand!'
      endif

      end function nulike_partintegrand
