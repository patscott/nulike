! This program runs the nulike routines to find limits on
! generic WIMPs, as discussed in the nulike paper. Longer
! explanations are provided in nulike_test.f.
!
! Author: Pat Scott p.scott@imperial.ac.uk
! Date: June 1 2015


      program nulike_test_wimp

      use iso_c_binding, only: c_bool

      implicit none
      !Nulike include
      include 'nulike.h'
      !DarkSUSY includes
      include 'dsio.h'
      include 'dsidtag.h'
      include 'dsntcom.h'
      include 'dsdirver.h'
      include 'dswacom.h'

      ! Nuclear scattering
      real*8 sigsdp
      ! Capture rate
      real*8 dsntcapsuntab, ca, tt_sun, annrate
      ! Neutrino likelihoods
      real*8 theoryError, diff_CL, zeroin
      logical uselogNorm
      logical BGLikePrecompute
      logical(c_bool) threadsafe
      character (len=nulike_clen) iclike2015
      character (len=nulike_clen) experiment(3), eventf
      character (len=nulike_clen) BGf, partiald, efareaf
      double precision :: ref_CL = 90.d0
      external diff_CL
      common/wimpcom/ref_CL, theoryError, experiment, threadsafe
      ! Book-keeping
      integer i,j
      logical :: first = .true.
      real*8, parameter :: dummyval = 0, mwimpmin = 10.d0, mwimpmax = 1000
      !These are the WIMP masses for which DarkSUSY contains WimpSim results; it obtains spectra for other masses by interpolating between these.
      real*8, parameter :: chosen_masses(20) = (/10.d0, 25.d0, 50.d0, 80.3d0, 91.2d0, 1.d2, 1.5d2, 1.76d2, 2.d2, 2.5d2, 3.5d2, 5.d2, 7.5d2, 1.d3, 1.5d3, 2.d3, 3.d3, 5.d3, 7.5d3, 1.d4/)
      integer, parameter :: mwimp_pts = 10
      logical, parameter :: talky = .false., chosen_ones_only = .true.

      ! See the header of src/init.f for detailed explanations of the following options.
      iclike2015 = 'data/IceCube/likelihood2015/'
      uselogNorm = .true.
      BGLikePrecompute = .true.
      threadsafe = .true.

      experiment(1) = 'IC-79 SL'
      eventf  = trim(iclike2015)//'IC79_Events_SL_llhInput_60Deg.txt'
      BGf     = trim(iclike2015)//'IC79_Background_distributions_SL.txt'
      efareaf = trim(iclike2015)//'IC79_Effective_Area_SL.txt'
      partiald= trim(iclike2015)//'IC79_Partial_Likelihoods_SL'
      call nulike_init(experiment(1), eventf, BGf, efareaf, partiald,
     & dummyval, uselogNorm, BGLikePrecompute)

      experiment(2) = 'IC-79 WL'
      eventf  = trim(iclike2015)//'IC79_Events_WL_llhInput_60Deg.txt'
      BGf     = trim(iclike2015)//'IC79_Background_distributions_WL.txt'
      efareaf = trim(iclike2015)//'IC79_Effective_Area_WL.txt'
      partiald= trim(iclike2015)//'IC79_Partial_Likelihoods_WL'
      call nulike_init(experiment(2), eventf, BGf, efareaf, partiald,
     & dummyval, uselogNorm, BGLikePrecompute)

      experiment(3) = 'IC-79 WH'
      eventf  = trim(iclike2015)//'IC79_Events_WH_llhInput_60Deg.txt'
      BGf     = trim(iclike2015)//'IC79_Background_distributions_WH.txt'
      efareaf = trim(iclike2015)//'IC79_Effective_Area_WH.txt'
      partiald= trim(iclike2015)//'IC79_Partial_Likelihoods_WH'
      call nulike_init(experiment(3), eventf, BGf, efareaf, partiald,
     & dummyval, uselogNorm, BGLikePrecompute)

      ! Initialise DarkSUSY
      if (.not. talky) then
        prtlevel=0
      else
        write(*,*)
      endif
      call dsinit
      prtlevel=0

      ! Initialise the WIMP model
      dswasetupcalled=.true.
      ! Set annihilation branching fractions
      wabr(1) = 0.d0        !H1 H1
      wabr(2) = 0.d0        !H1 H2
      wabr(3) = 0.d0        !H2 H2
      wabr(4) = 0.d0        !H3 H3
      wabr(5) = 0.d0        !H1 H3
      wabr(6) = 0.d0        !H2 H3
      wabr(7) = 0.d0        !H+ H-
      wabr(8) = 0.d0        !Z  H1
      wabr(9) = 0.d0        !Z  H2
      wabr(10) = 0.d0       !Z  H3
      wabr(11) = 0.d0       !W+ H- and W- H+
      wabr(12) = 0.d0       !Z  Z
      wabr(13) = 0.d0       !W+ W-
      wabr(14) = 0.d0       !nu_e nubar_e
      wabr(15) = 0.d0       !e+ e-
      wabr(16) = 0.d0       !nu_mu nubar_mu
      wabr(17) = 0.d0       !mu+ mu-
      wabr(18) = 0.d0       !nu_tau nubar_tau
      wabr(19) = 1.d0       !tau+ tau-
      wabr(20) = 0.d0       !u  ubar
      wabr(21) = 0.d0       !d  dbar
      wabr(22) = 0.d0       !c  cbar
      wabr(23) = 0.d0       !s  sbar
      wabr(24) = 0.d0       !t  tbar
      wabr(25) = 0.d0       !b  bbar
      wabr(26) = 0.d0       !g  g
      wabr(27) = 0.d0       !qqg (not implemented)
      wabr(28) = 0.d0       !gamma gamma
      wabr(29) = 0.d0       !Z  gamma
      ! Neutral Higgs decay branching fractions
      was0br(1,1) = 0.d0    !H1 H1
      was0br(2,1) = 0.d0    !H1 H2
      was0br(3,1) = 0.d0    !H2 H2
      was0br(4,1) = 0.d0    !H3 H3
      was0br(5,1) = 0.d0    !H1 H3
      was0br(6,1) = 0.d0    !H2 H3
      was0br(7,1) = 0.d0    !H+ H-
      was0br(8,1) = 0.d0    !Z  H1
      was0br(9,1) = 0.d0    !Z  H2
      was0br(10,1) = 0.d0   !Z  H3
      was0br(11,1) = 0.d0   !W+ H- and W- H+
      was0br(12,1) = 0.0266 !Z  Z
      was0br(13,1) = 0.216  !W+ W-
      was0br(14,1) = 0.d0   !nu_e nubar_e
      was0br(15,1) = 0.d0   !e+ e-
      was0br(16,1) = 0.d0   !nu_mu nbar_mu
      was0br(17,1) = .000221!mu+ mu-
      was0br(18,1) = 0.d0   !nu_tau nubar_tau
      was0br(19,1) = 0.0637 !tau+ tau-
      was0br(20,1) = 0.d0   !u  ubar
      was0br(21,1) = 0.d0   !d  dbar
      was0br(22,1) = 0.0267 !c  cbar
      was0br(23,1) = .000439!s  sbar
      was0br(24,1) = 0.d0   !t  tbar
      was0br(25,1) = 0.577  !b  bbar
      was0br(26,1) = 0.0855 !g  g
      was0br(27,1) = 0.d0   !qqg (not implemented)
      was0br(28,1) = 0.00229!gamma gamma
      was0br(29,1) = 0.00155!Z  gamma
      do i=2,3
         do j=1,29
            was0br(j,i)=0.d0
         enddo
      enddo
      ! Neutral Higgs masses
      was0m(1) = 125.d0
      was0m(2) = 0.d0
      was0m(3) = 0.d0
      ! Charged Higgs branching fractions
      do j=1,15
         wascbr(j)=0.d0
      enddo
      ! Charged Higgs masses
      wascm=0.d0

      ! Step through each of the requested WIMP masses
      do i = 1, merge(size(chosen_masses),mwimp_pts,chosen_ones_only)

        ! Set WIMP mass and annihilation cross-section, then write them out
        if (chosen_ones_only) then
          wamwimp = chosen_masses(i)
        else
          wamwimp = 10.**(log10(mwimpmin) + dble(i-1)/dble(mwimp_pts-1)*log10(mwimpmax/mwimpmin))
        endif
        wasv = 3.d-26
        theoryError = 5d-2 * merge(1.d0, dsqrt(wamwimp*1d-2), wamwimp .le. 100.d0)
        if (talky) write(*,*) '  WIMP mass = ', wamwimp
        if (talky) write(*,*) '  Annihilation cross-section = ',wasv,' cm^-3 s^-1'
        if (talky) write(*,*) '  Theory error = ',theoryError,'%'

        ! Calculate the SD proton scattering cross-section as the 90% CL upper limit
        if (talky) write (*,*) 'Calculating SD proton scattering cross section...'
        sigsdp = 10.d0**zeroin(-45.d0, -30.d0, diff_CL, 1.d-5)
        if (talky) write(*,*) '  sigsdp (pb) = ',sigsdp*1.0d36

        ! Recalculate capture and annihlation rates with DarkSUSY, using the
        ! default Gouldian calculation with a numerical integration over
        ! the velocity distribution, sped up by tabulation.
        if (talky) write(*,*) 'Calculating capture rate in the Sun'
        csu=dsntcapsuntab(wamwimp,0.d0,sigsdp)   ! Capture rate (s^-1)
        ca=wasv/6.6d28*(wamwimp/20.d0)**(3./2.)
        tausu=1.0d0/dsqrt(csu*ca)                ! Equilibration time (s)
        tt_sun=1.5d17*dsqrt(csu*ca)
        annrate=csu*0.5d0*tanh(tt_sun)**2        ! Annihiliation rate (s^-1)

        if (talky) then
          write(*,*) '  Capture rate in the Sun = ',csu,' s^-1'
          write(*,*) '  Annihilation rate in the Sun = ',annrate,' s^-1'
        else
          if (first) then
            first = .false.
            write(*,*) '#mWIMP (GeV)   sigmav (cm^-3 s^-1)',
     &                 ' sig_SDp (cm^2) C (s^-1)   A (s^-1)'
          endif
          write(*,'(5(1x,e14.8))') wamwimp, wasv, sigsdp, csu, annrate
        endif

      enddo

      write(*,*)
      write(*,*) 'The nulike WIMP benchmark tests completed successfully.'
      write(*,*)

      end program nulike_test_wimp


      ! Difference between the sought CL and the actual one for a given sigma_SDp
      real*8 function diff_CL(log10sigma)

      use iso_c_binding, only: C_NULL_PTR, c_bool, c_ptr, c_int

      implicit none
      include 'nulike.h'
      include 'dsntcom.h'
      include 'dswacom.h'

      real*8 ca, annrate, log10sigma, sigpred, bgpred
      real*8 dsntcapsuntab, tt_sun, nuyield_test, DGAMMA, DGAMIC
      real*8 lnLike(3), pval(3), refLike(3), dof, theoryError
      double precision :: ref_CL
      integer totobs, likechoice, i
      integer(c_int) speed
      logical(c_bool) pvalFromRef
      logical(c_bool) threadsafe
      character (len=nulike_clen) experiment(3)
      type(c_ptr) ptr
      external nuyield_test
      common/wimpcom/ref_CL, theoryError, experiment, threadsafe

      ! Set the cross-sections
      wasigsip = 0.d0
      wasigsdp = 10.d0**log10sigma

      ! Calculate capture and annihlation rates with DarkSUSY, using the
      ! default Gouldian calculation with a numerical integration over
      ! the velocity distribution, sped up by tabulation.
      csu=dsntcapsuntab(wamwimp,wasigsip,wasigsdp)   ! Capture rate (s^-1)
      ca=wasv/6.6d28*(wamwimp/20.d0)**(3./2.)
      tausu=1.0d0/dsqrt(csu*ca)                ! Equilibration time (s)
      tt_sun=1.5d17*dsqrt(csu*ca)
      annrate=csu*0.5d0*tanh(tt_sun)**2        ! Annihiliation rate (s^-1)

      ! See the header of src/nulike_bounds.f for detailed explanations of the following options.
      speed = 0
      pvalFromRef = .true.
      dof = 1.0  !Calculate exclusion assuming conditioning on everything except a single parameter (eg sigma_SDp).
      ptr = C_NULL_PTR
      likechoice = 4
      !All obtained by setting signal content to zero => total high-scale decoupling.
      if (likechoice .eq. 1) then
        pvalFromRef = .false.
        refLike(1) = -5.8244
        refLike(2) = -4.7500
        refLike(3) = -4.8768
      else
        refLike(1) =  -5015.6474
        refLike(2) =  -1813.4503
        refLike(3) = -11874.8689
      endif

      do i = 1, 3
        !Use nulike to get signal and background predictions, number of observed events, likelihood and p-value
        call nulike_bounds(experiment(i), wamwimp, annrate, nuyield_test, sigpred, bgpred,
     &   totobs, lnLike(i), pval(i), likechoice, theoryError, speed, pvalFromRef,
     &   refLike(i), dof, ptr, threadsafe)
      enddo

      diff_CL = ref_CL - 1.d2*(1.d0-DGAMIC(dof*0.5d0,max(0.d0, sum(refLike-lnLike)))/DGAMMA(dof*0.5d0))

      end


      ! Function returning neutrino flux at detector.
      real*8 function nuyield_test(log10E,ptype,context)
      use iso_c_binding, only: c_ptr
      implicit none
      real*8 log10E, dsntmuonyield
      integer ptype, istat
      external dsntmuonyield
      type(c_ptr) :: context, dummy
      if (.false.) dummy = context
      nuyield_test = 1.d-30 * dsntmuonyield(10.d0**log10E,10.d0,'su',3,1,ptype,istat)
      end
