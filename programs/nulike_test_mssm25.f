! This program runs the nulike routines for the MSSM-25 benchmarks
! discussed in the nulike paper.  Note that it must be used with 
! DarkSUSY 5.1.3 or later, as earlier versions contain a bug
! in dsgive_model25, and have poor interpolation routines in neutralino
! mass and neutrino flux.  Longer explanations are provided in
! nulike_test.f.
!
! Author: Pat Scott p.scott@imperial.ac.uk
! Date: May 30 2015


      program nulike_test_mssm25

      use iso_c_binding, only: C_NULL_PTR, c_bool, c_ptr, c_int

      implicit none
      !Nulike include
      include 'nulike.h'
      !DarkSUSY includes
      include 'dsio.h'
      include 'dsidtag.h'
      include 'dsntcom.h'
      include 'dsdirver.h'
      include 'dswacom.h'

      ! Relic density
      real*8 oh2,xf,dsrdomega
      ! Nuclear scattering
      real*8 sigsip,sigsin,sigsdp,sigsdn
      ! Capture rate
      real*8 dsntcapsuntab, ca, tt_sun, annrate, nuyield_test
      ! Neutrino likelihoods
      real*8 sigpred(3), bgpred(3), lnLike(3), pval(3) 
      real*8 theoryError, refLike(3), dof, DGAMMA, DGAMIC
      integer totobs(3), likechoice
      integer(c_int) speed
      logical uselogNorm
      logical(c_bool) pvalFromRef
      logical(c_bool) threadsafe
      logical BGLikePrecompute
      type(c_ptr) ptr
      character (len=nulike_clen) iclike2015
      character (len=nulike_clen) experiment(3), eventf
      character (len=nulike_clen) BGf, partiald, efareaf
      external nuyield_test
      ! Book-keeping
      integer unphys,hwarning,iend,ierr,iwar,nfc,i,j
      logical :: on_spoke_1 = .false., on_spoke_2 = .false.
      logical :: first = .true., finished = .false.
      real*8 :: x = 0.d0, y = 0.d0, excl, rs
      real*8, parameter :: dummyval=0.d0, xmin=0.d0, ymin=0.d0
      real*8, parameter :: xmax = 50, ymax = 25, omegacdmh2 = 0.12038d0 !From best fit Planck+WMAP (Planck 2013)  
      integer, parameter :: spoke_1_pts = 10, spoke_2_pts = 10   
      logical, parameter :: talky = .false.                      

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
      j=0

      ! Open a file with some SUSY models inside and loop over them.  See DarkSUSY's dstest for details.
      open (unit=11,file='programs/testmodels25.mod')
      ! Skip the header
      call read_model(11,-15,iend,ierr,0.d0,0.d0)

      do while (.not. finished)

        ! Switch depending on whether we're working on benchmarks or slopes
        if (on_spoke_1) then

          x = xmin + dble(j)/dble(spoke_1_pts) * (xmax - xmin)
          y = 0.d0
          if (talky) write (*,*) 'Spoke 1 of MODEL: ',idtag,' with x = ',x,' GeV.'
          call read_model(11,0,iend,ierr,x,y)
          if (j .ne. spoke_1_pts) then
            j = j + 1
          else 
            on_spoke_1 = .false.
            on_spoke_2 = .true.
            j = 1
          endif          
                    
        elseif (on_spoke_2) then
        
          x = 0.d0
          y = ymin + dble(j)/dble(spoke_2_pts) * (ymax - ymin)
          if (talky) write (*,*) 'Spoke 2 of MODEL: ',idtag,' with y = ',y,' GeV.'
          call read_model(11,0,iend,ierr,x,y)
          if (j .ne. spoke_2_pts) then
            j = j + 1
          else 
            finished = .true.
          endif          
        
        else

          ! Read in the model        
          call read_model(11,0,iend,ierr,0.d0,0.d0)
          if (ierr.ne.0) then
            write (*,*) 'Error while reading testmodels25.mod'
            write (*,*) 'The nulike MSSM-25 calculations cannot continue.'
            stop
          else if (iend.eq.1) then ! Reached the end of the model file; start on the spokes
            on_spoke_1 = .true.          
            j = 1
            cycle
          else
            if (talky) write(*,*)
            if (talky) write(*,*) 'MODEL: ',idtag,' read from file testmodels25.mod'
          endif

        endif

        ! Calculate the sparticle spectrum.  See DarkSUSY's dstest for details.
        call dssusy(unphys,hwarning)

        ! Write out neutralino mass and sigmav
        if (talky) write(*,*) '  Neutralino mass = ', wamwimp
        if (talky) write(*,*) '  Annihilation cross-section = ',wasv,' cm^-3 s^-1'

        ! Set the estimated relative theoretical error in neutrino flux calculation 
        theoryError = 5d-2 * merge(1.d0, dsqrt(wamwimp*1d-2), wamwimp .le. 100.d0)
        if (talky) write(*,*) '  Theory error = ',theoryError,'%'

        ! Calculate relic density.  See DarkSUSY's dstest for details.
        if (talky) write(*,*) 'Calculating omega h^2 with coannihilations...'
        oh2=dsrdomega(1,1,xf,ierr,iwar,nfc)
        if (talky) write(*,*) '  with coannihilations Oh2 = ',oh2

        ! Rescale the halo density according to the relic density.  See DarkSUSY's dstest for details.
        call dshmrescale_rho(oh2,omegacdmh2)

        ! Calculate the nuclear scattering cross-sections
        if (talky) write (*,*) 'Calculating scattering cross sections...'
        call dsddneunuc(sigsip,sigsin,sigsdp,sigsdn)
        if (talky) then
          write(*,*) '  sigsip (pb) = ',sigsip*1.0d36
          write(*,*) '  sigsin (pb) = ',sigsin*1.0d36
          write(*,*) '  sigsdp (pb) = ',sigsdp*1.0d36
          write(*,*) '  sigsdn (pb) = ',sigsdn*1.0d36
        endif

        ! Calculate capture and annihlation rates with DarkSUSY, using the
        ! default Gouldian calculation with a numerical integration over
        ! the velocity distribution, sped up by tabulation.
        if (talky) write(*,*) 'Calculating capture rate in the Sun'
        csu=dsntcapsuntab(wamwimp,sigsip,sigsdp) ! Capture rate (s^-1) (wamwimp = WIMP mass in GeV)
        ca=wasv/6.6d28*(wamwimp/20.d0)**(3./2.) 
        tausu=1.0d0/dsqrt(csu*ca)                ! Equilibration time (s)
        tt_sun=1.5d17*dsqrt(csu*ca)              
        annrate=csu*0.5d0*tanh(tt_sun)**2        ! Annihiliation rate (s^-1)
        if (talky) write(*,*) '  Capture rate in the Sun = ',csu,' s^-1'
        if (talky) write(*,*) '  Annihilation rate in the Sun = ',annrate,' s^-1'

        ! See the header of src/nulike_bounds.f for detailed explanations of the following options.
        speed = 0
        pvalFromRef = .true.
        dof = 1.0  !Calculate exclusion assuming conditioning on everything except a single parameter (eg mchi).
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

        if (talky) write(*,*) 'Calculating likelihood and p value.'
        do i = 1, 3
          !Use nulike to get signal and background predictions, number of observed events, likelihood and p-value
          call nulike_bounds(experiment(i), wamwimp, annrate, nuyield_test, sigpred(i), bgpred(i), 
     &     totobs(i), lnLike(i), pval(i), likechoice, theoryError, speed, pvalFromRef,
     &     refLike(i), dof, ptr, threadsafe)

          if (talky) then
            write(*,*) ' Individual results from ',trim(experiment(i))
            write(*,*) '  Predicted signal events:    ', sigpred(i)
            write(*,*) '  Total predicted events:     ', sigpred(i)+bgpred(i)
            write(*,*) '  Observed events:            ', totobs(i)
            write(*,*) '  log-likelihood:             ', lnLike(i)
            write(*,*) '  p-value:                    ', pval(i)
            write(*,*) '  Model ruled out at ', 1.d2-1.d2*pval(i),'% CL.'
          endif

        enddo
         
        excl = 1.d2-1.d2*DGAMIC(dof*0.5d0,max(0.d0, sum(refLike-lnLike)))/DGAMMA(dof*0.5d0)
        if (talky) then
          write(*,*) ' Combined results from WH, WL and SL datasets:'
          write(*,*) '  Predicted signal events:    ', sum(sigpred)
          write(*,*) '  Total predicted events:     ', sum(sigpred+bgpred)
          write(*,*) '  Observed events:            ', sum(totobs)
          write(*,*) '  log-likelihood:             ', sum(lnLike)
          write(*,*) '  p-value:                    ', minval(pval)
          write(*,*) '  Model ruled out at '
          write(*,*) '   by worst p value: ',  1.d2-1.d2*minval(pval),'% CL.'
          write(*,*) '   by composite likelihood: ', excl,'% CL.'
        else
          if (first) then
            first = .false.
            write(*,*) '#Model      x (GeV)        y (GeV)       ',
     &            ' mchi (GeV)     sigmav (cm^-3 s^-1) oh2      ',
     &            ' sig_SIp (cm^2) sig_SIn (cm^2) sig_SDp (cm^2)',
     &            ' sig_SDn (cm^2) sig_SIp`(cm^2) sig_SIn`(cm^2)',
     &            ' sig_SDp`(cm^2) sig_SDn`(cm^2) C (s^-1)      ',
     &            ' A (s^-1)       exclusion (%CL)'
          endif
          rs = min(oh2/omegacdmh2, 1.0)
          write(*,'(a12,16(1x,e14.8))') idtag, x, y, wamwimp, wasv, oh2, 
     &     sigsip, sigsin, sigsdp, sigsdn, rs*sigsip, rs*sigsin, rs*sigsdp,
     &     rs*sigsdn, csu, annrate, excl
        endif
        
      end do

      close (11)
      write(*,*)
      write(*,*) 'The nulike MSSM-25 benchmark tests completed successfully.'     
      write(*,*)

      end program nulike_test_mssm25


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


      ! Reads a model from a file into DarkSUSY
      subroutine read_model(lunit,nmodel,iend,ierr,x,y)
      implicit none
      include 'dsidtag.h'
      include 'dsmssm.h'
      integer nmodel,lunit,iend,ierr,i
      real*8 am1,am2,am3,amu,ama,atanbe,
     &  amsqL1,amsqL2,amsqL3,amsqRu,amsqRc,amsqRt,amsqRd,
     &  amsqRs,amsqRb,amslL1,amslL2,amslL3,amslRe,
     &  amslRmu,amslRtau,atm,abm,ataum,aemum,x,y
      real*8 am1_s,am2_s,am3_s,amu_s,ama_s,atanbe_s,
     &  amsqL1_s,amsqL2_s,amsqL3_s,amsqRu_s,amsqRc_s,amsqRt_s,amsqRd_s,
     &  amsqRs_s,amsqRb_s,amslL1_s,amslL2_s,amslL3_s,amslRe_s,
     &  amslRmu_s,amslRtau_s,atm_s,abm_s,ataum_s,aemum_s
      save am1_s,am2_s,am3_s,amu_s,ama_s,atanbe_s,
     &  amsqL1_s,amsqL2_s,amsqL3_s,amsqRu_s,amsqRc_s,amsqRt_s,amsqRd_s,
     &  amsqRs_s,amsqRb_s,amslL1_s,amslL2_s,amslL3_s,amslRe_s,
     &  amslRmu_s,amslRtau_s,atm_s,abm_s,ataum_s,aemum_s

 2000 format (1x,a12,25(1x,e14.8))
      ierr=0

      ! If x or y are non-zero, use them to rescale the previous point
      if (x .gt. epsilon(0.d0) .or. y .gt. epsilon(0.d0)) then

        am1 = am1_s - 12.d0*x - 16.d0*y
        am2 = am2_s + 64.d0*y
        am3 = am3_s + 11.d0*x + 64.d0*y
        amu = amu_s - 64*y
        ama = ama_s + 64*y
        atanbe = atanbe_s
        amsqL1 = amsqL1_s + 11*x + 9*y
        amsqL2 = amsqL2_s + 11*x + 9*y
        amsqL3 = amsqL3_s + 11*x + 64*y
        amsqRu = amsqRu_s + 11*x + 9*y
        amsqRc = amsqRc_s + 11*x + 9*y
        amsqRt = amsqRt_s + 64*y
        amsqRd = amsqRd_s + 64*y
        amsqRs = amsqRs_s + 64*y
        amsqRb = amsqRb_s + 64*y
        amslL1 = amslL1_s + 11*x + 64*y
        amslL2 = amslL2_s + 11*x + 64*y
        amslL3 = amslL3_s + 64*y
        amslRe = amslRe_s + 64*y
        amslRmu = amslRmu_s + 64*y
        amslRtau = amslRtau_s + 64*y
        atm = atm_s - 64*y
        abm = abm_s + 64*y
        ataum = ataum_s + 64*y
        aemum = aemum_s

      else 

        ! When nmodel<0, skip -nmodel lines
        if (nmodel.lt.0) then
           do i=1,-nmodel
              read (lunit,*,end=1000,err=3000)
           enddo
           return
        endif
        ! If nmodel>0, read n-th model (assumes no headers)
        if (nmodel.gt.0) then
           do i=1,nmodel-1
              read (lunit,*,end=1000,err=3000)
           enddo
        endif
        ! If nmodel=0, read next model
        read (lunit,2000,end=1000,err=1000) 
     &   idtag, am1,am2,am3,amu,ama,atanbe,
     &   amsqL1,amsqL2,amsqL3,amsqRu,amsqRc,amsqRt,amsqRd,
     &   amsqRs,amsqRb,amslL1,amslL2,amslL3,amslRe,
     &   amslRmu,amslRtau,atm,abm,ataum,aemum
        ! Save for next time, in case this is the last.
        am1_s = am1
        am2_s = am2
        am3_s = am3
        amu_s = amu
        ama_s = ama
        atanbe_s = atanbe
        amsqL1_s = amsqL1
        amsqL2_s = amsqL2
        amsqL3_s = amsqL3
        amsqRu_s = amsqRu
        amsqRc_s = amsqRc
        amsqRt_s = amsqRt
        amsqRd_s = amsqRd
        amsqRs_s = amsqRs
        amsqRb_s = amsqRb
        amslL1_s = amslL1
        amslL2_s = amslL2
        amslL3_s = amslL3
        amslRe_s = amslRe
        amslRmu_s = amslRmu
        amslRtau_s = amslRtau
        atm_s = atm
        abm_s =  abm
        ataum_s =  ataum
        aemum_s = aemum
        ! modify/set additional parameters
        !  higloop=5  ! 5 = Full FeynHiggs;  6 = FeynHiggsFast
        !  W mass for unitarity of tree-level annihilation amplitudes
      endif
      
      call dsgive_model25(am1,am2,am3,amu,ama,atanbe,
     &  amsqL1,amsqL2,amsqL3,amsqRu,amsqRc,amsqRt,amsqRd,
     &  amsqRs,amsqRb,amslL1,amslL2,amslL3,amslRe,
     &  amslRmu,amslRtau,atm,abm,ataum,aemum)

      return

 1000 continue
      iend=1
      write (*,*) 'End of model file reached.'
      return 

 3000 continue
      ierr=1
      write(*,*) 'read_model threw an error...'
      stop

      end
