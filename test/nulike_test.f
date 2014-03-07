! This program tests the nulike routines, and is adapted from dstest in 
! DarkSUSY.
!
! Author: Pat Scott patscott@phsyics.mcgill.ca
! Date: March 6 2014 


      program nulike_test

      implicit none
      real*8 oh2,xf,dsrdomega                                 ! relic density
      real*8 sigsip,sigsin,sigsdp,sigsdn                      ! nuclear scattering
      real*8 dsntcapsun, dsntcapsunnum, dsntcapsuntab, ca     ! capture rate 
      real*8 tt_sun, annrate, nuyield                         ! capture rate
      real*8 eth,thmax,rateea,ratesu,energy,theta             ! neutrino telescopes
      integer rtype,ptype                                     ! neutrino telescopes
      real*8 sigpred, bgpred, lnLike, pval, refLike, dof      ! neutrino likelihood
      real*8 theoryError,phi_cut                              ! neutrino likelihood
      integer totobs, likechoice                              ! neutrino likelihood
      logical uselogNorm, pvalFromRef                         ! neutrino likelihood
      logical BGLikePrecompute                                ! neutrino likelihood
      character (len=256) iclike2012, iclike2014, experiment  ! neutrino likelihood
      character (len=256) eventf, edispf, BGf, efaorvf        ! neutrino likelihood
      integer unphys,hwarning,iend,ierr,istat,iwar,nfc        ! bookkeeping
      external nuyield

      include 'dsio.h'
      include 'dsidtag.h'
      include 'dsntcom.h'
      include 'dsdirver.h'
      include 'dswacom.h'

      ! Initialise DarkSUSY
      call dsinit

      ! Choose a print level for DarkSUSY output to stdout
      prtlevel=1

      ! Before starting anything, we initialise the neutrino telescope data and
      ! the nulike routines.  See the header of src/init.f for detailed 
      ! explanations of the following options.

      ! Choose the data files containing neutrino telescope events, energy dispersion,
      ! observed background distributions and effective area/median angular 
      ! resolution.  

      ! The likelihood2012 folder contains data files designed for use with the
      ! 2012 likelihood (Scott, Savage, Edsjö & IceCube Collab 2012, 
      ! JCAP 11:057, arXiv:1207.0810).
      ! The likelihood2014 folder contains data files designed for use with the
      ! 2014 likelihood (IceCube Collab 2014, JCAP xx:xxx, arXiv:14xx.xxxx).
      iclike2012 = 'data/IceCube/likelihood2012/'
      iclike2014 = 'data/IceCube/likelihood2014/'

      ! Here we use the IC-22 data that ships with nulike.
      !  experiment = 'IC-22'
      !  eventf  = trim(iclike2012)//'events_10deg_IC22.dat'
      !  edispf  = trim(iclike2012)//'energy_histograms_IC22.dat'
      !  BGf     = trim(iclike2012)//'BG_distributions_IC22.dat'
      !  efareaf = trim(iclike2012)//'nuEffArea_IC22.dat'

      ! Here we use the IC79 WH data that ships with nulike.
      !  experiment = 'IC-79 WH'
      !  eventf  = trim(iclike2014)//'IC79_Events_WH_10degrees.dat'
      !  edispf  = trim(iclike2014)//'IC79_energy_histograms_WH.dat'
      !  BGf     = trim(iclike2014)//'IC79_Background_distributions_WH.dat'
      !  efaorvf = trim(iclike2014)//'IC79_Effective_Volume_WH.dat'

      ! Here we use the IC-86 simulation that ships with nulike
      experiment = 'IC-86 (predicted)'
      eventf  = trim(iclike2012)//'events_20deg_IC86_sim_nosig.dat'
      edispf  = trim(iclike2012)//'energy_histograms_IC86_sim_dummy.dat'
      BGf     = trim(iclike2012)//'BG_distributions_IC86_sim.dat'
      efaorvf = trim(iclike2012)//'nuEffArea_IC86_sim.dat'

      ! Set the analysis cut in degrees around the solar position
      ! phi_cut = 10.d0 !for IC22 or IC79 WH
      phi_cut = 20.d0 !for IC86

      ! Set the estimated relative theoretical error in neutrino flux calculation 
      theoryError = 0.05d0

      ! Choose a log-normal or a Gaussian distribution for the systematic error on
      ! the number of neutrino events
      uselogNorm = .true.

      ! Choose whether to precompute the background p-value in the case where Poissonian 
      ! statistics is used.
      BGLikePrecompute = .true.

      ! Initialise the IceCube data and calculations. 
      call nulike_init(eventf, edispf, BGf, efaorvf, phi_cut,
     & theoryError, uselogNorm, BGLikePrecompute)

      ! Open a file with some SUSY models inside.  See DarkSUSY's dstest for details.
      open (unit=11,file='test/testmodels.mod')
      read (11,*) ! this skips the header
      iend=0

      ! Loop over the points in the model file
      do

        ! Read in the model
        call read_model(11,0,iend,ierr)
        if (ierr.ne.0) then
          write (*,*) 'Error while reading testmodels.mod'
          write (*,*) 'The nulike test cannot continue.'
          stop
        else if (iend.eq.1) then ! Reached the end of the model file; exit the do loop
          exit
        else
          write(*,*)
          write(*,*) 'Model parameters read from file testmodels.mod'
        endif
        write(*,*) 
        write(*,*) '***** MODEL: ',idtag,' *****'

        ! Tone down the output from here on.
        prtlevel=0

        ! Calculate the sparticle spectrum.  See DarkSUSY's dstest for details.
        call dssusy(unphys,hwarning)

        ! Calculate relic density.  See DarkSUSY's dstest for details.
        write(*,*) 'Calculating omega h^2 with coannihilations...'
        oh2=dsrdomega(1,1,xf,ierr,iwar,nfc)
        write(*,*) '  with coannihilations Oh2 = ',oh2

        ! Rescale the halo density according to the relic density.  See DarkSUSY's dstest for details.
        call dshmrescale_rho(oh2,0.025d0)

        ! Calculate the nuclear scattering cross-sections
        write (*,*) 'Calculating scattering cross sections...'
        call dsddneunuc(sigsip,sigsin,sigsdp,sigsdn)
        write(*,*) '  sigsip (pb) = ',sigsip*1.0d36
        write(*,*) '  sigsin (pb) = ',sigsin*1.0d36
        write(*,*) '  sigsdp (pb) = ',sigsdp*1.0d36
        write(*,*) '  sigsdn (pb) = ',sigsdn*1.0d36

        ! Calculate capture and annihlation rates with DarkSUSY, using the
        ! default Gouldian calculation with a numerical integration over
        ! the velocity distribution, sped up by tabulation.
        write(*,*) 'Calculating capture rate in the Sun'
        csu=dsntcapsuntab(wamwimp,sigsip,sigsdp) ! Capture rate (s^-1) (wamwimp = WIMP mass in GeV)
        ca=wasv/6.6d28*(wamwimp/20.d0)**(3./2.) 
        tausu=1.0d0/dsqrt(csu*ca)                ! Equilibration time (s)
        tt_sun=1.5d17*dsqrt(csu*ca)              
        annrate=csu*0.5d0*tanh(tt_sun)**2        ! Annihiliation rate (s^-1)
        write(*,*) '  Capture rate in the Sun = ',csu,' s^-1'
        write(*,*) '  Annihilation rate in the Sun = ',annrate,' s^-1'


        !  Likelihoods and p values from IceCube
        write(*,*) 'Calculating likelihood and p value for '//trim(experiment)

        !  See the header of src/nulike_bounds.f for
        !  more detailed explanations of the following options.

        ! Choose what to include in the actual likelihood calculations.  Note that this
        ! only applies to the likelihood calculation; the p value is always calculated 
        ! considering only the number count, not the spectral or angular information.
        ! If you're only interested in p-values, not likelihoods, use likechoice = 1, as
        ! this is the fastest.
        !likechoice = 1		! Number of events only
        !likechoice = 2		! Number of events and event arrival angles
        !likechoice = 3		! Number of events and energy estimator (nchan = number of lit DOMs)
        likechoice = 4		! Number of events, event arrival angles and energy estimator

        ! Choose whether to calculate the p value relative to a reference value of 
        ! the likelihood or to the background
        pvalFromRef = .false.

        ! Set the reference value of ln(Likelihood).  (Ignored if pvalFromRef = F)
        refLike = -500.d0

        ! Set the number of degrees of freedom to use in the p-value calculation 
        ! (Ignored if pvalFromRef = F).
        dof = 8.d0
   
        ! Finally use nulike to get signal and background predictions, number of observed events, likelihood and p-value
        call nulike_bounds(wamwimp, annrate, nuyield, sigpred, bgpred, 
     &   totobs, lnLike, pval, likechoice, pvalFromRef, refLike, dof)
     
        write(*,*) '  Predicted signal events:    ', sigpred
        write(*,*) '  Total predicted events:     ', sigpred+bgpred
        write(*,*) '  Observed events:            ', totobs
        write(*,*) '  log-likelihood:             ', lnLike
        write(*,*) '  p-value:                    ', pval
        write(*,*) '  Model ruled out at ', 1.d2-1.d2*pval,'% CL.'
        
      end do

      close (11)
      write(*,*)
      write(*,*) 'The nulike test program completed successfully.'     
      write(*,*)

      end program nulike_test


      ! Function returning neutrino flux at detector for 2012 likelihood.
      real*8 function nuyield(log10E,ptype)
      implicit none
      real*8 log10E, dsntmuonyield
      integer ptype, istat
      external dsntmuonyield
      nuyield = 1.d-30 * dsntmuonyield(10.d0**log10E,10.d0,'su',3,1,ptype,istat)
      end


      ! Reads a model from a file into DarkSUSY
      subroutine read_model(lunit,nmodel,iend,ierr)
      implicit none
      include 'dsidtag.h'
      include 'dsmssm.h'
      integer nmodel,lunit,iend,ierr
      real*8 at,ab,mqtild
      integer i
      character*40 message
 2000 format (1x,a12,7(1x,e14.8))
      ierr=0
      ! When nmodel<0, skip -nmodel lines
      if (nmodel.lt.0) then
         do i=1,-nmodel
            read (lunit,*,end=1000,err=3000)
         enddo
         return
      endif
      ! If nmodel>0, read n-th model (assumes header line)
      if (nmodel.gt.0) then
         do i=1,nmodel
            read (lunit,*,end=1000,err=3000)
         enddo
      endif
      ! If nmodel=0, read next model
      read (lunit,2000,end=1000,err=1000) 
     &     idtag,mu,m2,ma,tanbe,mqtild,at,ab
      ! modify/set additional parameters
      !  higloop=5  ! 5 = Full FeynHiggs;  6 = FeynHiggsFast
      !  W mass for unitarity of tree-level annihilation amplitudes
      call dsgive_model(mu,m2,ma,tanbe,mqtild,at,ab)
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
