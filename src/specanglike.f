***********************************************************************
*** nulike_specanglike returns the contribution of a single event to the 
*** unbinned likelihood, based on values of the arrival direction and 
*** the energy estimator for the event, as well as the theoretical 
*** energy spectrum of incoming neutrinos.
*** This routine is used only with the 2014 likelihood.
***
*** input:  event_number
***         theta_S      total predicted number of signal events
***                       within analysis window (cut cone)
***         f_S          signal fraction; percentage of predicted counts
***                       expected to be due to signal rather than back-
***                       ground.
***         annrate      Annihilation rate (s^-1) 
***         logmw        log_10(m_WIMP / GeV)
***         logEmin      log_10(Emin/GeV), where Emin is the lower energy
***                       boundary of the analysis energy range  
***         nuyield      external double function that returns
***                       the differential neutrino flux
***                       at the detector in units of m^-2 GeV^-1 
***                       annihilation^-1
***        context       A c_ptr passed in to nuyield when it is called
***
*** output:              ln(Likelihood / chan^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Jun 8, 2014
***********************************************************************


      double precision function nulike_specanglike(event_number,theta_S,
     & f_S,annrate,logmw,logEmin,nuyield,context)

      use iso_c_binding, only: c_ptr
      use Precision_Model
      use CUI

      implicit none
      include 'nulike_internal.h'

      integer event_number
      real*8 theta_S, f_S, annrate, logmw, logEmin, nuyield
      real*8 bgpartial, sigpartial, cosphi, ee, eps
      real*8 nulike_bgangpdf, nulike_bgspec, SAbsErr, SVertices(1,2)
      integer IER, SRgType
      type(c_ptr) context
      external nuyield, nulike_specangintegrand
      parameter (eps = 5.d-3, SRgType = HyperQuad)

      interface
        function nulike_specangintegrand(NumFun,X) result(Value)
          use iso_c_binding, only: c_ptr
          include 'nulike_internal.h'
          integer, intent(in) :: NumFun
          real*8, intent(in) :: X(:)
          real*8 :: Value(NumFun)
        end function nulike_specangintegrand
      end interface

      ! Retrieve the event info
      cosphi = events_cosphi(event_number,analysis)
      ee = events_ee(event_number,analysis)

      ! Include the background component
      bgpartial = nulike_bgangpdf(cosphi) * nulike_bgspec(ee,2014)

      ! Include the signal component
      if (logEmin .gt. logmw) then
        sigpartial = 0.d0
      else
        eventnumshare = event_number
        IER = 0
        SVertices(1,:) = (/logEmin, logmw/)
        call CUBATR(1,nulike_specangintegrand,SVertices,SRgType,
     &   sigpartial,SAbsErr,IER,MaxPts=5000000,EpsRel=eps,Job=2,Key=2)
        if (IER .ne. 0) then
          write(*,*) 'Error raised by CUBATR in nulike_specanglike: ', IER 
          stop
        endif
        call CUBATR()
        sigpartial = exp_time(analysis) / theta_S * annrate * dlog(10.d0) * sigpartial
      endif

      ! Combine background and signal components
      nulike_specanglike = f_S * sigpartial + (1.d0-f_S) * bgpartial

      ! Take ln of total likelihood
      nulike_specanglike = log(nulike_specanglike)

      end function nulike_specanglike
