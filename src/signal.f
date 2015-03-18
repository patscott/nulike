***********************************************************************
*** nulike_signal computes the predicted number of neutrino events due
*** to neutralino annihilation.
***        
*** Input:      nuyield         External double function that returns
***                             the differential neutrino flux
***                             at the detector in units of m^-2 GeV^-1
***                             annihilation^-1
***             context         A c_ptr passed in to nuyield when it is
***                             called
***             annrate         Annihilation rate (s^-1) 
***             logmw           log_10(m_WIMP / GeV)
***             like            Likelihood version (2012 or 2014)
***
*** Output:     theta_S         predicted number of signal events.
*** 
*** Author: Pat Scott (p.scott@imperial.ac.uk)
*** Date: Apr 22, 2011
*** Modified: March 6 2014
*** Modified: Jun 3, 6, 8 2014
***********************************************************************


      double precision function nulike_signal(nuyield, context, annrate, logmw, like)

      use iso_c_binding, only: c_ptr
      use Precision_Model
      use CUI

      implicit none
      include 'nulike_internal.h'

      real*8 integral, logmw, nuyield, upperLimit, theta_Snu, theta_Snubar, annrate
      real*8 eps2012, eps2014, SAbsErr, SVertices(1,2)
      integer like, IER, SRgType
      type(c_ptr) context
      parameter (eps2012 = 1.d-2, eps2014 = 3.d-2, SRgType = Simplex)
      external nuyield
      
      interface
        function nulike_sigintegrand(NumFun,X) result(Value)
          integer, intent(in) :: NumFun
          real*8, intent(in) :: X(:)
          real*8 :: Value(NumFun)
        end function nulike_sigintegrand
      end interface

      interface
        function nulike_specangintegrand(NumFun,X) result(Value)
          integer, intent(in) :: NumFun
          real*8, intent(in) :: X(:)
          real*8 :: Value(NumFun)
        end function nulike_specangintegrand
      end interface

      ! Short-circuit if the mass is too low to produce any observable events.
      if (logmw .lt. sens_logE(1,1,analysis)) then
        nulike_signal = 0.d0
        return
      endif
      
      ! Set the global context pointers unable to be passed through CUBPACK
      context_shared = context
      nuyield_ptr => nuyield

      ! Switch according to likelihood version.
      select case (like)

      ! 2012 likelihood, as per arXiv:1207.0810
      case (2012)

        if (logmw .lt. sens_logE(2,nSensBins(analysis),analysis)) then
          upperLimit = logmw
        else
          upperLimit = sens_logE(2,nSensBins(analysis),analysis)
        endif

        ! Neutrinos
        ptypeshare = 1
        IER = 0
        SVertices(1,:) = (/sens_logE(1,1,analysis), upperLimit/)
        call CUBATR(1,nulike_sigintegrand,SVertices,SRgType,
     &   integral,SAbsErr,IER,MaxPts=5000000,EpsRel=eps2012,Job=2,Key=2)
        if (IER .ne. 0) then
          write(*,*) 'Error raised by CUBATR in nulike_signal: ', IER 
          stop
        endif
        call CUBATR()
        theta_Snu = integral * dlog(10.d0) * exp_time(analysis) * annrate

        ! Anti-neutrinos
        ptypeshare = 2
        IER = 0
        SVertices(1,:) = (/sens_logE(1,1,analysis), upperLimit/)
        call CUBATR(1,nulike_sigintegrand,SVertices,SRgType,
     &   integral,SAbsErr,IER,MaxPts=5000000,EpsRel=eps2012,Job=2,Key=2)
        if (IER .ne. 0) then
          write(*,*) 'Error raised by CUBATR in nulike_signal: ', IER 
          stop
        endif
        call CUBATR()
        theta_Snubar = integral * dlog(10.d0) * exp_time(analysis) * annrate

        ! Total
        nulike_signal = theta_Snu + theta_Snubar 

      !2014 likelihood, as per arXiv:150x.xxxxx
      case (2014)

        eventnumshare = 0 ! Use effective area from previous tabulation.
        IER = 0
        SVertices(1,:) = (/precomp_log10E(1,analysis), 
     &   min(precomp_log10E(nPrecompE(analysis),analysis),logmw)/)
        call CUBATR(1,nulike_specangintegrand,SVertices,SRgType,
     &   integral,SAbsErr,IER,MaxPts=5000000,EpsRel=eps2014,Job=2,Key=2)
        if (IER .ne. 0) then
          write(*,*) 'Error raised by CUBATR in nulike_signal: ', IER 
          stop
        endif
        call CUBATR()
        nulike_signal = integral * dlog(10.d0) * exp_time(analysis) * annrate

      case default
        write(*,*) "Unrecognised likelihood version in nulike_signal."
        write(*,*) "Quitting..."
        stop

      end select


      end function nulike_signal

