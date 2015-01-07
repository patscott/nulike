!=======================================================================
!  Integrate function f between a and b using Simpson's rule.
!
!  Input:  integrand f, auxilary function aux
!          integration limits a and b
!
!  Author: Joakim Edsjo (edsjo@physto.se) 96-05-16
!          2000-07-19 Paolo Gondolo added eps as argument
!          2007-07-27 Pat Scott added j > 6 condition to prevent
!           spurious early convergence, removed os initial definition.
!          2014-03-05 Pat Scott removed from DarkSUSY and added aux.
!
!  Adapted from DarkSUSY r919 (v5.1.1+)
!  Initially based on Paolo Gondolo's wxint.f routine.
!=======================================================================

      real*8 function nulike_simpson1(f,aux,a,b,eps)
      implicit none

      real*8 f,aux,a,b,tot,eps,st,os,ost,del,sum,x
      integer jmax,it,l,j,nfcn,jdid
      external f,aux
      parameter (jmax=25)

      nulike_simpson1=0.d0
      del=b-a
      ost=0.5*del*(f(a,aux)+f(b,aux))
      x=0.5*(b+a)
      st=0.5*(ost+del*f(x,aux))
      ost=st
      it=1
      nfcn=3
      do j=3,jmax
        it=2*it
        del=0.5*del
        x=a+0.5*del
        sum=0.0
        do l=1,it
          sum=sum+f(x,aux)
          nfcn=nfcn+1
          x=x+del
        enddo
        st=0.5*(st+del*sum)
        tot=(4.0*st-ost)/3.0
        jdid=j
        if (j.gt.6) then
           if (abs(tot-os).le.eps*abs(os)) then
              nulike_simpson1=tot
              return
           endif
        endif
        os=tot
        ost=st
      enddo
      write(*,*) '  Error in nulike_simpson1: too many steps.'
      write(*,*) '  Integral set to zero.'
      nulike_simpson1=0.0d0

      end
