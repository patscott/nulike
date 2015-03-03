************************************************************************
      real*8 function expbesi0(t,u,w)
************************************************************************
c     exp[-(t^2 + u^2)/2w] * I_0(tu/w)
      implicit none
      real*8 t,u,w,x,y,z,ax
      double precision p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      data p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     *     1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     *     0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,
     *     0.2635537d-1,-0.1647633d-1,0.392377d-2/
      x = t*u/w
      ax=abs(x)
      if (ax.lt.3.75) then
         y=(x/3.75)**2
         expbesi0=(p1+y*(p2+y*(p3+y*(p4+y*
     *             (p5+y*(p6+y*p7)))))) * exp(-0.5d0*(t*t+u*u)/w)
      else
         y=3.75/ax
         z=t-u
         expbesi0=exp(-0.5d0*z*z/w)*(q1+y*(q2+y*(q3+y*(q4+y*
     *            (q5+y*(q6+y*(q7+y*(q8+y*q9))))))))/sqrt(ax)
      endif
      return
      end
