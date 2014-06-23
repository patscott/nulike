
      real*8 function d1mach(i)
      integer i
c
c  double-precision machine constants
c  d1mach( 1) = b**(emin-1), the smallest positive magnitude.
c  d1mach( 2) = b**emax*(1 - b**(-t)), the largest magnitude.
c  d1mach( 3) = b**(-t), the smallest relative spacing.
c  d1mach( 4) = b**(1-t), the largest relative spacing.
c  d1mach( 5) = log10(b)
c
      integer small(2)
      integer large(2)
      integer right(2)
      integer diver(2)
      integer log10(2)
      integer sc, cray1(38), j
      common /d9mach/ cray1
      save small, large, right, diver, log10, sc
      real*8 dmach(5)
      equivalence (dmach(1),small(1))
      equivalence (dmach(2),large(1))
      equivalence (dmach(3),right(1))
      equivalence (dmach(4),diver(1))
      equivalence (dmach(5),log10(1))
c  this version adapts automatically to most current machines.
c  r1mach can handle auto-double compiling, but this version of
c  d1mach does not, because we do not have quad constants for
c  many machines yet.
c  to compile on older machines, add a c in column 1
c  on the next line
      data sc/0/
c  and remove the c from column 1 in one of the sections below.
c  constants for even older machines can be obtained by
c          mail netlib@research.bell-labs.com
c          send old1mach from blas
c  please send corrections to dmg or ehg@bell-labs.com.
c
c     machine constants for the honeywell dps 8/70 series.
c      data small(1),small(2) / o402400000000, o000000000000 /
c      data large(1),large(2) / o376777777777, o777777777777 /
c      data right(1),right(2) / o604400000000, o000000000000 /
c      data diver(1),diver(2) / o606400000000, o000000000000 /
c      data log10(1),log10(2) / o776464202324, o117571775714 /, sc/987/
c
c     machine constants for pdp-11 fortrans supporting
c     32-bit integers.
c      data small(1),small(2) /    8388608,           0 /
c      data large(1),large(2) / 2147483647,          -1 /
c      data right(1),right(2) /  612368384,           0 /
c      data diver(1),diver(2) /  620756992,           0 /
c      data log10(1),log10(2) / 1067065498, -2063872008 /, sc/987/
c
c     machine constants for the univac 1100 series.
c      data small(1),small(2) / o000040000000, o000000000000 /
c      data large(1),large(2) / o377777777777, o777777777777 /
c      data right(1),right(2) / o170540000000, o000000000000 /
c      data diver(1),diver(2) / o170640000000, o000000000000 /
c      data log10(1),log10(2) / o177746420232, o411757177572 /, sc/987/
c
c     on first call, if no data uncommented, test machine types.
      if (sc .ne. 987) then
         dmach(1) = 1.d13
         if (      small(1) .eq. 1117925532
     *       .and. small(2) .eq. -448790528) then
*           *** ieee big endian ***
            small(1) = 1048576
            small(2) = 0
            large(1) = 2146435071
            large(2) = -1
            right(1) = 1017118720
            right(2) = 0
            diver(1) = 1018167296
            diver(2) = 0
            log10(1) = 1070810131
            log10(2) = 1352628735
         else if ( small(2) .eq. 1117925532
     *       .and. small(1) .eq. -448790528) then
*           *** ieee little endian ***
            small(2) = 1048576
            small(1) = 0
            large(2) = 2146435071
            large(1) = -1
            right(2) = 1017118720
            right(1) = 0
            diver(2) = 1018167296
            diver(1) = 0
            log10(2) = 1070810131
            log10(1) = 1352628735
         else if ( small(1) .eq. -2065213935
     *       .and. small(2) .eq. 10752) then
*               *** vax with d_floating ***
            small(1) = 128
            small(2) = 0
            large(1) = -32769
            large(2) = -1
            right(1) = 9344
            right(2) = 0
            diver(1) = 9472
            diver(2) = 0
            log10(1) = 546979738
            log10(2) = -805796613
         else if ( small(1) .eq. 1267827943
     *       .and. small(2) .eq. 704643072) then
*               *** ibm mainframe ***
            small(1) = 1048576
            small(2) = 0
            large(1) = 2147483647
            large(2) = -1
            right(1) = 856686592
            right(2) = 0
            diver(1) = 873463808
            diver(2) = 0
            log10(1) = 1091781651
            log10(2) = 1352628735
         else if ( small(1) .eq. 1120022684
     *       .and. small(2) .eq. -448790528) then
*           *** convex c-1 ***
            small(1) = 1048576
            small(2) = 0
            large(1) = 2147483647
            large(2) = -1
            right(1) = 1019215872
            right(2) = 0
            diver(1) = 1020264448
            diver(2) = 0
            log10(1) = 1072907283
            log10(2) = 1352628735
         else if ( small(1) .eq. 815547074
     *       .and. small(2) .eq. 58688) then
*           *** vax g-floating ***
            small(1) = 16
            small(2) = 0
            large(1) = -32769
            large(2) = -1
            right(1) = 15552
            right(2) = 0
            diver(1) = 15568
            diver(2) = 0
            log10(1) = 1142112243
            log10(2) = 2046775455
         else
            dmach(2) = 1.d27 + 1
            dmach(3) = 1.d27
            large(2) = large(2) - right(2)
            if (large(2) .eq. 64 .and. small(2) .eq. 0) then
               cray1(1) = 67291416
               do 10 j = 1, 20
                  cray1(j+1) = cray1(j) + cray1(j)
 10               continue
               cray1(22) = cray1(21) + 321322
               do 20 j = 22, 37
                  cray1(j+1) = cray1(j) + cray1(j)
 20               continue
               if (cray1(38) .eq. small(1)) then
*                  *** cray ***
                  call i1mcry(small(1), j, 8285, 8388608, 0)
                  small(2) = 0
                  call i1mcry(large(1), j, 24574, 16777215, 16777215)
                  call i1mcry(large(2), j, 0, 16777215, 16777214)
                  call i1mcry(right(1), j, 16291, 8388608, 0)
                  right(2) = 0
                  call i1mcry(diver(1), j, 16292, 8388608, 0)
                  diver(2) = 0
                  call i1mcry(log10(1), j, 16383, 10100890, 8715215)
                  call i1mcry(log10(2), j, 0, 16226447, 9001388)
               else
                  write(*,9000)
                  stop 779
                  end if
            else
               write(*,9000)
               stop 779
               end if
            end if
         sc = 987
         end if
*    sanity check
      if (dmach(4) .ge. 1.0d0) stop 778
      if (i .lt. 1 .or. i .gt. 5) then
         write(*,*) 'd1mach(i): i =',i,' is out of bounds.'
         stop
         end if
      d1mach = dmach(i)
      return
 9000 format(/' adjust d1mach by uncommenting data statements'/
     *' appropriate for your machine.')
* /* standard c source for d1mach -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*double d1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return dbl_min;
*	  case 2: return dbl_max;
*	  case 3: return dbl_epsilon/flt_radix;
*	  case 4: return dbl_epsilon;
*	  case 5: return log10((double)flt_radix);
*	  }
*	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
*	exit(1); return 0; /* some compilers demand return values */
*}
      end
      subroutine i1mcry(a, a1, b, c, d)
**** special computation for old cray machines ****
      integer a, a1, b, c, d
      a1 = 16777216*b + c
      a = 16777216*a1 + d
      end
