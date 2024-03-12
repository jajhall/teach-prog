      program picard
      implicit none
      integer mx_k
      integer k
      double precision phi0, phi1, phi2
      double precision x, x0, x1, x2
      double precision r1, r2
      data  x0/1.0/x1/1.0/x2/1.0/
      data mx_k/21/

      x = 1.3247179572447460260
      r1 = 1d0
      r2 = 1d0
      do 10, k=1, mx_k
         write(1,9000)k-1, x1, (x1-x), (x1-x)/r1, x2, (x2-x), (x2-x)/r2
         r1 = x1-x
         r2 = x2-x
         x1 = phi1(x1)
         x2 = phi2(x2)
 10   continue
 9000 format(i2, 2(1x, f10.7, 2(e11.4)))
      end

      double precision function phi0(x)
      implicit none
      double precision x
      phi0 = x**3-1
      return
      end
      double precision function phi1(x)
      implicit none
      double precision x
      double precision power
      power = 1d0/3d0
      phi1 = (x+1)**power
      return
      end
      double precision function phi2(x)
      implicit none
      double precision x
      phi2 = (2d0*x*x*x+1d0)/(3d0*x*x-1d0)
      return
      end

