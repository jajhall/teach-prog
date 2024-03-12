      program picard
      implicit none
      integer mx_k
      integer k

      double precision x0, x1, x2, x3, x4, x5, x6, phi
      data x0/-0.01/, x1/0.0/, x2/0.01/, x3/0.5/
      data x4/0.99/, x5/1.0/, x6/1.01/
      data mx_k/21/

      do 10, k=1, mx_k
         write(1,9000)k-1, x0, x1, x2, x3, x4, x5, x6
         x0 = phi(x0)
         x1 = phi(x1)
         x2 = phi(x2)
         x3 = phi(x3)
         x4 = phi(x4)
         x5 = phi(x5)
         x6 = phi(x6)
 10   continue
 9000 format(i2, 1x, f11.4, 6(1x,f10.7))
      end

      double precision function phi(x)
      implicit none
      double precision x

      phi = (-x+2d0)*x
      return
      end
