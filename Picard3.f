      program picard
      implicit none
      integer mx_k
      integer k
      double precision phi
      double precision x, x1, r
      double precision err
      data  x1/0.1/
      data mx_k/21/

      read*, r
      x = 1d0/r
      err = 1d0
      do 10, k=1, mx_k
         write(1,9000)k-1, x1, abs(x1-x), abs(x1-x)/(err*err*err)
         err = abs(x1-x)
         x1 = phi(x1, r)
 10   continue
 9000 format(i2, 1x, f10.8, 2(1x, e11.4))
      end

      double precision function phi(x, r)
      implicit none
      double precision x, r
      phi = x*(3d0+r*x*(-3d0+r*x))
c      phi = x*(2d0-r*x)
      return
      end
