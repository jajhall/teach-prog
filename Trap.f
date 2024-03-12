      program numint
      implicit none
      integer n
      double precision a, b, h, trap, simp, er, tru_int
      double precision xj
      double precision f
      integer alg, j

      print*, ' Enter Trap/Simp 0/1'
      read*, alg
      print*, ' Enter I'
      read*, tru_int
      if (alg .eq. 0) then
         write(*, 9000)
      else
         write(*, 9010)
      endif
 10   continue
      print*, ' Enter n'
      read*, n
      if (n .le. 0) stop
      a = 0.0
      b = 1.0
      h = (b-a)/n
      
      if (alg .eq. 0) then
         trap = 0d0
         do 100, j = 1, n-1
            xj = a + j*h
            trap = trap + f(xj)
 100     continue
         trap = trap*2.0
         trap = trap + f(a)
         trap = trap + f(b)
         trap = trap*(h/2.0)
         er = tru_int - trap
         write(*, 9100)h, tru_int, trap, er, er/h, er/h**2, er/h**3
      else
         write(*, 9010)
      endif
 9000 format('h I T(h) E(h) E(h)/h E(h)/h^2 E(h)/h^3')
 9010 format('h I S(h) E(h) E(h)/h E(h)/h^2 E(h)/h^3 E(h)/h^4 E(h)/h^5')
 9100 format(10(1x, g15.8))
      goto 10
      end
      
      double precision function f(x)
      implicit none
      double precision x

      f = 1.0/(1.0+x*x)
c      f = sqrt(1.0-x*x)
      return
      end
