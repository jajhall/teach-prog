      program numint
      implicit none
      integer n_ivl
      double precision a, b, h, trap, simp, tru_int
      double precision tru_simp, int_simp_er, simp_er, int_tru_simp_er
      double precision tru_trap, int_trap_er, trap_er, int_tru_trap_er
      double precision m2, m4
      external f
      double precision f
      integer alg, og_n_dpl, n_dpl

      a = 0.0
      b = 1.0
      print*, ' Enter Trap/Simp 1/0'
      read*, alg
      print*, ' Enter I'
c      read*, tru_int
c      tru_int = atan(1.0)
      tru_int = 1d0-exp(-1d0)
      print*, ' Enter M2, M4'
      read*, M2, M4
      print*, ' Enter number of decimal places'
      read*, og_n_dpl
      if (og_n_dpl .lt. 0) then
         if (alg .eq. 1) then
            write(1, 9000)
         else
            write(1, 9001)
         endif
      else
         if (alg .eq. 1) then
            write(1, 9010)
         else
            write(1, 9011)
         endif
      endif
 10   continue
      print*, ' Enter number of intervals '
      read*, n_ivl
      if (n_ivl .le. 0) stop
      h = (b-a)/float(n_ivl)

      n_dpl = og_n_dpl
      if (alg .eq. 1) then
         call trap_int(a, b, n_ivl, f, n_dpl, trap)
         int_trap_er = tru_int - trap
         if (n_dpl .lt. 0) then
            write(1, 9100)h, tru_int, trap, int_trap_er,
     &           m2*h**2/12d0, int_trap_er/h**2
            goto 1000
         endif
         int_trap_er = abs(int_trap_er)
         n_dpl = -1
         call trap_int(a, b, n_ivl, f, n_dpl, tru_trap)
         trap_er = abs(tru_trap-trap)
         int_tru_trap_er = abs(tru_int-tru_trap)
         write(1, 9100)h, 
     &        tru_trap,
     &        trap, 
     &        trap_er, 
     &        1d0/(2d0*1d1**og_n_dpl), 
     &        int_tru_trap_er,
     &        m2*h**2/12d0,
     &        int_trap_er, 
     &        m2*h**2/12d0+1d0/(2d0*1d1**og_n_dpl)
      else
         call simp_int(a, b, n_ivl, f, n_dpl, simp)
         int_simp_er = tru_int - simp
         if (n_dpl .lt. 0) then
            write(1, 9100)h, 
     &           tru_int, 
     &           simp, 
     &           int_simp_er,
     &           m4*h**4/180d0, 
     &           int_simp_er/h**4
            goto 1000
         endif
         int_simp_er = abs(int_simp_er)
         n_dpl = -1
         call simp_int(a, b, n_ivl, f, n_dpl, tru_simp)
         simp_er = abs(tru_simp-simp)
         int_tru_simp_er = abs(tru_int-tru_simp)
         write(1, 9100)
     &        h, 
     &        tru_simp,
     &        simp, 
     &        simp_er, 
     &        1d0/(2d0*1d1**og_n_dpl), 
     &        int_tru_simp_er,
     &        m4*h**4/180d0,
     &        int_simp_er, 
     &        m4*h**4/180d0+1d0/(2d0*1d1**og_n_dpl)
      endif
 1000 continue
      call flush(1)
 9000 format('  h               I               T(h)    ',
     &     '        E(h)            M2*h^2/12       E(h)/h^2')
 9001 format('  h               I               S(h)    ',
     &     '        E(h)            M4*h^4/180      E(h)/h^4')
 9010 format(
     &     '  h             ',
     &     '  TruT(h)       ',
     &     '  T(h)          ',
     &     '  |TruT(h)-T(h)|',
     &     '  0.5*10^{-n}   ',
     &     '  |I-TruT(h)|   ',
     &     '  M2*h^2/12     ',
     &     '  |I-T(h)|      ',
     &     '  M2*h^2/12+0.5*10^{-n}')
 9011 format(
     &     '  h             ',
     &     '  TruS(h)       ',
     &     '  S(h)          ',
     &     '  |TruS(h)-S(h)|',
     &     '  0.5*10^{-n}   ',
     &     '  |I-TruS(h)|   ',
     &     '  M4*h^2/180    ',
     &     '  |I-S(h)|      ',
     &     '  M4*h^4/180+0.5*10^{-n}')
 9100 format(10(1x, g15.8))
      goto 10
      end

      double precision function f(x, n_dpl)
      implicit none
      double precision x
      integer n_dpl
      double precision v, f_v_t_n_dpl
      double precision v_t_n_dpl

      if (n_dpl .ge. 0) then
         v = exp(-x)
         f_v_t_n_dpl = v_t_n_dpl(v, n_dpl)
         f = f_v_t_n_dpl
      else
         f = exp(-x)
      endif
c      f = 1.0/(1.0+x*x)
c      f = sqrt(1.0-x*x)
      return
      end

