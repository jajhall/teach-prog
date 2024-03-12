      subroutine trap_int(a, b, n_ivl, f, n_dpl, trap)
      implicit none
      integer n_ivl, n_dpl
      double precision a, b, trap
      double precision f
      integer j
      double precision h, xj
      
      if (n_ivl .le. 0) goto 8000
      h = (b-a)/float(n_ivl)
      trap = f(a, n_dpl)
      do 10, j = 1, n_ivl-1
         xj = a + j*h
         trap = trap + 2.0*f(xj, n_dpl)
 10   continue
      trap = trap + f(b, n_dpl)
      trap = trap*(h/2.0)
 7000 continue
      return
 8000 continue
      print*,' trap_int called with n_ivl .le. 0'
      goto 7000
      end

      subroutine simp_int(a, b, n_ivl, f, n_dpl, simp)
      implicit none
      integer n_ivl, n_dpl
      double precision a, b, simp
      double precision f
      external f
      integer j
      double precision h, xj
      
      if (n_ivl .le. 0) goto 8000
      if (mod(n_ivl, 2) .ne. 0) goto 8010
      h = (b-a)/float(n_ivl)
      simp = f(a, n_dpl)
      do 10, j = 1, (n_ivl/2)-1
         xj = a + (2*j)*h
         simp = simp + 2.0*f(xj, n_dpl)
         xj = a + (2*j-1)*h
         simp = simp + 4.0*f(xj, n_dpl)
 10   continue
      xj = a + (n_ivl-1)*h
      simp = simp + 4.0*f(xj, n_dpl)
      simp = simp + f(b, n_dpl)
      simp = simp*(h/3.0)
 7000 continue
      return
 8000 continue
      print*,' simp_int called with n_ivl .le. 0'
      goto 7000
 8010 continue
      print*,' simp_int called with n_ivl odd'
      goto 7000
      end

      double precision function v_t_n_dpl(v, n)
      implicit none
      double precision v
      integer n
      double precision nw_v
      character*16 ch16
      character*7 ch7_fmt
      ch7_fmt = '(f16.0)'
      if (n .ge. 10) then
         print*,'n>9'
         stop
      endif
      if (n .ge. 0) ch7_fmt = '(f16.'//char(48+n)//')'
      write(ch16, ch7_fmt, err=8000)v
      read(ch16, *)nw_v
      v_t_n_dpl = nw_v
 7000 continue
      return
 8000 continue
      print*, 'Write error with ', v, n, ch16, ch7_fmt
      goto 7000
      end
