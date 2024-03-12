      program numdiff
      implicit none
      double precision x0, h, cent, fwd, tru_drv
      double precision tru_cent, drv_cent_er, cent_er, drv_tru_cent_er
      double precision tru_fwd, drv_fwd_er, fwd_er, drv_tru_fwd_er
      double precision m2, m3
      external f
      double precision f
      integer alg, og_n_dpl, n_dpl

      x0 = 1.0
      print*, ' Enter Fwd/Cent 1/0'
      read*, alg
c      alg = 1
      print*, ' Enter fd(x0)'
c      read*, tru_drv
      tru_drv = -2d0*exp(-1d0)
      print*, ' Enter M2, M3'
c      read*, M2, M3
      m2 = exp(-1d0)
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
      print*, ' Enter h '
      read*, h
      if (h .le. 0) stop

      m3 = h*exp(-1d0+h)
      n_dpl = og_n_dpl
      if (alg .eq. 1) then
         call fwd_diff_drv(x0, h, f, n_dpl, fwd)
         drv_fwd_er = tru_drv - fwd
         if (n_dpl .lt. 0) then
            write(1, 9100)h, tru_drv, fwd, drv_fwd_er,
     &           m2*h/2d0, drv_fwd_er/h
            goto 1000
         endif
         drv_fwd_er = abs(drv_fwd_er)
         n_dpl = -1
         call fwd_diff_drv(x0, h, f, n_dpl, tru_fwd)
         fwd_er = abs(tru_fwd-fwd)
         drv_tru_fwd_er = abs(tru_drv-tru_fwd)
         write(1, 9100)
     &        h, 
     &        tru_fwd,
     &        fwd, 
     &        fwd_er,
     &        1d0/(1d1**og_n_dpl*h),
     &        drv_tru_fwd_er,
     &        m2*h/2d0, 
     &        drv_fwd_er, 
     &        m2*h/2d0 + 1d0/(1d1**og_n_dpl*h)
      else
         call cent_diff_drv(x0, h, f, n_dpl, cent)
         drv_cent_er = tru_drv - cent
         if (n_dpl .lt. 0) then
            write(1, 9100)h, tru_drv, cent, drv_cent_er,
     &           m3*h*h/6d0, drv_cent_er/h**2
            goto 1000
         endif
         drv_cent_er = abs(drv_cent_er)
         n_dpl = -1
         call cent_diff_drv(x0, h, f, n_dpl, tru_cent)
         cent_er = abs(tru_cent-cent)
         drv_tru_cent_er = abs(tru_drv-tru_cent)
         write(1, 9100)
     &        h, 
     &        tru_cent,
     &        cent, 
     &        cent_er,
     &        1d0/(1d1**og_n_dpl*2*h),
     &        drv_tru_cent_er,
     &        m3*h*h/6d0, 
     &        drv_cent_er, 
     &        m3*h*h/6d0 + 1d0/(1d1**og_n_dpl*2*h)
      endif
 1000 continue
      call flush(1)
 9000 format('  h               I               T(h)    ',
     &     '        E(h)            M2*h/2          E(h)/h')
 9001 format('  h               I               S(h)    ',
     &     '        E(h)            M4*h^4/180      E(h)/h^4')
 9010 format(
     &     '  h             ',
     &     '  TruD(h)       ',
     &     '  D(h)          ',
     &     '  |TruD(h)-D(h)|',
     &     '  10^{-n}/h     ',
     &     '  |fd-TruD(h)|  ',
     &     '  M2*h/2        ',
     &     '  |fd-D(h)|     ',
     &     '  M2*h/2+10^-n/h')
 9011 format(
     &     '  h             ',
     &     '  TruD(h)       ',
     &     '  D(h)          ',
     &     '  |TruD(h)-D(h)|',
     &     '  10^{-n}/h     ',
     &     '  |fd-TruD(h)|  ',
     &     '  M3*h^2/6      ',
     &     '  |fd-D(h)|     ',
     &     '  M3*h^2/6+10^-n/h')
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
c         v = exp(-x)
         v = (2+x)*exp(-x)
         f_v_t_n_dpl = v_t_n_dpl(v, n_dpl)
         f = f_v_t_n_dpl
      else
c         f = exp(-x)
         f = (2+x)*exp(-x)
      endif
      return
      end

      subroutine fwd_diff_drv(x0, h, f, n_dpl, fwd)
      implicit none
      integer n_dpl
      double precision x0, h, f, fwd
      double precision f0, f1

      f0 = f(x0, n_dpl)
      f1 = f(x0+h, n_dpl)
      fwd = (f1-f0)/h
      return
      end
      
      subroutine cent_diff_drv(x0, h, f, n_dpl, cent)
      implicit none
      integer n_dpl
      double precision x0, h, f, cent
      double precision f0, f1

      f0 = f(x0-h, n_dpl)
      f1 = f(x0+h, n_dpl)
      cent = (f1-f0)/(2*h)
      return
      end
      
      
