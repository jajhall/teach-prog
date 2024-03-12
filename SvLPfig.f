      program lpfig
      implicit none
      include 'LPFIG.INC'
      double precision y_paper_ht
      integer r_n
      double precision scl_fac
      character*64 filename

      lp_rd_cn = 1
      filename = 'test.lp'
      print*,' Enter Filename'
      read*, filename
      open(unit = lp_rd_cn, file = filename)
      call lpfig_rd_lp
      close(lp_rd_cn)

      scl_fac = 1d0
      print*, 'Enter scale factor'
      read*, scl_fac
      call lpfig_iz_com(scl_fac)

      wr_cn = 1
      open(unit = wr_cn, file = 'test.fig')
      axis_mn_x = -2
      axis_mx_x = 6
      print*, 'Enter min and max for x1'
      read*, axis_mn_x, axis_mx_x
      axis_mn_x = min(axis_mn_x, 0d0)
      axis_mx_x = max(axis_mx_x, 0d0)

      axis_mn_y = -1
      axis_mx_y = 6
      print*, 'Enter min and max for x2'
      read*, axis_mn_y, axis_mx_y
      axis_mn_y = min(axis_mn_y, 0d0)
      axis_mx_y = max(axis_mx_y, 0d0)

      xscl = (mx_paper_wd-lh_axis_margin-rh_axis_margin)/
     &     float(axis_mx_x-axis_mn_x)
      yscl = -xscl
      y_paper_ht = -yscl*float(axis_mx_y-axis_mn_y)+
     &     top_axis_margin+bottom_axis_margin
      if (y_paper_ht .gt. mx_paper_ht) then
         yscl = -(mx_paper_ht-top_axis_margin-bottom_axis_margin)/
     &        float(axis_mx_y-axis_mn_y)
         xscl = -yscl
      endif
      x0 = lh_axis_margin -  xscl*axis_mn_x
      y0 = top_axis_margin - yscl*axis_mx_y
      if (mod(x0, cm2xfig) .gt. 0) x0 = x0 - mod(x0, cm2xfig) + cm2xfig
      if (mod(y0, cm2xfig) .gt. 0) y0 = y0 - mod(y0, cm2xfig) + cm2xfig

      lp_mn_x = axis_mn_x - lh_axis_margin/xscl
      lp_mx_x = axis_mx_x + rh_axis_margin/xscl
      lp_mn_y = axis_mn_y + bottom_axis_margin/yscl
      lp_mx_y = axis_mx_y - top_axis_margin/yscl
      lp_cs_wd = cs_wd/xscl
      fig_mn_x = x0 + xscl*axis_mn_x - lh_axis_margin
      fig_mx_x = x0 + xscl*axis_mx_x + rh_axis_margin
      fig_mn_y = y0 + yscl*axis_mx_y - top_axis_margin
      fig_mx_y = y0 + yscl*axis_mn_y + bottom_axis_margin

      call lpfig_wr_hdr
      call lpfig_wr_axes

      do 10, r_n = 1, n_r
         call lpfig_plot_cs(mtx(1, r_n), mtx(2, r_n), rhs(r_n))
 10   continue
      call lpfig_plot_cs(-1d0, 0d0, 0d0)
      call lpfig_plot_cs(0d0, -1d0, 0d0)
      stop
      end

      subroutine lpfig_rd_lp
      implicit none
      include 'LPFIG.INC'
      integer n_c, mx_mn
      integer r_n, c_n
      double precision f0
      read(lp_rd_cn, *)n_r, n_c, mx_mn
      if (n_r .gt. mx_n_r .or. n_c .ne. 2) then
         print*, 'Parameter error'
         stop
      endif
      read(lp_rd_cn, *)(co(c_n), c_n = 1, 2), f0
      read(lp_rd_cn, *)
     &     ((mtx(c_n, r_n), c_n = 1, 2), rhs(r_n), r_n = 1, n_r)
      return
      end

      subroutine lpfig_wr_hdr
      implicit none
      include 'LPFIG.INC'
      write(wr_cn, 9000)
      return
 9000 format('#FIG 3.2'/'Portrait'/'Center'/'Metric'/
     &     'A4      '/'100.00'/'Single'/'-2'/'1200 2'/)
      end

      subroutine lpfig_wr_axes
      implicit none
      include 'LPFIG.INC'
      character*1 backslash
      backslash = char(92)

      write(wr_cn, 9000)fig_mn_x, y0, fig_mx_x, y0
      write(wr_cn, 9010)fig_mx_x, y0+300, 'x', backslash
      write(wr_cn, 9012)fig_mx_x+50, y0+350, '1', backslash
      write(wr_cn, 9000)x0, fig_mx_y, x0, fig_mn_y
      write(wr_cn, 9010)x0-250, fig_mn_y, 'x', backslash
      write(wr_cn, 9012)x0-200, fig_mn_y+50, '2', backslash
      return
 9000 format('2 1 0 2 0 7 50 0 -1 0.000 0 0 -1 1 0 2'/
     &     '1 1 2.00 120.00 240.00'/4(2x, i6))
 9010 format('4 1 0 50 0 3 32 0.0000 0 225 210',
     &     2(2x, i8), 2(1x, a1), '001')
 9012 format('4 0 0 50 0 1 18 0.0000 0 225 210',
     &     2(2x, i8), 2(1x, a1), '001')
      end

      subroutine lpfig_plot_cs(a1, a2, b)
      implicit none
      include 'LPFIG.INC'
      double precision a1, a2, b
      double precision mu, c, d, dy
      double precision x1, y1, x2, y2, x3, y3, x4, y4

      if (a2 .eq. 0d0) then
         if (a1 .eq. 0d0) goto 8000
         c = b/a1
         if (c .lt. lp_mn_x .or. c .gt. lp_mx_x) goto 8000
         if (a1 .lt. 0d0) then
c
c     Plot x1 >= c
c
            call lpfig_wr_axis_v(.true., c)
            call lpfig_wr_cs(c, lp_mn_y, c, lp_mx_y,
     &           c-lp_cs_wd, lp_mx_y, c-lp_cs_wd, lp_mn_y)
         else
c
c     Plot x1 <= c
c
            call lpfig_wr_axis_v(.true., c)
            call lpfig_wr_cs(c, lp_mn_y, c, lp_mx_y,
     &           c+lp_cs_wd, lp_mx_y, c+lp_cs_wd, lp_mn_y)
         endif
      else if (a1 .eq. 0d0) then
         c = b/a2
         if (c .lt. lp_mn_y .or. c .gt. lp_mx_y) goto 8000
         if (a2 .lt. 0d0) then
c
c     Plot x2 >= c
c
            call lpfig_wr_axis_v(.false., c)
            call lpfig_wr_cs(lp_mn_x, c, lp_mx_x, c,
     &           lp_mx_x, c-lp_cs_wd, lp_mn_x, c-lp_cs_wd)
         else
c
c     Plot x2 <= c
c
            call lpfig_wr_axis_v(.false., c)
            call lpfig_wr_cs(lp_mn_x, c, lp_mx_x, c,
     &           lp_mx_x, c+lp_cs_wd, lp_mn_x, c+lp_cs_wd)
         endif
      else
         mu = -a1/a2
         c = b/a2
         if (c .le. lp_mx_y .and. c .ge. lp_mn_y) then
            call lpfig_wr_axis_v(.false., c)
         endif
         d = b/a1
         if (d .le. lp_mx_x .and. d .ge. lp_mn_x) then
            call lpfig_wr_axis_v(.true., d)
         endif
         x1 = lp_mx_x
         call lpfig_lin_f_p(mu, c, x1, y1)
         x2 = lp_mn_x
         call lpfig_lin_f_p(mu, c, x2, y2)
         dy = lp_cs_wd*sqrt(1d0+mu*mu)*sign(1d0, a2)
         x3 = lp_mn_x
         call lpfig_lin_f_p(mu, c+dy, x3, y3)
         x4 = lp_mx_x
         call lpfig_lin_f_p(mu, c+dy, x4, y4)
         call lpfig_wr_cs(x1, y1, x2, y2, x3, y3, x4, y4)
      endif
 7000 continue
      return
 8000 continue
      print*, 'Failed to plot constraint with data ', a1, a2, b
      goto 7000
 9000 format('2 1 0 2 0 7 50 0 -1 0.000 0 0 -1 1 0 2'/
     &     '1 1 2.00 120.00 240.00'/4(2x, i6))
      end

      subroutine lpfig_wr_cs(x1, y1, x2, y2, x3, y3, x4, y4)
      implicit none
      include 'LPFIG.INC'
      double precision x1, y1, x2, y2, x3, y3, x4, y4
      integer fig_x1, fig_y1, fig_x2, fig_y2
      integer fig_x3, fig_y3, fig_x4, fig_y4

      fig_x1 = x0 + xscl*x1
      fig_x2 = x0 + xscl*x2
      fig_x3 = x0 + xscl*x3
      fig_x4 = x0 + xscl*x4
      fig_y1 = y0 + yscl*y1
      fig_y2 = y0 + yscl*y2
      fig_y3 = y0 + yscl*y3
      fig_y4 = y0 + yscl*y4
      write(wr_cn, 9000)fig_x1, fig_y1, fig_x2, fig_y2
      write(wr_cn, 9010)fig_x1, fig_y1, fig_x2, fig_y2,
     &     fig_x3, fig_y3, fig_x4, fig_y4, fig_x1, fig_y1
      return
 9000 format('2 1 0 1 0 7 50 0 16 0.000 0 0 -1 0 0 2'/
     &     4(2x, i8))
 9010 format('2 1 0 0 0 7 100 0 16 0.000 0 0 -1 0 0 5'/
     &     10(2x, i8))
      end

      subroutine lpfig_wr_axis_v(x_axis, v)
      implicit none
      include 'LPFIG.INC'
      logical x_axis
      double precision v
      double precision frac_v
      integer fig_x_v, fig_y_v
      integer fig_dx_v, fig_dy_v
      integer f_ch
      character*1 backslash
      character*6 ch6_v

      if (v .eq. 0) goto 7000
      backslash = char(92)
      if (x_axis) then
         fig_x_v = x0 + xscl*v
         fig_y_v = y0
         fig_dy_v = fig_y_v + 0.1*cm2xfig
         write(wr_cn, 9000)fig_x_v, fig_y_v, fig_x_v, fig_dy_v
      else
         fig_x_v = x0
         fig_y_v = y0 + yscl*v
         fig_dx_v = fig_x_v-0.1*cm2xfig
         write(wr_cn, 9000)fig_x_v, fig_y_v, fig_dx_v, fig_y_v
      endif
      frac_v = v - aint(v)
      if (v .lt. -999999.0) then
         ch6_v = '******'
      else if (v .gt. 999999) then
         ch6_v = '******'
      else if (frac_v .ne. 0d0) then
         if (v .lt. -999.999) then
            write(ch6_v, '(i6)')int(v)
         else if (v .lt. -99.9999) then
            write(ch6_v, '(f6.1)')v
         else if (v .lt. -9.99999) then
            write(ch6_v, '(f6.2)')v
         else if (v .lt. -.999999) then
            write(ch6_v, '(f6.3)')v
         else if (v .gt. 999.999) then
            write(ch6_v, '(i6)')int(v)
         else if (v .gt. 99.9999) then
            write(ch6_v, '(f6.1)')v
         else if (v .gt. 9.99999) then
            write(ch6_v, '(f6.2)')v
         else if (v .gt. 9.99999) then
            write(ch6_v, '(f6.3)')v
         else
            write(ch6_v, '(f6.4)')v
         endif
      else
         write(ch6_v, '(i6)')int(v)
      endif
      do 10, f_ch = 1, 6
         if (ch6_v(f_ch:f_ch) .ne. ' ') goto 15
 10   continue
      f_ch = 6
 15   continue
      if (x_axis) then
         fig_dy_v = fig_y_v + 0.8*cm2xfig
         write(wr_cn, 9010)fig_x_v, fig_dy_v, ch6_v(f_ch:6), backslash
      else
         fig_dx_v = fig_x_v - cm2xfig
         fig_dy_v = fig_y_v + 0.2*cm2xfig
         write(wr_cn, 9010)fig_dx_v, fig_dy_v, ch6_v(f_ch:6), backslash
      endif
 7000 continue
      return
 9000 format('2 1 0 1 0 7 25 0 -1 0.000 0 0 -1 0 0 2'/
     &     4(2x, i8))
 9010 format('4 1 0 25 0 0 18 0.0000 4 0 0 ', 2(2x, i8), 1x, a, a1, 
     &     '001')
      end

      subroutine lpfig_lin_f_p(mu, c, x, y)
      implicit none
      include 'LPFIG.INC'
      double precision mu, c, x, y

      y = c + mu*x
      if (y .gt. lp_mx_y) then
         y = lp_mx_y
         x = (y-c)/mu
      else if (y .lt. lp_mn_y) then
         y = lp_mn_y
         x = (y-c)/mu
      endif
      return
      end

      subroutine lpfig_iz_com(scl_fac)
      implicit none
      include 'LPFIG.INC'
      double precision scl_fac

      cm2xfig = 450
      mx_paper_wd = 19d0*cm2xfig*scl_fac
      mx_paper_ht = 27d0*cm2xfig*scl_fac
      lh_axis_margin = 1d0*cm2xfig*scl_fac
      rh_axis_margin = 1d0*cm2xfig*scl_fac
      top_axis_margin = 1d0*cm2xfig*scl_fac
      bottom_axis_margin = 1d0*cm2xfig*scl_fac
      cs_wd = 0.5*cm2xfig*scl_fac
      return
      end

      
