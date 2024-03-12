      program lptex
      implicit none
      include 'LP.INC'
      integer r_n, c_n
      character*64 filename
      character*500 ch500
      character*1 backslash
      character*12 rl_t_TEX_ch10
      integer p

      backslash = char(92)
      lp_rd_cn = 1
      filename = 'test.lp'
      print*,' Enter Filename'
      read*, filename
      open(unit = lp_rd_cn, file = filename)
      call lptex_rd_lp
      close(lp_rd_cn)

      p = 0
      write(ch500(p+1:p+18+2*(n_c+1)), 9000)
     &     backslash, ('cr', c_n = 1, n_c+1)
      p = p + 18+2*(n_c+1)
      ch500(p+1:p+1) = '}'
      p = p + 1
      write(3, 9100)ch500(1:p)
c
      p = 0
      write(ch500(p+1:p+25), 9010)backslash
      p = p + 25
      call lptex_wr_lhs_cf(ch500, p, co)
      ch500(p+1:p+2) = backslash//backslash
      p = p + 2
      write(3, 9100)ch500(1:p)
c
      p = 0
      write(ch500(p+1:p+25), 9020)backslash
      p = p + 25
      call lptex_wr_lhs_cf(ch500, p, mtx(1, 1))
      ch500(p+1:p+5) = '&'//backslash//'le&'
      p = p + 5
      if (rhs(1) .eq. 1d0) then
         ch500(p+1:p+12) = '         1  '
      else
         ch500(p+1:p+12) = rl_t_TEX_ch10(rhs(1))
      endif
      p = p + 12
      ch500(p+1:p+2) = backslash//backslash
      p = p + 2
      write(3, 9100)ch500(1:p)

      do 10, r_n = 2, n_r
         p = 0
         write(ch500(p+1:p+25), 9030)
         p = p + 25
         call lptex_wr_lhs_cf(ch500, p, mtx(1, r_n))
         ch500(p+1:p+5) = '&'//backslash//'le&'
         p = p + 5
         if (rhs(r_n) .eq. 1d0) then
            ch500(p+1:p+12) = '         1  '
         else
            ch500(p+1:p+12) = rl_t_TEX_ch10(rhs(r_n))
         endif
         p = p + 12
         ch500(p+1:p+2) = backslash//backslash
         p = p + 2
         write(3, 9100)ch500(1:p)
 10   continue
      p = 0
      write(ch500(p+1:p+23), 9040)backslash, 2*n_c+1
      p = p + 23
      do 20, c_n = 1, n_c
         write(ch500(p+1:p+16), 9045)c_n, backslash, backslash
         p = p + 16
 20   continue
      p = p - 6
      ch500(p+1:p+2) = ',}'
      p = p + 2
      write(3, 9100)ch500(1:p)

      write(3, 9060)backslash
      stop
 9000 format('$$', a1, 'begin{array}{ll', 100a2)
 9010 format('       {', a1, 'rm maximize}&f&=')
 9020 format('     {', a1, 'rm subject~to}& & ')
 9030 format('                     & & ')
 9040 format('&&&', a1, 'multicolumn{', i2, '}{r}{')
 9045 format('x_{', i2, '}', a1, 'ge0', a1, 'quad ')
 9050 format('x_{', i2, '}', a1, 'ge0', a1, 'quad ')
 9060 format(a1, 'end{array}',/'$$')
 9100 format(a)
      end

      subroutine lptex_rd_lp
      implicit none
      include 'LP.INC'
      integer r_n, c_n
      read(lp_rd_cn, *)n_r, n_c, mx_mn
      if (n_r .gt. mx_n_r .or. n_c .gt. mx_n_c) then
         print*, 'Parameter error'
         stop
      endif
      read(lp_rd_cn, *)(co(c_n), c_n = 1, n_c), f0
      do r_n = 1, n_r
         read(lp_rd_cn, *)
     &        (mtx(c_n, r_n), c_n = 1, n_c), rhs(r_n)
      enddo
      return
      end

      subroutine lptex_wr_lhs_cf(ch500, p, lhs_cf)
      implicit none
      include 'LP.INC'
      integer c_n
      double precision lhs_cf(n_c)
      character*500 ch500
      integer p
      character*12 ch10_cf(mx_n_c)
      character*1 ch1_sgn(mx_n_c)
      character*12 rl_t_TEX_ch10
      double precision cf
      integer f_nz_c

      f_nz_c = -1
      do 10, c_n = 1, n_c
         cf = lhs_cf(c_n)
         if (cf .gt. 0d0) then
            ch1_sgn(c_n) = '+'
            ch10_cf(c_n) = rl_t_TEX_ch10(cf)
            if (f_nz_c .lt. 0) f_nz_c = c_n
         else if (cf .lt. 0d0) then
            ch1_sgn(c_n) = '-'
            cf = -cf
            ch10_cf(c_n) = rl_t_TEX_ch10(cf)
            if (f_nz_c .lt. 0) f_nz_c = c_n
         else
            ch1_sgn(c_n) = ' '
            ch10_cf(c_n) = '            '
         endif
 10   continue
      if (ch1_sgn(f_nz_c) .eq. '+') ch1_sgn(f_nz_c) = ' '
      cf = lhs_cf(1)
      if (cf .ne. 0d0) then
         write(ch500(p+1:p+18), 9001)ch1_sgn(1), ch10_cf(1), 1
         p = p + 18
      else
         ch500(p+1:p+18) = '&                '
         p = p + 18
      endif
      do 20, c_n = 2, n_c
         cf = lhs_cf(c_n)
         if (cf .ne. 0d0) then
            write(ch500(p+1:p+2), 9000)ch1_sgn(c_n)
            p = p + 2
            write(ch500(p+1:p+18), 9002)ch10_cf(c_n), c_n
            p = p + 18
         else
            ch500(p+1:p+2+18) = '& &                '
            p = p + 2
            p = p + 18
         endif
 20   continue
      return
 9000 format('&', a1)
 9001 format('&', a1, a10, 'x_{', i2, '}')
 9002 format('&', 1x, a10, 'x_{', i2, '}')
      end

      character*12 function rl_t_TEX_ch10(rl)
      implicit none
      logical rational_display
c      parameter (rational_display = .false.)
      parameter (rational_display = .true.)
      double precision eps
      parameter (eps=1d-8)
      double precision rl
      double precision rl_v
      character*12 ch10
      integer i_v, den

      if (rl .le. -1d32) then
         ch10 = ' -'//char(92)//'infty   '
      else if (rl .ge. 1d32) then
         ch10 = '  '//char(92)//'infty  '
      else if (rl .eq. 0d0) then
         ch10 = '          '
      else if (rl .eq. 1d0) then
         ch10 = '          '
      else
         if (rational_display) then
            rl_v = rl*(1d0+eps)
            i_v = rl_v
            rl_v = dble(i_v)
            if (abs(rl-rl_v) .lt. eps) then
               write(ch10, 9200)i_v
               goto 1000
            endif
            do 100, den = 2, 10
               rl_v = rl*(1d0+eps)*float(den)
               i_v = rl_v
               rl_v = dble(i_v)/dble(den)
               if (abs(rl-rl_v) .lt. eps) then
                  write(ch10, 9201)i_v, den
                  goto 1000
               endif
 100        continue
         endif
         if (abs(rl) .lt. 1d-1) then
            write(ch10, 9100)rl
         else if (abs(rl) .lt. 1d5) then
            write(ch10, 9000)rl
         else
            write(ch10, 9100)rl
         endif
      endif
 1000 continue
      rl_t_TEX_ch10 = ch10
      return
 9000 format(f10.4)
 9100 format(e10.4)
 9200 format(  i10)
 9201 format(i7, '/', i2)
      end
