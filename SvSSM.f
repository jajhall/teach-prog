c     Solves max mx_mn*c^Tx
c                  st  A^Tx <= b
c                         x >= 0
      implicit none
      integer mx_n_c, mx_n_r
      parameter (mx_n_c = 10, mx_n_r = 10)
      double precision tbu(0:mx_n_c+mx_n_r, 0:mx_n_r)
      double precision rao(0:mx_n_r)
      double precision aa, aa_1, aa_2, aa_mn, pv, mx_pv, mx_rdc_co, mu
      double precision xp_tau, xp_tl
      double precision pr_act(0:mx_n_c+mx_n_r)
      double precision pr_act_rec(0:mx_n_c+mx_n_r, 6)
      integer n_c, n_r, vr_in_r(0:mx_n_r)
      integer c_n, r_n, vr_n, vr_t_en_bs, vr_t_lv_bs, r_o_vr_t_lv_bs
      integer n_it, rec_n, rd_vr_t_en_bs
      integer mx_n_it
      double precision mx_mn
      character*64 filename

      print*,' Enter Filename'
      read*, filename
      open (unit = 1, file = filename)
      read(1,*)n_r, n_c, mx_mn
      if (n_r .gt. mx_n_r .or. n_c .gt. mx_n_c .or.
     &     abs(mx_mn) .ne. 1d0) stop
      read(1,*) (tbu(c_n,   0), c_n = 1, n_c), tbu(0, 0)
      read(1,*)((tbu(c_n, r_n), c_n = 1, n_c), tbu(0, r_n),
     &     r_n = 1, n_r)
      mx_n_it = 99999
c      read(1,*)mx_n_it
c      tbu(0, 0) = 0d0
      vr_in_r(0) = 0
      do 30, r_n = 1, n_r
         vr_in_r(r_n) = n_c+r_n
 30   continue
      do 40, c_n = 0, n_c
         tbu(c_n, 0) = mx_mn*tbu(c_n, 0)
         pr_act(c_n) = 0d0
         do 41, rec_n = 1, 6
            pr_act_rec(c_n, rec_n) = 0d0
 41      continue
 40   continue
      do 50, r_n = 1, n_r
         pr_act(vr_in_r(r_n)) = tbu(0, r_n)
 50   continue
      n_c = n_c + n_r
      n_it = 1
c      print*,' Enter EXPAND tau '
c      xp_tau = 4.9d-11
c      xp_tau = 1d-10
      xp_tau = 0d0
c      read*, xp_tau
      if (xp_tau .le. 0d0) then
         xp_tau = 0d0
         xp_tl = 0d0
      else
         read(1,*)xp_tl
c     xp_tl = 0.5d-6
c     xp_tl = 0d0
      endif
 100  continue
      rec_n = mod(n_it, 6) + 1
      do 110, c_n = 1, n_r
         vr_n = vr_in_r(c_n)
         do 120, r_n = 1, n_r
            if (r_n .eq. c_n) then
               tbu(vr_n, r_n) = 1d0
            else
               tbu(vr_n, r_n) = 0d0
            endif
 120     continue
 110  continue
c      call ssm_wr_tbu(2, 0, n_r, n_c, mx_n_c+mx_n_r, 
c     &     n_it, vr_in_r, tbu, pr_act, 0d0, rao, aa)
      do 130, c_n = 1, n_c
         if (pr_act(c_n) .lt. -(xp_tl+xp_tau)) then
            write(2, 9410)c_n, pr_act(c_n), -xp_tl
            write(6, 9410)c_n, pr_act(c_n), -xp_tl
         endif
         pr_act_rec(c_n, rec_n) = pr_act(c_n) - pr_act_rec(c_n, rec_n)
 130  continue
c      write(2, 9300)rec_n, (pr_act_rec(c_n, rec_n), c_n = 1, n_c)
      do 135, c_n = 1, n_c
         pr_act_rec(c_n, rec_n) = pr_act(c_n)
 135  continue
      do 140, r_n = 1, n_r
         if (tbu(0, r_n) .ne. pr_act(vr_in_r(r_n))) then
            print*,' ERROR: tbu(0, r_n) .ne. pr_act(vr_in_r(r_n)) '
            print*,r_n, tbu(0, r_n), vr_in_r(r_n), pr_act(vr_in_r(r_n))
            stop
         endif
 140  continue
      vr_t_en_bs = 0
      mx_rdc_co = 0d0
      do 150, c_n = 1, n_c
         if (tbu(c_n, 0) .gt. mx_rdc_co) then
            vr_t_en_bs = c_n
            mx_rdc_co = tbu(c_n, 0)
         endif
 150  continue
      write(2, 9000)rec_n, vr_t_en_bs, mx_rdc_co
 200  continue
      rd_vr_t_en_bs = vr_t_en_bs
c      read*, rd_vr_t_en_bs
      if (rd_vr_t_en_bs .ge. 0) vr_t_en_bs = rd_vr_t_en_bs
      if (vr_t_en_bs .eq. 0) goto 1000
      if (vr_t_en_bs .lt. 0 .or. vr_t_en_bs .gt. n_c) then
         print*, ' Invalid column choice '
         goto 200
      endif
      xp_tl = xp_tl + xp_tau
      r_o_vr_t_lv_bs = 0
      aa_1 = 1d32
      do 210, r_n = 1, n_r
         if (tbu(vr_t_en_bs, r_n) .le. 0d0) then
            rao(r_n) =  1d32
         else
            rao(r_n) = (tbu(0, r_n)+xp_tl)/tbu(vr_t_en_bs, r_n)
            if (rao(r_n) .lt. aa_1) then
               r_o_vr_t_lv_bs = r_n
               aa_1 = rao(r_n)
            endif
         endif
 210  continue
      call ssm_wr_tbu(2, 1, n_r, n_c, mx_n_c+mx_n_r,
     &     n_it, vr_in_r, tbu, pr_act, xp_tl, rao, aa_1)
      call ssm_wr_tbu(3, 1, n_r, n_c, mx_n_c+mx_n_r,
     &     n_it, vr_in_r, tbu, pr_act, xp_tl, rao, aa_1)
c      write(3, 9200)vr_t_en_bs, tbu(vr_t_en_bs, 0)
      if (xp_tau .gt. 0) then
         write(2, 9100)rec_n, r_o_vr_t_lv_bs, aa_1
         r_o_vr_t_lv_bs = 0
         aa_2 = 1d32
         mx_pv = 0d0
         do 220, r_n = 1, n_r
            if (tbu(vr_t_en_bs, r_n) .le. 0d0) then
               rao(r_n) =  1d32
            else
               rao(r_n) = tbu(0, r_n)/tbu(vr_t_en_bs, r_n)
               if (rao(r_n) .le. aa_1) then
                  pv = abs(tbu(vr_t_en_bs, r_n))
                  if (pv .gt. mx_pv) then
                     r_o_vr_t_lv_bs = r_n
                     aa_2 = rao(r_n)
                     mx_pv = pv
                  endif
               endif
            endif
 220     continue
         if (r_o_vr_t_lv_bs .ne. mod(n_it+1, 2)+1) then
            write(2,*)n_it, vr_t_en_bs, r_o_vr_t_lv_bs, xp_tl
            write(2, 9990)
     &           (tbu(0, 1)+xp_tl)/tbu(vr_t_en_bs, 1),
     &           (tbu(0, 2)+xp_tl)/tbu(vr_t_en_bs, 2), aa_1
            write(2, 9990)rao(1), rao(2), aa_2, aa_mn, aa
            write(2, 9990)((tbu(c_n, r_n), c_n=1,n_c),
     &           tbu(0, r_n), r_n=1,n_r)
            write(2, 9990)(tbu(c_n, 0), c_n=1,n_c)
 9990       format(7(2x,g14.8))
         endif
         if (mod(n_it, 2) .eq. 0)
     &        print*, (tbu(0, 1)+xp_tl)/tbu(vr_t_en_bs, 1),
     &        tbu(0, 2)/tbu(vr_t_en_bs, 2),
     &        (tbu(0, 1)+xp_tl)/tbu(vr_t_en_bs, 1)-
     &           tbu(0, 2)/tbu(vr_t_en_bs, 2)
         call ssm_wr_tbu(2, 1, n_r, n_c, mx_n_c+mx_n_r,
     &        n_it, vr_in_r, tbu, pr_act, 0d0, rao, aa_2)
         call ssm_wr_tbu(3, 1, n_r, n_c, mx_n_c+mx_n_r,
     &        n_it, vr_in_r, tbu, pr_act, xp_tl, rao, aa_1)
         aa_mn = xp_tau/mx_pv
         if (aa_2 .ge. aa_mn) then
            aa = aa_2
         else
            aa = aa_mn
         endif
         write(2, 9110)rec_n, r_o_vr_t_lv_bs, aa_2, mx_pv, aa
c         write(3, 9210)vr_in_r(r_o_vr_t_lv_bs), vr_t_en_bs, aa
      else
         aa = aa_1
         write(2, 9120)rec_n, r_o_vr_t_lv_bs, aa
c         write(3, 9210)vr_in_r(r_o_vr_t_lv_bs), vr_t_en_bs, aa
      endif
 300  continue
      r_n = r_o_vr_t_lv_bs
c      read*, r_n
      if (r_n .eq. 0) goto 1000
      if (r_n .lt. 0 .or. r_n .gt. n_r) then
         print*, ' Invalid row choice '
         goto 300
      endif
      pv = tbu(vr_t_en_bs, r_o_vr_t_lv_bs)
      if (r_n .ne. r_o_vr_t_lv_bs) then
         r_o_vr_t_lv_bs = r_n
         aa = rao(r_o_vr_t_lv_bs)
      endif
      vr_t_lv_bs = vr_in_r(r_o_vr_t_lv_bs)
      do 310, r_n = 1, n_r
         vr_n = vr_in_r(r_n)
         pr_act(vr_n) = pr_act(vr_n) - aa*tbu(vr_t_en_bs, r_n)
         tbu(0, r_n) = pr_act(vr_n)
 310  continue
      tbu(0, 0) = tbu(0, 0) + aa*tbu(vr_t_en_bs, 0)
      pr_act(vr_t_en_bs) = pr_act(vr_t_en_bs) + aa
      tbu(0, r_o_vr_t_lv_bs) = pr_act(vr_t_en_bs)
      do 320, c_n = 1, n_c
         tbu(c_n, r_o_vr_t_lv_bs) = tbu(c_n, r_o_vr_t_lv_bs)/pv
 320  continue
      do 340, r_n = 0, n_r
         if (r_n .eq. r_o_vr_t_lv_bs) goto 340
         mu = tbu(vr_t_en_bs, r_n)
         do 330, c_n = 1, n_c
            tbu(c_n, r_n) = tbu(c_n, r_n) - tbu(c_n, r_o_vr_t_lv_bs)*mu
 330     continue
 340  continue
      vr_in_r(r_o_vr_t_lv_bs) = vr_t_en_bs
      n_it = n_it + 1
      if (n_it .le. mx_n_it) goto 100
 1000 continue
      call ssm_wr_tbu(2, 0, n_r, n_c, mx_n_c+mx_n_r, 
     &     n_it, vr_in_r, tbu, pr_act, 0d0, rao, aa)
      call ssm_wr_tbu(3, 0, n_r, n_c, mx_n_c+mx_n_r,
     &     n_it, vr_in_r, tbu, pr_act, 0d0, rao, aa)
      stop
 9000 format(' CHUZC ', i1, ': vr_t_en_bs = ', i4,
     &     ' has reduced cost = ', g10.4,
     &     ' Enter vr_t_en_bs')
 9100 format(' EXPAND ', i1, ' pass 1: r_o_vr_t_lv_bs = ', i4,
     &     ' gives aa_1 = ', g10.4)
 9110 format(' EXPAND ', i1, ' pass 2: r_o_vr_t_lv_bs = ', i4,
     &     ' gives aa_2 = ', g10.4, ' for mx_pv = ', g10.4,
     &     ' so aa = ', g10.4)
 9120 format(' CHUZR ', i1, ': r_o_vr_t_lv_bs = ', i4,
     &     ' gives aa = ', g10.4,
     &     ' Enter r_o_vr_t_lv_bs ')
 9200 format(/'Consider increasing nonbasic variable $x_', i1,
     &     '$ which has maximum objective coefficient of ', g11.4)
 9210 format(/'The first basic variable to be zeroed is $x_', i1,
     &     '$ when nonbasic variable $x_', i1, '$ reaches ', g11.4)
 9300 format(' DELTA ', i1, ':  ', 6(3x, g10.4))
 9410 format(' ERROR: Variable ', i4, ' has activity ', g10.4,
     &     ' which violates the tolerance ', g10.4)
      end

      subroutine ssm_wr_tbu(cn_n, n_aux, n_r, n_c, mx_n_c,
     &     n_it, vr_in_r, tbu, pr_act, tl, rao, aa)
      implicit none
      integer cn_n, n_aux, n_r, n_c, mx_n_c, n_it, vr_in_r(0:n_r)
      double precision tbu(0:mx_n_c, 0:n_r), pr_act(0:n_c), rao(0:n_r)
      double precision tl, aa
      character*12 rl_t_TEX_ch12
      integer r_n, ln, c_n, p
      character*250 tbu_r_li
      character*12 TEX_ch(1:10)
      character*1 backslash
      logical w_rao

      w_rao = n_aux .gt. 0

      backslash = char(92)
      if (cn_n .eq. 3) then
         write(cn_n, '(a1, a7/)')backslash, 'medskip'

         p = 29+n_c
         write(tbu_r_li(1:p), 9200)backslash, backslash, 
     &        ('r', c_n = 1, n_c)
         if (w_rao) then
            write(tbu_r_li(p+1:p+17), 9201)'r', backslash, n_c+1
         else
            write(tbu_r_li(p+1:p+17), 9201)' ', backslash, n_c+1
         endif
         p = p + 17
         write(cn_n, '(a)')tbu_r_li(1:p)

         p = 14*n_c
         write(tbu_r_li(1:p), 9210)(c_n, c_n = 1, n_c)
         if (w_rao) then
            write(tbu_r_li(p+1:p+38), 9211)'&    Ratio  ', 
     &           backslash, backslash, backslash, n_c+1
         else
            write(tbu_r_li(p+1:p+38), 9211)'            ', 
     &           backslash, backslash, backslash, n_c+1
         endif
         p = p + 38
         write(cn_n, '(a)')tbu_r_li(1:p)

         do 20, r_n = 1, n_r
            do 10, c_n = 1, n_c
               TEX_ch(c_n) = rl_t_TEX_ch12(tbu(c_n, r_n))
 10         continue
            TEX_ch(n_c+1) = rl_t_TEX_ch12(tbu(0, r_n))
            if (w_rao) then
               TEX_ch(n_c+2) = rl_t_TEX_ch12(rao(r_n))
               p = 14*(n_c+2)
               write(tbu_r_li(1:p), 9220)(TEX_ch(c_n), c_n = 1, n_c+2)
            else
               p = 14*(n_c+1)
               write(tbu_r_li(1:p), 9220)(TEX_ch(c_n), c_n = 1, n_c+1)
            endif
            write(tbu_r_li(p+1:p+2), 9221)backslash, backslash
            p = p + 2
            if (r_n .eq. n_r) then
               write(tbu_r_li(p+1:p+11), 9222)backslash, n_c+1
               p = p + 11
            endif
            write(cn_n, '(a)')tbu_r_li(1:p)
 20      continue
         do 30, c_n = 1, n_c
            TEX_ch(c_n) = rl_t_TEX_ch12(tbu(c_n, 0))
 30      continue
         TEX_ch(n_c+1) = rl_t_TEX_ch12(-tbu(0, 0))
         p = 14*(n_c+1)
         write(tbu_r_li(1:p), 9220)(TEX_ch(c_n), c_n = 1, n_c+1)
         write(tbu_r_li(p+1:p+13), 9223)
     &        backslash, backslash, backslash, n_c+1
         p = p + 13
         write(cn_n, '(a)')tbu_r_li(1:p)
         write(cn_n, 9240)backslash, backslash
      else
         write(cn_n, *)
         call ssm_g_tbu_hd_li(n_aux, n_c, n_it, tbu_r_li, ln)
         write(cn_n, 9010)tbu_r_li(1:ln)
         do 110, r_n = 1, n_r
            call ssm_g_tbu_r_li(n_aux, 
     &           vr_in_r(r_n), n_c, tbu(0, r_n), tl, rao(r_n),
     &           tbu_r_li, ln)
            write(cn_n, 9010)tbu_r_li(1:ln)
 110     continue
         call ssm_g_tbu_r_li(n_aux, 0, n_c, tbu(0, 0), tl, aa,
     &        tbu_r_li, ln)
         write(cn_n, 9010)tbu_r_li(1:ln)
         call ssm_g_tbu_r_li(0, -n_it, n_c, pr_act, tl, aa,
     &        tbu_r_li, ln)
         write(cn_n, 9010)tbu_r_li(1:ln)
      endif
      return
 9010 format(a)
 9200 format(a1, 'centerline{', a1, 'begin{tabular}{|',10a1)
 9201 format('|r|', a1, 'l}', a1, 'cline{1-', i1, '}')

 9210 format(10(5x, i3, 5x, '&'))
 9211 format('    RHS      ', a12, 3a1, 'cline{1-', i1, '}')

 9220 format(12(1x, a12, '&'))
 9221 format(2a1)
 9222 format( a1, 'cline{1-', i1, '}')
 9223 format(3a1, 'cline{1-', i1, '}')
 9230 format(12(1x, a12, '&'), '\\')
 9240 format(a1, 'end{tabular}}', a1, 'medskip')
      end

      subroutine ssm_g_tbu_hd_li(n_aux, n_c, n_it, tbu_r_li, ln)
      implicit none
      integer n_aux, n_c, n_it, ln
      character*250 tbu_r_li
      character*4 i_t_ch4
      integer c_n

      ln = 1
      write(tbu_r_li(ln:ln+12), 9000)n_it
      ln = ln + 13
      do 10, c_n = 1, n_c-1
         tbu_r_li(ln:ln+2) = '   '
         tbu_r_li(ln+3:ln+6) = i_t_ch4(c_n)
         tbu_r_li(ln+7:ln+12) = '      '
         ln = ln + 13
 10   continue
      tbu_r_li(ln:ln+2) = '   '
      tbu_r_li(ln+3:ln+6) = i_t_ch4(c_n)
      tbu_r_li(ln+7:ln+12) = '     |'
      ln = ln + 13
      tbu_r_li(ln:ln+12) = '     RHS    |'
      ln = ln + 13
      if (n_aux .eq. 1) then
         tbu_r_li(ln:ln+12) = '    Ratio   |'
         ln = ln + 13
      endif
      ln = ln - 1
      return
 9000 format(' | It', i6, ' |')
      end
      
      subroutine ssm_g_tbu_r_li(n_aux, r_n, n_c, tbu_r, tl, rao,
     &     tbu_r_li, ln)
      implicit none
      integer n_aux, r_n, n_c, ln
      double precision tbu_r(0:n_c), tl, rao
      character*250 tbu_r_li
      character*4 i_t_ch4
      character*10 rl_t_ch10
      integer c_n
      double precision v

      ln = 1
      if (r_n .lt. 0) then
         tbu_r_li(ln:ln+12) = ' | Pr_Act '//char(49+mod(-r_n, 6))//' |'
      else if (r_n .eq. 0) then
         tbu_r_li(ln:ln+12) = ' | Rdc_Cost |'
         v = tbu_r(0)
      else
         tbu_r_li(ln:ln+2) = ' | '
         tbu_r_li(ln+3:ln+6) = i_t_ch4(r_n)
         tbu_r_li(ln+7:ln+12) = '     |'
         v = tbu_r(0) + tl
      endif
      ln = ln + 13
      do 10, c_n = 1, n_c-1
         tbu_r_li(ln:ln) = ' '
         tbu_r_li(ln+1:ln+10) = rl_t_ch10(tbu_r(c_n))
         tbu_r_li(ln+11:ln+12) = '  '
         ln = ln + 13
 10   continue
      tbu_r_li(ln:ln) = ' '
      tbu_r_li(ln+1:ln+10) = rl_t_ch10(tbu_r(c_n))
      tbu_r_li(ln+11:ln+12) = ' |'
      ln = ln + 13
      if (r_n .ge. 0) then
         tbu_r_li(ln:ln) = ' '
         tbu_r_li(ln+1:ln+10) = rl_t_ch10(v)
         tbu_r_li(ln+11:ln+12) = ' |'
         ln = ln + 13
      endif
      if (n_aux .eq. 1) then
         tbu_r_li(ln:ln) = ' '
         tbu_r_li(ln+1:ln+10) = rl_t_ch10(rao)
         tbu_r_li(ln+11:ln+12) = ' |'
         ln = ln + 13
      endif
      ln = ln - 1
      return
      end
      
      character*4 function i_t_ch4(i)
      implicit none
      integer i
      character*4 ch4
      write(ch4, 9000)i
      i_t_ch4 = ch4
      return
 9000 format(i4)
      end

      character*10 function rl_t_ch10(rl)
      implicit none
      double precision rl
      character*10 ch10
      if (rl .le. -1d32) then
         ch10 = ' -infinity'
      else if (rl .ge. 1d32) then
         ch10 = ' +infinity'
      else if (rl .eq. 0d0) then
         ch10 = '     .    '
      else if (rl .eq. 1d0) then
         ch10 = '    1.0   '
      else if (abs(rl) .lt. 1d-1) then
         write(ch10, 9100)rl
      else if (abs(rl) .lt. 1d5) then
         write(ch10, 9000)rl
      else
         write(ch10, 9100)rl
      endif
      rl_t_ch10 = ch10
      return
 9000 format(f10.4)
 9100 format(e10.4)
      end

      character*12 function rl_t_TEX_ch12(rl)
      implicit none
      double precision rl
      double precision rl_v
      character*12 ch12
      integer i_v, den

      if (rl .le. -1d32) then
         ch12 = '$ -'//char(92)//'infty  $ '
      else if (rl .ge. 1d32) then
         ch12 = '$  '//char(92)//'infty  $'
      else if (rl .eq. 0d0) then
         ch12 = '$    0     $'
      else if (rl .eq. 1d0) then
         ch12 = '$    1     $'
      else
         i_v = rl
         if (rl .eq. dble(i_v)) then
            write(ch12, 9200)i_v
            goto 1000
         endif
         do 100, den = 2, 10
            rl_v = rl*float(den)
            i_v = rl_v
            rl_v = dble(i_v)/dble(den)
            if (abs(rl-rl_v) .lt. 1d-8) then
               write(ch12, 9201)i_v, den
               goto 1000
            endif
 100     continue
         if (abs(rl) .lt. 1d-1) then
            write(ch12, 9100)rl
         else if (abs(rl) .lt. 1d5) then
            write(ch12, 9000)rl
         else
            write(ch12, 9100)rl
         endif
      endif
 1000 continue
      rl_t_TEX_ch12 = ch12
      return
 9000 format('$', f10.4, '$')
 9100 format('$', e10.4, '$')
 9200 format('$',   i10, '$')
 9201 format('$',i7, '/', i2, '$')
      end
