c     Solves max mx_mn*c^Tx
c                  st  A^Tx <= b
c                         x >= 0
      implicit none
      include 'SSM.INC'
      double precision inf
      parameter (inf = 1d31)
      double precision tbu(0:mx_n_tbu_c, 0:mx_n_r)
      double precision rao(0:mx_n_r)
      double precision aa, aa_1, aa_2, aa_mn, pv, mx_pv
      double precision mx_ph1_rdc_co, mx_rdc_co, mu
      double precision xp_tau, xp_tl
      double precision ph_1_rdc_co(0:mx_n_tbu_c)
      double precision pr_act(0:mx_n_tbu_c)
      double precision pr_act_rec(0:mx_n_tbu_c, 6)
      double precision vr_ub(0:mx_n_tbu_c)
      logical at_lb(0:mx_n_tbu_c)
      integer mv_dir
      integer ml_n_c, n_c, n_r, vr_in_r(0:mx_n_r)
      integer c_n, r_n, lc_r_n, vr_n
      integer vr_t_en_bs, vr_t_lv_bs, r_o_vr_t_lv_bs
      integer n_it, rec_n, rd_vr_t_en_bs
      integer mx_n_it
      double precision mx_mn
      logical ph_1
      logical er_fd
      logical ck_vr_ub
      logical rp_bs_cg
      parameter (rp_bs_cg = .false.)
      character*64 filename
      character*80 ch80_li

      cn3_tbu_rational = .false.

      filename = 'test.lp'
      print*,' Enter Filename'
      read*, filename
      open (unit = 1, file = filename)
      read(1,*)n_r, n_c, mx_mn
      if (n_r .gt. mx_n_r .or. n_c+2*n_r .gt. mx_n_tbu_c .or.
     &     abs(mx_mn) .ne. 1d0) stop
      ml_n_c = n_c
      read(1,*) (tbu(c_n,   0), c_n = 1, n_c), tbu(0, 0)
      do r_n = 1, n_r
         read(1,*)(tbu(c_n, r_n), c_n = 1, n_c), tbu(0, r_n)
      enddo
      do 5, c_n = 1, mx_n_tbu_c
         vr_ub(c_n) = inf
         at_lb(c_n) = .true.
 5    continue
      ck_vr_ub = .false.
      read(1, *, end=10)(vr_ub(c_n), c_n = 1, n_c)
      ck_vr_ub = .true.
 10   continue
      mx_n_it = 99999
      print*,'Read model'
c      read(1,*)mx_n_it
c      tbu(0, 0) = 0d0
      do 40, c_n = 0, n_c
         tbu(c_n, 0) = mx_mn*tbu(c_n, 0)
         pr_act(c_n) = 0d0
         do 41, rec_n = 1, 6
            pr_act_rec(c_n, rec_n) = 0d0
 41      continue
         ph_1_rdc_co(c_n) = 0d0
 40   continue
      ph_1 = .false.
      vr_in_r(0) = 0
      ph_1_rdc_co(0) = 0d0
      do 45, r_n = 1, n_r
         if (tbu(0, r_n) .lt. 0d0) then
            ph_1 = .true.
            n_c = n_c + 1
            vr_in_r(r_n) = n_c + n_r
c
c     Add the slack column explicitly as a negated identity column
c
            do 43, lc_r_n = 1, n_r
               tbu(ml_n_c+r_n, lc_r_n) = 0d0
 43         continue
            tbu(ml_n_c+r_n, r_n) = -1d0
            ph_1_rdc_co(ml_n_c+r_n) = -1d0
            do 44, c_n = 1, ml_n_c
               tbu(c_n, r_n) = -tbu(c_n, r_n)
               ph_1_rdc_co(c_n) = ph_1_rdc_co(c_n) + tbu(c_n, r_n)
 44         continue
            tbu(0, r_n) = -tbu(0, r_n)
            ph_1_rdc_co(0) = ph_1_rdc_co(0) - tbu(0, r_n)
            ph_1_rdc_co(n_c+n_r) = 0d0
            pr_act(n_c) = 0d0
         else
            vr_in_r(r_n) = ml_n_c + r_n
            ph_1_rdc_co(ml_n_c+r_n) = 0d0
         endif
 45   continue
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
c     Zero the rest of the tableau
      do c_n = ml_n_c+1, n_c
         do r_n = 0, n_r
            tbu(c_n, r_n) = 0
         enddo
      enddo

 100  continue
       
      rec_n = mod(n_it, 6) + 1
      if (rp_bs_cg) then
         if (n_it .gt. 0 .and. mod(n_it, 2) .eq. 1) then
            write(3, '(a1, a)')char(92), 'newpage'
         endif
      endif
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
c      call ssm_wr_tbu(2, 0, n_r, n_c, 
c     &     n_it, vr_in_r, tbu, pr_act, ph_1_rdc_co, 0d0, rao, aa)
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
      if (ph_1) then
         mx_ph1_rdc_co = 0d0
         mx_rdc_co = -inf
         do 151, c_n = 1, n_c
            if (at_lb(c_n)) then
               if (ph_1_rdc_co(c_n) .gt. mx_ph1_rdc_co) then
                  vr_t_en_bs = c_n
                  mx_ph1_rdc_co = ph_1_rdc_co(c_n)
                  mx_rdc_co = tbu(c_n, 0)
               else if (ph_1_rdc_co(c_n) .ge. mx_ph1_rdc_co .and.
     &                 tbu(c_n, 0) .gt. mx_rdc_co) then
                  vr_t_en_bs = c_n
                  mx_rdc_co = tbu(c_n, 0)
               endif
            else
               if (-ph_1_rdc_co(c_n) .gt. mx_ph1_rdc_co) then
                  vr_t_en_bs = c_n
                  mx_ph1_rdc_co = -ph_1_rdc_co(c_n)
                  mx_rdc_co = -tbu(c_n, 0)
               else if (-ph_1_rdc_co(c_n) .ge. mx_ph1_rdc_co .and.
     &                 -tbu(c_n, 0) .gt. mx_rdc_co) then
                  vr_t_en_bs = c_n
                  mx_rdc_co = -tbu(c_n, 0)
               endif
            endif
 151     continue
      else
         mx_rdc_co = 0d0
         do 152, c_n = 1, n_c
            if (at_lb(c_n)) then
               if (tbu(c_n, 0) .gt. mx_rdc_co) then
                  vr_t_en_bs = c_n
                  mx_rdc_co = tbu(c_n, 0)
               endif
            else
               if (-tbu(c_n, 0) .gt. mx_rdc_co) then
                  vr_t_en_bs = c_n
                  mx_rdc_co = -tbu(c_n, 0)
               endif
            endif
 152     continue
      endif
      mx_rdc_co = 1d3*ph_1_rdc_co(vr_t_en_bs) + tbu(vr_t_en_bs, 0)
      write(2, 9000)rec_n, vr_t_en_bs, mx_rdc_co
 200  continue
      rd_vr_t_en_bs = vr_t_en_bs
      write(*, '(a, i2)')
     &     'Enter variable to enter basis - SSM chooses ', vr_t_en_bs
c      read(5, *, end=201, err=201)rd_vr_t_en_bs
      if (rd_vr_t_en_bs .ge. 0) vr_t_en_bs = rd_vr_t_en_bs
 201  continue
      if (vr_t_en_bs .eq. 0) goto 1000
      if (vr_t_en_bs .lt. 0 .or. vr_t_en_bs .gt. n_c) then
         print*, ' Invalid column choice '
         goto 200
      endif
      write(ch80_li, '(a, i2)')'Entering variable is ', vr_t_en_bs
      print*, ch80_li
      if (rp_bs_cg) write(3, '(/a)')ch80_li
      mx_rdc_co = 1d3*ph_1_rdc_co(vr_t_en_bs) + tbu(vr_t_en_bs, 0)
      if (mx_rdc_co .gt. 0d0) then
         mv_dir =  1
         if (.not. at_lb(vr_t_en_bs)) then
            print*, 'STRANGE: mv_dir = 1',
     &           ' but not at lower bound for variable ', vr_t_en_bs
            stop
         endif
      else
         mv_dir = -1
         if (at_lb(vr_t_en_bs)) then
            print*, 'STRANGE: mv_dir = -1',
     &           ' but at lower bound for variable ', vr_t_en_bs
            stop
         endif
      endif
      xp_tl = xp_tl + xp_tau
      vr_in_r(0) = vr_t_en_bs
      if (at_lb(vr_t_en_bs)) then
         if (vr_ub(vr_t_en_bs) .lt. inf) then
            rao(0) = vr_ub(vr_t_en_bs)!-vr_lb(vr_t_en_bs)
         else       
            rao(0) = 1d32
         endif
      else
!         if (vr_lb(vr_t_en_bs) .gt. -inf) then
            rao(0) = vr_ub(vr_t_en_bs)!-vr_lb(vr_t_en_bs)
!         else       
!            rao(0) = 1d32
!         endif
      endif
      r_o_vr_t_lv_bs = 0
      aa_1 = rao(0)
      do 210, r_n = 1, n_r
         pv = mv_dir*tbu(vr_t_en_bs, r_n)
         rao(r_n) =  1d32
         if (abs(pv) .le. 1d-8) goto 210
         vr_n = vr_in_r(r_n)
         if (pv .le. 0d0) then
            if (vr_ub(vr_n) .lt. inf) 
     &           rao(r_n) = -(vr_ub(vr_n)+xp_tl-tbu(0, r_n))/pv
         else
            rao(r_n) = (tbu(0, r_n)+xp_tl)/pv
         endif
         if (rao(r_n) .lt. aa_1) then
            r_o_vr_t_lv_bs = r_n
            aa_1 = rao(r_n)
         endif
 210  continue
      call ssm_wr_tbu(2, 1, n_r, n_c, 
     &     n_it, vr_in_r, tbu, pr_act, ph_1_rdc_co, xp_tl, rao, aa_1)
      call ssm_wr_tbu(3, 1, n_r, n_c, 
     &     n_it, vr_in_r, tbu, pr_act, ph_1_rdc_co, xp_tl, rao, aa_1)
      call ssm_wr_tbu(6, 1, n_r, n_c, 
     &     n_it, vr_in_r, tbu, pr_act, ph_1_rdc_co, xp_tl, rao, aa_1)
      call flush(2)
c      write(3, 9200)vr_t_en_bs, tbu(vr_t_en_bs, 0)
c      if (xp_tau .gt. 0d0) then
      if (xp_tl .gt. 0d0) then
         write(2, 9100)rec_n, r_o_vr_t_lv_bs, aa_1
         if (at_lb(vr_t_en_bs)) then
            if (vr_ub(vr_t_en_bs) .lt. inf) then
               rao(0) = vr_ub(vr_t_en_bs)+xp_tl !-vr_lb(vr_t_en_bs)
            else       
               rao(0) = 1d32
            endif
         else
!         if (vr_lb(vr_t_en_bs) .gt. -inf) then
            rao(0) = vr_ub(vr_t_en_bs)+xp_tl !-vr_lb(vr_t_en_bs)
!         else       
!            rao(0) = 1d32
!         endif
         endif
         r_o_vr_t_lv_bs = 0
         mx_pv = 1d0
         aa_2 = rao(0)


         do 220, r_n = 1, n_r
            pv = mv_dir*tbu(vr_t_en_bs, r_n)
            vr_n = vr_in_r(r_n)
            rao(r_n) =  1d32
            if (abs(pv) .le. 1d-8) goto 220
            if (pv .le. 0d0) then
               if (vr_ub(vr_n) .lt. inf) 
     &              rao(r_n) = -(vr_ub(vr_n)+xp_tl-tbu(0, r_n))/pv
            else
               rao(r_n) = (tbu(0, r_n)+xp_tl)/pv
            endif
            if (rao(r_n) .le. aa_1) then
               pv = abs(pv)
               if (pv .gt. mx_pv) then
                  r_o_vr_t_lv_bs = r_n
                  aa_2 = rao(r_n)
                  mx_pv = pv
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
         call ssm_wr_tbu(2, 1, n_r, n_c, 
     &        n_it, vr_in_r, tbu, pr_act, ph_1_rdc_co, 0d0, rao, aa_2)
         call ssm_wr_tbu(3, 1, n_r, n_c, 
     &        n_it, vr_in_r, tbu, pr_act, ph_1_rdc_co, xp_tl, rao, aa_1)
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
      aa = mv_dir*aa
 300  continue
      r_n = r_o_vr_t_lv_bs
      write(*, '(a, i2)')
     &     'Enter row of leaving variable - SSM chooses ', r_n
c      read(5, *, end=301, err=301)r_n
      if (r_n .lt. 0 .or. r_n .gt. n_r) then
         print*, ' Invalid row choice '
         goto 300
      endif
      if (r_n .ne. r_o_vr_t_lv_bs) then
         r_o_vr_t_lv_bs = r_n
         aa = mv_dir*rao(r_o_vr_t_lv_bs)
      endif
 301  continue
      pv = tbu(vr_t_en_bs, r_o_vr_t_lv_bs)
      vr_t_lv_bs = vr_in_r(r_o_vr_t_lv_bs)
      write(ch80_li, '(a, i2, a, i2)')
     &     'Leaving variable is ', vr_t_lv_bs, 
     &     ' in row ', r_o_vr_t_lv_bs
      print*, ch80_li
      if (rp_bs_cg) write(3, '(/a)')ch80_li
      do 310, r_n = 1, n_r
         vr_n = vr_in_r(r_n)
         pr_act(vr_n) = pr_act(vr_n) - aa*tbu(vr_t_en_bs, r_n)
         tbu(0, r_n) = pr_act(vr_n)
 310  continue
      tbu(0, 0) = tbu(0, 0) + aa*tbu(vr_t_en_bs, 0)
      pr_act(vr_t_en_bs) = pr_act(vr_t_en_bs) + aa
      if (vr_t_en_bs .eq. vr_t_lv_bs) then
         at_lb(vr_t_en_bs) = .not. at_lb(vr_t_en_bs)
         goto 349
      else
         at_lb(vr_t_lv_bs) = pr_act(vr_t_lv_bs) .le. 1d-4
      endif
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
 349  continue
      if (ph_1) then
         mu = ph_1_rdc_co(vr_t_en_bs)
         ph_1_rdc_co(0) = ph_1_rdc_co(0) + mv_dir*aa*mu
         do 350, c_n = 1, n_c
            ph_1_rdc_co(c_n) =
     &           ph_1_rdc_co(c_n) - tbu(c_n, r_o_vr_t_lv_bs)*mu
 350     continue
      endif
      vr_in_r(r_o_vr_t_lv_bs) = vr_t_en_bs
      n_it = n_it + 1
c
c     Remove the artificial variables if feasible
c
      if (ph_1 .and. ph_1_rdc_co(0) .ge. -1d-12) then
         do 355, r_n = 1, n_r
            if (vr_in_r(r_n) .gt. ml_n_c + n_r) then
               print*,' Feasible but artificial ',
     &              vr_in_r(r_n)-(ml_n_c + n_r), ' still basic'
               er_fd = .true.
            endif
 355     continue
         if (er_fd) goto 1000
         call ssm_wr_tbu(2, 0, n_r, n_c, 
     &        n_it, vr_in_r, tbu, pr_act, ph_1_rdc_co, 0d0, rao, aa)
         call ssm_wr_tbu(3, 0, n_r, n_c, 
     &        n_it, vr_in_r, tbu, pr_act, ph_1_rdc_co, 0d0, rao, aa)
         ph_1 = .false.
         n_c = ml_n_c + n_r
         do 356, c_n = 0, n_c
            ph_1_rdc_co(c_n) = 0d0
 356     continue
      endif
      if (n_it .le. mx_n_it) goto 100
 1000 continue
      call ssm_wr_tbu(2, 0, n_r, n_c, 
     &     n_it, vr_in_r, tbu, pr_act, ph_1_rdc_co, 0d0, rao, aa)
      call ssm_wr_tbu(3, 0, n_r, n_c, 
     &     n_it, vr_in_r, tbu, pr_act, ph_1_rdc_co, 0d0, rao, aa)
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

      subroutine ssm_wr_tbu(cn_n, n_aux, n_r, n_c, 
     &     n_it, vr_in_r, tbu, pr_act, ph_1_rdc_co, tl, rao, aa)
      implicit none
      include 'SSM.INC'
      integer cn_n, n_aux, n_r, n_c,n_it, vr_in_r(0:n_r)
      double precision tbu(0:mx_n_tbu_c, 0:n_r)
      double precision pr_act(0:n_c), rao(0:n_r)
      double precision ph_1_rdc_co(0:mx_n_tbu_c)
      double precision tl, aa
      character*12 rl_t_TEX_ch12
      integer r_n, ln, c_n, p
      character*1000 tbu_r_li
      character*12 TEX_ch(1:mx_n_tbu_c)
      double precision du_act(0:mx_n_tbu_c)
      character*1 backslash
      logical w_rao

      w_rao = n_aux .gt. 0

      backslash = char(92)
      do 1, c_n = 0, n_c
         du_act(c_n) = 1d3*ph_1_rdc_co(c_n) + tbu(c_n, 0)
 1    continue
      
      if (cn_n .eq. 3) then
         write(cn_n, '(a1, a7/)')backslash, 'medskip'

         p = 29+n_c
         write(tbu_r_li(1:p), 9200)backslash, backslash, 
     &        ('r', c_n = 1, n_c)
         if (w_rao) then
            write(tbu_r_li(p+1:p+18), 9201)'r', backslash, n_c+1
         else
            write(tbu_r_li(p+1:p+18), 9201)' ', backslash, n_c+1
         endif
         p = p+18
         write(cn_n, '(a)')tbu_r_li(1:p)

         p = 14*n_c
         write(tbu_r_li(1:p), 9210)(c_n, c_n = 1, n_c)
         if (w_rao) then
            write(tbu_r_li(p+1:p+39), 9211)'&    Ratio  ', 
     &           backslash, backslash, backslash, n_c+1
         else
            write(tbu_r_li(p+1:p+39), 9211)'            ', 
     &           backslash, backslash, backslash, n_c+1
         endif
         p = p+39
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
               write(tbu_r_li(p+1:p+12), 9222)backslash, n_c+1
               p = p+12
            endif
            write(cn_n, '(a)')tbu_r_li(1:p)
 20      continue
         do 30, c_n = 1, n_c
            TEX_ch(c_n) = rl_t_TEX_ch12(du_act(c_n))
 30      continue
         TEX_ch(n_c+1) = rl_t_TEX_ch12(-du_act(0))
         p = 14*(n_c+1)
         write(tbu_r_li(1:p), 9220)(TEX_ch(c_n), c_n = 1, n_c+1)
         write(tbu_r_li(p+1:p+14), 9223)
     &        backslash, backslash, backslash, n_c+1
         p = p+14
         write(cn_n, '(a)')tbu_r_li(1:p)
c
c     Write out primal activities
c
         do 40, c_n = 1, n_c
            TEX_ch(c_n) = rl_t_TEX_ch12(pr_act(c_n))
 40      continue
         TEX_ch(n_c+1) = rl_t_TEX_ch12(0d0)
         p = 14*(n_c+1)
         write(tbu_r_li(1:p), 9220)(TEX_ch(c_n), c_n = 1, n_c+1)
         write(tbu_r_li(p+1:p+14), 9223)
     &        backslash, backslash, backslash, n_c+1
         p = p+14
         write(cn_n, '(a)')tbu_r_li(1:p)
c
c     End the tabular environment
c
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
         call ssm_g_tbu_r_li(n_aux, 0, n_c, du_act, tl, aa,
     &        tbu_r_li, ln)
         write(cn_n, 9010)tbu_r_li(1:ln)
         call ssm_g_tbu_r_li(0, -n_it, n_c, pr_act, tl, aa,
     &        tbu_r_li, ln)
         write(cn_n, 9010)tbu_r_li(1:ln)
      endif
      return
 9010 format(a)
 9200 format(a1, 'centerline{', a1, 'begin{tabular}{|',20a1)
 9201 format('|r|', a1, 'l}', a1, 'cline{1-', i2, '}')

 9210 format(20(5x, i3, 5x, '&'))
 9211 format('    RHS      ', a12, 3a1, 'cline{1-', i2, '}')

 9220 format(22(1x, a12, '&'))
 9221 format(2a1)
 9222 format( a1, 'cline{1-', i2, '}')
 9223 format(3a1, 'cline{1-', i2, '}')
 9230 format(22(1x, a12, '&'), '\\')
 9240 format(a1, 'end{tabular}}', a1, 'medskip')
      end

      subroutine ssm_g_tbu_hd_li(n_aux, n_c, n_it, tbu_r_li, ln)
      implicit none
      integer n_aux, n_c, n_it, ln
      character*1000 tbu_r_li
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
      character*1000 tbu_r_li
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
      include 'SSM.INC'
      double precision eps
      parameter (eps=1d-8)
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
         if (cn3_tbu_rational) then
            rl_v = rl*(1d0+eps)
            i_v = rl_v
            rl_v = dble(i_v)
            if (abs(rl-rl_v) .lt. eps) then
               write(ch12, 9200)i_v
               goto 1000
            endif
            do 100, den = 2, 10
               rl_v = rl*(1d0+eps)*float(den)
               i_v = rl_v
               rl_v = dble(i_v)/dble(den)
               if (abs(rl-rl_v) .lt. eps) then
                  write(ch12, 9201)i_v, den
                  goto 1000
               endif
 100        continue
         endif
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
c 9000 format('$', 4x, f6.2, '$')
c 9100 format('$', 4x, f6.2, '$')
 9200 format('$',   i10, '$')
 9201 format('$',i7, '/', i2, '$')
      end
