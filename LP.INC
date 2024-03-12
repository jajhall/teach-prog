      integer mx_n_r, mx_n_c
      parameter (mx_n_r = 10, mx_n_c = 10)
      double precision co(mx_n_c), mtx(mx_n_c, mx_n_r), rhs(mx_n_r)
      double precision f0
      integer n_r, n_c, mx_mn, lp_rd_cn
      common/lpfig_lp_com/co, mtx, rhs, f0, n_r, n_c, mx_mn, lp_rd_cn
