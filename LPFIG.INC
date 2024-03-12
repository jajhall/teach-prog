      integer cm2xfig
      double precision mx_paper_wd
      double precision mx_paper_ht
      double precision lh_axis_margin
      double precision rh_axis_margin
      double precision top_axis_margin
      double precision bottom_axis_margin
      double precision cs_wd
      common/lp_fig_pm_com/mx_paper_wd,
     &     mx_paper_ht,
     &     lh_axis_margin,
     &     rh_axis_margin,
     &     top_axis_margin,
     &     bottom_axis_margin,
     &     cs_wd,
     &     cm2xfig

      integer mx_n_r
      parameter (mx_n_r = 10)
      double precision co(2), mtx(2, mx_n_r), rhs(mx_n_r)
      integer n_r
      common/lpfig_lp_com/co, mtx, rhs, n_r
      double precision xscl, yscl
      double precision lp_cs_wd
      double precision lp_mn_x, lp_mx_x, lp_mn_y, lp_mx_y
      double precision axis_mn_x, axis_mx_x, axis_mn_y, axis_mx_y
      integer x0, y0
      integer fig_mn_x, fig_mx_x, fig_mn_y, fig_mx_y
      integer lp_rd_cn, wr_cn
      logical sans_serif
      common/lpfig_com/xscl, yscl, lp_cs_wd,
     &     lp_mn_x, lp_mx_x, lp_mn_y, lp_mx_y,
     &     axis_mn_x, axis_mx_x, axis_mn_y, axis_mx_y,
     &     x0, y0,
     &     fig_mn_x, fig_mx_x, fig_mn_y, fig_mx_y, 
     &     lp_rd_cn, wr_cn,
     &     sans_serif

