c     Solves max mx_mn*c^Tx
c                  st  A^Tx <= b
c                         x >= 0
      implicit none
      integer mx_n_c, mx_n_r
      parameter (mx_n_c = 10, mx_n_r = 2)
      double precision tbu(0:mx_n_c+mx_n_r, 0:mx_n_r)
      double precision pr_act(0:mx_n_c+mx_n_r)
      double precision bs_mtx(mx_n_r, mx_n_r)
      integer non_bc1, non_bc2, non_bc3, bc1, bc2
      double precision det, obj, mx_obj
      logical fs
      logical non_bc_vr(0:mx_n_c+mx_n_r)
      integer n_c, n_r
      integer c_n, r_n
      integer vtx_n, opt_vtx_n
      double precision mx_mn
      character*64 filename

      print*,' Enter Filename'
      read*, filename
      open (unit = 1, file = filename)
      read(1,*)n_r, n_c, mx_mn
      if (n_r .gt. mx_n_r .or. n_c .gt. mx_n_c .or.
     &     abs(mx_mn) .ne. 1d0) then
         if (n_r .gt. mx_n_r) 
     &        print*, 'n_r .gt. mx_n_r', n_r, mx_n_r
         if (n_c .gt. mx_n_c) 
     &        print*, 'n_c .gt. mx_n_c', n_c, mx_n_c
         if (abs(mx_mn) .ne. 1d0) 
     &        print*, 'abs(mx_mn) .ne. 1d0', abs(mx_mn)
         stop
      endif
      read(1,*) (tbu(c_n,   0), c_n = 1, n_c)
      read(1,*)((tbu(c_n, r_n), c_n = 1, n_c), tbu(0, r_n),
     &     r_n = 1, n_r)

      do 6, c_n = 1, n_r
         do 5, r_n = 1, n_r
            tbu(n_c+c_n, r_n) = 0d0
 5       continue
         tbu(n_c+c_n, c_n) = 1d0
 6    continue
      do 10, c_n = 1, n_r+n_c
         non_bc_vr(c_n) = .false.
 10   continue
      vtx_n = 0
      mx_obj = -1d31
      opt_vtx_n = 0
      write(2, 9000)
      do 90, non_bc1 = 1, n_r+n_c
         do 80, non_bc2 = non_bc1+1, n_r+n_c
            do 70, non_bc3 = non_bc2+1, n_r+n_c
               vtx_n = vtx_n + 1
               non_bc_vr(non_bc1) = .true.
               non_bc_vr(non_bc2) = .true.
               non_bc_vr(non_bc3) = .true.
               pr_act(non_bc1) = 0d0
               pr_act(non_bc2) = 0d0
               pr_act(non_bc3) = 0d0
               do 20, bc1 = 1, n_r+n_c
                  if (.not.non_bc_vr(bc1)) goto 25
 20            continue
 25            continue
               do 30, bc2 = 1, n_r+n_c
                  if (.not.non_bc_vr(bc2) .and. bc2 .ne. bc1) goto 35
 30            continue
 35            continue
               do 40, r_n = 1, n_r
                  bs_mtx(1, r_n) = tbu(bc1, r_n)
                  bs_mtx(2, r_n) = tbu(bc2, r_n)
 40            continue
               det = bs_mtx(1, 1)*bs_mtx(2, 2)-bs_mtx(1, 2)*bs_mtx(2, 1)
               if (det .eq. 0d0) then
                  write(2, 9100)vtx_n, 
     &                 non_bc1, non_bc2, non_bc3, bc1, bc2
               else
                  pr_act(bc1) = 
     &                 (bs_mtx(2, 2)*tbu(0, 1)-
     &                 bs_mtx(2, 1)*tbu(0, 2))/det
                  pr_act(bc2) = 
     &                 (-bs_mtx(1, 2)*tbu(0, 1)+
     &                 bs_mtx(1, 1)*tbu(0, 2))/det
                  fs = pr_act(bc1) .ge. 0d0 .and.
     &                 pr_act(bc2) .ge. 0d0
                  obj = 0d0
                  do 50, c_n = 1, n_c
                     obj = obj + pr_act(c_n)*tbu(c_n, 0)
 50               continue
                  if (fs .and. obj .gt. mx_obj) then
                     opt_vtx_n = vtx_n
                     mx_obj = obj
                  endif
                  if (fs) then
                     write(2, 9100)vtx_n, 
     &                    non_bc1, non_bc2, non_bc3, bc1, bc2,
     &                    (pr_act(c_n), c_n = 1, n_c+n_r),
     &                    obj
                  else
                     write(2, 9101)vtx_n, 
     &                    non_bc1, non_bc2, non_bc3, bc1, bc2,
     &                    (pr_act(c_n), c_n = 1, n_c+n_r)
                  endif
               endif
               non_bc_vr(non_bc3) = .false.
               non_bc_vr(non_bc2) = .false.
               non_bc_vr(non_bc1) = .false.
 70         continue
 80      continue
 90   continue
      write(2, *)'Optimal vertex ', opt_vtx_n, ' has objective ', mx_obj
      stop
 9000 format(' ')
 9100 format(6i4, 5(2x, g11.4), 2x, g11.4, 2x, 'Fs')
 9101 format(6i4, 5(2x, g11.4), 15x, 'Ifs')
      end
