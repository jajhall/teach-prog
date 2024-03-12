      implicit none
      double precision v, v_t_n_dpl
      double precision nw_v
      integer n

 10   continue
      print*,' Enter v, n'
      read*, v, n
      if (n .lt. 0) stop

      nw_v = v_t_n_dpl(v, n)
      print*, ' v = ', v, ' is ', nw_v, ' to ', n, ' dpl'
      goto 10
      end
      
      double precision function v_t_n_dpl(v, n)
      implicit none
      double precision v
      integer n
      double precision nw_v
      character*16 ch16
      character*6 ch6_fmt

      ch6_fmt = '(f5.1)'
      ch6_fmt = '(f5.'//char(48+n)//')'
      write(ch16, ch6_fmt)v
      read(ch16, *)nw_v
      v_t_n_dpl = nw_v
      return
      end
         

      
          
