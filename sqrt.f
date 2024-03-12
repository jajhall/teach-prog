      program root
      implicit none
      double precision v0, e0
      double precision v, e, d
      real su_tt, tt(2), su_tt0, tt0(2), etime
      integer k, mx_k

      e0 = 0d0
      print*,'Enter mx_k, v0, d'
      read*, mx_k, v0, d
      e = e0
      v = v0
      su_tt0 = etime(tt0)
      do 10, k = 1, mx_k
         v = v + d
         e = e + sqrt(v)
 10   continue
      su_tt = etime(tt)
      print*, tt(1)-tt0(1), (tt(1)-tt0(1))/mx_k
      e = e0
      v = v0
      su_tt0 = etime(tt0)
      do 20, k = 1, mx_k
         v = v + d
         e = e + v
 20   continue
      su_tt = etime(tt)
      print*, tt(1)-tt0(1), (tt(1)-tt0(1))/mx_k
      stop
      end
