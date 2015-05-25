      PROGRAM test_fbspline
c23456789
      implicit none
      real*8 x0, x1
      integer*8 spline
      integer num, x0code, x1code
      real x0val, x1val
      real a(11), val
      real*8 x

      x0 = 0.0D0
      x1 = 1.0D0
      num = 11
      a(1) = 1.5
      a(2) = -1.3
      a(3) = 2.3
      a(4) = 3.1
      a(5) = 1.8
      a(6) = 0.9
      a(7) = -0.2
      a(8) = -1.0
      a(9) = -1.2
      a(10) = -0.8
      a(11) = 1.5
      x0code = 0
      x1code = 0
     
      call fcreate_ubspline_1d_s (x0, x1, num, x0code, x0val,
     +                            x1code, x1val, a, spline)
      x = 0.0D0
      do while (x .le. 1.00001D0) 
        call feval_ubspline_1d_s (spline, x, val)
        write (*,*) x, val
        x = x + 0.001D0
      enddo

      stop
      end
