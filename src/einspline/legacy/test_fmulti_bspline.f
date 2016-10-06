!     //////////////////////////////////////////////////////////////////////////////////////
!     // This file is distributed under the University of Illinois/NCSA Open Source License.
!     // See LICENSE file in top directory for details.
!     //
!     // Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
!     //
!     // File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign 
!     //                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
!     //
!     // File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
!     //////////////////////////////////////////////////////////////////////////////////////


      PROGRAM test_fbspline
c23456789
      implicit none
      real*8 x0, x1
      integer*8 spline
      integer num, x0code, x1code, num_splines
      real x0val, x1val
      real a(11), b(11), val(2)
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

      b( 1) = 3.1
      b( 2) = -1.7
      b( 3) = 2.1
      b( 4) = 6.7
      b( 5) = 1.4
      b( 6) = 3.9
      b( 7) = -0.3
      b( 8) = -1.5
      b( 9) = -9.1
      b(10) = 5.1
      b(11) = 3.1
      x0code = 0
      x1code = 0
      num_splines = 2
     
      call fcreate_multi_ubspline_1d_s (x0, x1, num, x0code, x0val,
     +     x1code, x1val, num_splines, spline)
      call fset_multi_ubspline_1d_s (spline, 0, a);
      call fset_multi_ubspline_1d_s (spline, 1, b);
      x = 0.0D0
      do while (x .le. 1.00001D0) 
        call feval_multi_ubspline_1d_s (spline, x, val)
        write (*,*) x, val(1), val(2)
        x = x + 0.001D0
      enddo

      stop
      end
