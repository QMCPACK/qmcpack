      SUBROUTINE transpose_1(nx,ny,first_x,last_x,input,output)
      integer nx,ny,first_x,last_x
      double complex input(ny,nx),output(nx,ny)
      integer j,k
      !!faster input and slower putput
      DO j=first_x+1,last_x
        DO k=1,ny
          output(j,k)=input(k,j)
        END DO
      END DO
      RETURN
      END

      SUBROUTINE trans_b(nx_th,ny,a,ny_th,nx,b,i0,j0)
      integer nx_th,ny,ny_th,nx
      integer i0,j0
      double complex a(0:ny-1,0:nx_th-1)
      double complex b(0:nx-1,0:ny_th-1)
      integer i,j,ii,jj
      ii=i0
      DO i=0,ny_th-1
        jj=j0
        DO j=0,nx_th-1
        b(jj,i)=a(ii,j)
        jj=jj+1
        END DO
        ii=ii+1
      END DO
      RETURN
      END

      SUBROUTINE transpose_xy(nx,ny,m,j_start,j_end,a,b)
      integer nx,ny,m
      integer j_start, j_end
      double complex a(ny,nx,m), b(nx,ny,m)
      integer i,j,k
      write(*,*) nx, ny, m, j_start,j_end
      DO i=1,m
        DO j=j_start+1,j_end
          DO k=1,nx
            b(k,j,i)=a(j,k,i)
          END DO
        END DO
      END DO
      RETURN
      END

      SUBROUTINE transpose_yx(nx,ny,m,j_start,j_end,a,b)
      integer nx,ny,m
      integer j_start, j_end
      double complex a(nx,ny,m), b(ny,nx,m)
      integer i,j,k
      DO i=1, m
        DO j=j_start+1,j_end
          DO k=1,ny
            b(k,j,i)=a(j,k,i)
          END DO
        END DO
      END DO
      RETURN
      END
