!//////////////////////////////////////////////////////////////////////////////////////
!// This file is distributed under the University of Illinois/NCSA Open Source License.
!// See LICENSE file in top directory for details.
!//
!// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
!//
!// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
!//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
!//
!// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
!//////////////////////////////////////////////////////////////////////////////////////
    
    



         subroutine makefun(iopt,iocc,indt,indtm,indpar,indorb,indshell
     1,nelskip,z,dd,zeta,r,rmu,distp,kion,nion,rion,iflagnorm,c) 
         implicit none 
         integer iopt,indt,i,k,nelskip,indpar,iocc(*),indorbp
     1,indorb,indshell,indshellp,ic,indtmax,kion(*),k_1
     1,nion,skip,iflagnorm,ngauss,indparp,indtm,npower
         real*8 z(nelskip,*),dd(*),zeta(*),rmu(3,0:indt,*),r(0:indt,*)
     1,distp(0:(indt+1)*nion,*),peff,fun,fun0,fun2,rp1,rp2,dd1,dd2
     1,dd3,dd4,dd5,cost1d,cost2d,cost3d,rion(3,*),c,pi,dsqrt,funp,fun2p
     1,peff2,arg,c0,c1,cost1f,cost2f,cost3f,cost4f,cost
!     3d without cusp condition
!     cost1d=0.5 cost2d=dsqrt(3.d0)*cost1d cost3d=2.d0*cost2d
      parameter(cost1d=0.5d0,cost2d=0.8660254038d0,cost3d=1.732050808d0)
!     parameter(cost1d=1.d0,cost2d=1.d0,cost3d=1.d0)
!     4f norm coeff
!     cost1f=0.5 cost2f=dsqrt(6.d0)/2.d0*cost1f cost3f=dsqrt(15.d0)*cost1f
!     cost4f=dsqrt(10.d0)/2.d0*cost1f
      parameter(cost1f=0.5d0,cost2f=0.6123724357d0,cost3f=1.936491673d0
     &,cost4f=0.790569415d0)
      parameter(pi=3.141592653589793d0)
!
!        indorb are the number of orbitals occupied before calling 
!        this subroutine 
!
!        indpar is the number of variational parameters used 
!        before calling this subroutine 
!
!        indshell is the index of the last  occupied orbital 
!        in the shell, characterized by occupation number iocc(indshell+1)...
!        
!        z(i,indt+4)   contains the laplacian of the orbital i
!        z(i,indt+mu) contains the gradient of the orbital i (mu=1,2,3)

      indtmax=0
      if(indt.eq.1) indtmax=1
      if(indt.eq.0) then 
              skip=1
      else 
              skip=indt+1
      endif

      select case(iopt)
      
      case(1)        ! 1s single Z NO CUSP!

         k_1=kion(1)

         indshellp=indshell+1
         

         if(iocc(indshellp).eq.1) then

            dd1=dd(indpar+1)
            if(iflagnorm.gt.2) then  
            c=dd1*dsqrt(dd1)/dsqrt(pi)
            endif

            indorbp=indorb+1
               do k=indtmax,indtm
                  distp(k,1)=c*dexp(-dd1*r(k,k_1))
               enddo

            do i=1,indtm
               z(indorbp,i)=distp(i,1)
            enddo

            if(indt.ne.1) then
               fun=-dd1*distp(0,1)


               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)/r(0,k_1)
                  enddo

                 z(indorbp,indt+4)=(-2.d0*dd1/r(0,k_1)+dd1**2)
     1            *distp(0,1)
        
               else

                  do i=1,4
                     z(indorbp,indt+i)=0.d0
                  enddo
               endif

            endif

            indorb=indorbp
            
         endif
         
         indpar=indpar+1
         indshell=indshellp


         case(2)     ! 1s double Z with cusp cond
!

         k_1=kion(1)

         indshellp=indshell+1

         if(iocc(indshellp).eq.1) then

            indorbp=indorb+1

            dd1=dd(indpar+1)
            dd2=dd(indpar+2)
            peff=(zeta(k_1)-dd1)/(dd2-zeta(k_1))

            if(iflagnorm.gt.2) then 
       c=1.d0/2.d0/dsqrt(2.d0*pi*(1.d0/(2.d0*dd1)**3
     &    +2.d0*peff/(dd1+dd2)**3+peff**2/(2.d0*dd2)**3))
            endif 

               do k=indtmax,indtm
                distp(k,1)=c*dexp(-dd1*r(k,k_1))
                distp(k,2)=c*dexp(-dd2*r(k,k_1))
              enddo

            do i=1,indtm
               z(indorbp,i)=distp(i,1)+peff*distp(i,2)
            enddo

            if(indt.ne.1) then
                fun=(-dd1*distp(0,1)-dd2*distp(0,2)*peff)/r(0,k_1)

               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

                 z(indorbp,indt+4)=(-2.d0*dd1/r(0,k_1)+dd1**2)
     1            *distp(0,1)+peff*(-2.d0*dd2/r(0,k_1)+dd2**2)
     1            *distp(0,2)
            
               else

                  do i=1,4
                     z(indorbp,indt+i)=0.d0
                  enddo
               endif

            endif

            indorb=indorbp
            
         endif
         
         indpar=indpar+2
         indshell=indshellp
 

         case(3)      ! 1s double Z NO CUSP
c

         k_1=kion(1)

         indshellp=indshell+1

         if(iocc(indshellp).eq.1) then

            indorbp=indorb+1

            dd1=dd(indpar+1)
            dd2=dd(indpar+2)
            peff=dd(indpar+3)

            if(iflagnorm.gt.2) then 
        c=1.d0/2.d0/dsqrt(2.d0*pi*(1.d0/(2.d0*dd1)**3
     &    +2.d0*peff/(dd1+dd2)**3+peff**2/(2.d0*dd2)**3))
            endif 

            do i=indpar+1,indpar+2
               do k=indtmax,indtm
                  distp(k,i-indpar)=c*dexp(-dd(i)*r(k,k_1))
               enddo
            enddo

            do i=1,indtm
               z(indorbp,i)=distp(i,1)+peff*distp(i,2)
            enddo

            if(indt.ne.1) then
               fun=-dd1*distp(0,1)-peff*dd2*distp(0,2)

               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)/r(0,k_1)
                  enddo

                 z(indorbp,indt+4)=(-2.d0*dd1/r(0,k_1)+dd1**2)
     1            *distp(0,1)+peff*(-2.d0*dd2/r(0,k_1)+dd2**2)
     1            *distp(0,2)
            
               else

                  do i=1,4
                     z(indorbp,indt+i)=0.d0
                  enddo
               endif

            endif

            indorb=indorbp
            
         endif
         
         indpar=indpar+3
         indshell=indshellp

          case(4)     ! 2s 2pz Hybryd single Z 
                                ! normalized
         k_1=kion(1)
         indshellp=indshell+1


         if(iocc(indshellp).eq.1) then 

            indorbp=indorb+1
            dd1=dd(indpar+1)
            dd2=dd(indpar+2)

            if(iflagnorm.gt.2) then 
            c=dd1**2.5d0/dsqrt(3.d0*pi*(1.d0+dd2**2/3.d0)) 
            endif 
            
            do k=indtmax,indtm
            distp(k,1)=c*dexp(-dd1*r(k,k_1))
            enddo


            do i=1,indtm
            z(indorbp,i)=(r(i,k_1)+dd2*rmu(3,i,k_1))*distp(i,1)
            enddo            

            if(indt.ne.1) then
      
             fun=distp(0,1)*(1.d0-dd1*r(0,k_1))
             funp=-dd2*dd1*distp(0,1)*rmu(3,0,k_1)
             fun2=distp(0,1)*(dd1**2*r(0,k_1)-2.d0*dd1)
             fun2p=dd1**2*dd2*distp(0,1)*rmu(3,0,k_1)

               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
      z(indorbp,indt+i)=(fun+funp)*rmu(i,0,k_1)/r(0,k_1)
                  enddo
                  z(indorbp,indt+3)=z(indorbp,indt+3)+dd2*distp(0,1)
             z(indorbp,indt+4)=(2.d0*fun+4.d0*funp)/r(0,k_1)
     1         +(fun2+fun2p)
               else

                  do i=1,4
                     z(indorbp,indt+i)=0.d0 
                  enddo
               endif

            endif

            indorb=indorbp

         endif
         indpar=indpar+2
         indshell=indshellp
  
      
          case(5)     ! 2s single Z NO CUSP
                                ! normalized
         k_1=kion(1)
         indshellp=indshell+1


         if(iocc(indshellp).eq.1) then 

            indorbp=indorb+1
            dd1=dd(indpar+1)
            
            do k=indtmax,indtm
               distp(k,1)=dexp(-dd1*r(k,k_1))
            enddo

            if(iflagnorm.gt.2) then 
            c=dd1**2.5d0/dsqrt(3.d0*pi) 
            endif 

            do i=1,indtm
            z(indorbp,i)=c*r(i,k_1)*distp(i,1)
            enddo            

            if(indt.ne.1) then
      
             fun=distp(0,1)*(1.d0-dd1*r(0,k_1))
             fun2=distp(0,1)*(dd1**2*r(0,k_1)-2.d0*dd1)

               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
                     z(indorbp,indt+i)=c*fun*rmu(i,0,k_1)/r(0,k_1)
                  enddo
                  z(indorbp,indt+4)=c*2.d0*fun/r(0,k_1)+c*fun2
               else

                  do i=1,4
                     z(indorbp,indt+i)=0.d0 
                  enddo
               endif

            endif

            indorb=indorbp

         endif
         indpar=indpar+1
         indshell=indshellp
  

       case(6)     ! 2s double  Z NO CUSP
                                ! normalized
         k_1=kion(1)
         indshellp=indshell+1


         if(iocc(indshellp).eq.1) then 

            indorbp=indorb+1
            dd1=dd(indpar+1)
            dd2=dd(indpar+2)
            peff=dd(indpar+3)
            
            if(iflagnorm.gt.2) then 
c            c=                                              WRONG
c    &0.5*dsqrt((1.d0/dd1**5/32.d0+2.d0*peff/(dd1+dd2)**5
c    &+peff**2/dd2**5/32.d0)/(24.d0*pi))

            c=1.d0/dsqrt((3.d0*pi)*
     & (1.d0/dd1**5+ 64.d0*peff/(dd1+dd2)**5+peff**2/dd2**5)) 

            endif 

            do k=indtmax,indtm
            distp(k,1)=c*dexp(-dd1*r(k,k_1))
            distp(k,2)=c*dexp(-dd2*r(k,k_1))
            enddo

            do i=1,indtm
            z(indorbp,i)=r(i,k_1)*(distp(i,1)+distp(i,2)*peff)
            enddo            

            if(indt.ne.1) then
      
             fun=distp(0,1)*(1.d0-dd1*r(0,k_1))
     1           +peff*distp(0,2)*(1.d0-dd2*r(0,k_1))
             fun2=distp(0,1)*(dd1**2*r(0,k_1)-2.d0*dd1)+peff*distp(0,2)
     1 *(dd2**2*r(0,k_1)-2.d0*dd2)

               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)/r(0,k_1)
                  enddo
                  z(indorbp,indt+4)=(2.d0*fun/r(0,k_1)+fun2)
               else

                  do i=1,4
                     z(indorbp,indt+i)=0.d0 
                  enddo
               endif

            endif

            indorb=indorbp

         endif
         indpar=indpar+3
         indshell=indshellp

      case(7)     ! 2s double Z NO CUSP
                                ! normalized IS WRONG!!!
         k_1=kion(1)
         indshellp=indshell+1


         if(iocc(indshellp).eq.1) then 

            indorbp=indorb+1
            dd1=dd(indpar+1)
            dd2=dd(indpar+2)
            peff=dd(indpar+3)
            
            do k=indtmax,indtm
               distp(k,1)=dexp(-dd1*r(k,k_1))
               distp(k,2)=dexp(-dd2*r(k,k_1))
            enddo
            if(iflagnorm.gt.2) then 
            c=
     &1/dsqrt(1/(3.D0/4.D0/dd1**5+peff**2/dd2**3/4+12*peff/
     & (dd1+dd2)**4))*1.d0/dsqrt(4.0*pi)
            endif 
            
            do i=1,indtm
            z(indorbp,i)=c*(distp(i,1)+r(i,k_1)*distp(i,2)*peff)
            enddo            

            if(indt.ne.1) then
      
             fun=-dd1*distp(0,1)+peff*distp(0,2)*(1.d0-dd2*r(0,k_1))
             fun2=distp(0,1)*dd1**2
     &+peff*distp(0,2)*(dd2**2*r(0,k_1)-2.d0*dd2)

               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
                     z(indorbp,indt+i)=fun*c*rmu(i,0,k_1)/r(0,k_1)
                  enddo
                  z(indorbp,indt+4)=c*(2.d0*fun/r(0,k_1)+fun2)
               else

                  do i=1,4
                  z(indorbp,indt+i)=0.d0 
                  enddo
               endif

            endif

            indorb=indorbp

         endif
         indpar=indpar+3
         indshell=indshellp



          case(8)     ! 2s double Z WITH CUSP
                                    ! normalized
! exp(-dd1*r) + (dd1-zeta) * r * exp(-dd2*r)

         k_1=kion(1)
         indshellp=indshell+1


         if(iocc(indshellp).eq.1) then 

            indorbp=indorb+1
            dd1=dd(indpar+1)
            dd2=dd(indpar+2)
            peff=dd1-zeta(k_1)
            
            do k=indtmax,indtm
               distp(k,1)=dexp(-dd1*r(k,k_1))
               distp(k,2)=dexp(-dd2*r(k,k_1))
            enddo

            if(iflagnorm.gt.2) then 
            c=1.d0/dsqrt(1/4.d0/dd1**3+12*peff/(dd1+dd2)**4+
     &3*peff**2/4/dd2**5)/dsqrt(4.0*pi)
            endif 
             
            do i=1,indtm
            z(indorbp,i)=c*(distp(i,1)+r(i,k_1)*distp(i,2)*peff)
            enddo            

            if(indt.ne.1) then
      
             fun=-dd1*distp(0,1)+peff*distp(0,2)*(1.d0-dd2*r(0,k_1))
             fun2=distp(0,1)*dd1**2
     &+peff*distp(0,2)*(dd2**2*r(0,k_1)-2.d0*dd2)

               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
                  z(indorbp,indt+i)=fun*c*rmu(i,0,k_1)/r(0,k_1)
                  enddo
                  z(indorbp,indt+4)=c*(2.d0*fun/r(0,k_1)+fun2)
               else

                  do i=1,4
                     z(indorbp,indt+i)=0.d0 
                  enddo
               endif

            endif

            indorb=indorbp

         endif
         indpar=indpar+2
         indshell=indshellp

      case(10)      ! 3s single zeta
c R(r)=r**2*exp(-z1*r)

      
         k_1=kion(1)   
         indshellp=indshell+1

         if(iocc(indshellp).eq.1) then 
                 
            indorbp=indorb+1
            dd1=dd(indpar+1)
            if(iflagnorm.gt.2) then 
            c=dsqrt((2*dd1)**7/720.d0/pi)/2.d0
            endif 
            
            do k=indtmax,indtm
            distp(k,1)=c*dexp(-dd1*r(k,k_1))
            enddo

            do i=1,indtm
            z(indorbp,i)=distp(i,1)*r(i,k_1)**2
            enddo            

            if(indt.ne.1) then
               fun=(2.d0-dd1*r(0,k_1))*distp(0,1)
               fun2=(2.d0-4*dd1*r(0,k_1)+(dd1*r(0,k_1))**2)*distp(0,1)
                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo
                  z(indorbp,indt+4)=2.d0*fun+fun2
            endif

            indorb=indorbp

         endif
         indpar=indpar+1
         indshell=indshellp         


      case(11)     ! 3s double zeta
c R(r)=r**2*(exp(-z1*r)+p*exp(-z2*r))

      
         k_1=kion(1)   
         indshellp=indshell+1

         if(iocc(indshellp).eq.1) then 
                 
            indorbp=indorb+1
            dd1=dd(indpar+1)
            dd2=dd(indpar+2)
            peff=dd(indpar+3)
            if(iflagnorm.gt.2) then 
      c=1.d0/2.d0/dsqrt(pi*720.d0*(1.d0/(2.d0*dd1)**7+
     &      2.d0*peff/(dd1+dd2)**7+peff**2/(2.d0*dd2)**7))
            endif

            do k=indtmax,indtm
            distp(k,1)=c*dexp(-dd1*r(k,k_1))
            distp(k,2)=c*dexp(-dd2*r(k,k_1))
            enddo

            do i=1,indtm
            z(indorbp,i)=(distp(i,1)+peff*distp(i,2))*r(i,k_1)**2
            enddo            

            if(indt.ne.1) then
               rp1=r(0,k_1)**2
c              the first derivative 
               fun=distp(0,1)*(2.d0*r(0,k_1)-dd1*rp1)
     &          +peff*distp(0,2)*(2.d0*r(0,k_1)-dd2*rp1)
c
c              the second derivative
               fun2=distp(0,1)*(2.d0-4.d0*dd1*r(0,k_1)+dd1**2*rp1)
     &  +peff*distp(0,2)*(2.d0-4.d0*dd2*r(0,k_1)+dd2**2*rp1)
c
               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)/r(0,k_1)
                  enddo

                  z(indorbp,indt+4)=2.d0*fun/r(0,k_1)+fun2
               else
                  do i=1,4
                     z(indorbp,indt+i)=0.d0 
                  enddo
               endif

            endif

            indorb=indorbp

         endif
         indpar=indpar+3
         indshell=indshellp         


      case(12)      ! 4s single zeta
c R(r)=r**3*exp(-z1*r)
c         
         indshellp=indshell+1
         k_1=kion(1)


        if(iocc(indshellp).eq.1) then 

            indorbp=indorb+1
            dd1=dd(indpar+1)
            if(iflagnorm.gt.2) then 
            c=dsqrt((2*dd1)**9/40320.d0/pi)/2.d0
            endif 

            do k=indtmax,indtm
            distp(k,1)=c*dexp(-dd1*r(k,k_1))
            enddo

            do i=1,indtm
            z(indorbp,i)=distp(i,1)*r(i,k_1)**3
            enddo            

            if(indt.ne.1) then
               rp1=r(0,k_1)**3
               rp2=r(0,k_1)**2
c
cc              the first derivative 
               fun=distp(0,1)*(3.d0*rp2-dd1*rp1)
cc
cc              the second derivative
               fun2=distp(0,1)*(6.d0*r(0,k_1)-6.d0*dd1*rp2+dd1**2*rp1)
cc
               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)/r(0,k_1)
                  enddo

                  z(indorbp,indt+4)=2.d0*fun/r(0,k_1)+fun2
               else
                  do i=1,4
                     z(indorbp,indt+i)=0.d0 
                  enddo
               endif

            endif
c
            indorb=indorbp
c
         endif
         indpar=indpar+1
         indshell=indshellp         
c


      case(13)      ! 4s double zeta
c R(r)=r**3*(exp(-z1*r)+z3*exp(-z2*r))
c         
         indshellp=indshell+1
         k_1=kion(1)
c
c
         if(iocc(indshellp).eq.1) then 
c
            indorbp=indorb+1
            dd1=dd(indpar+1)
            dd2=dd(indpar+2)
            dd3=dd(indpar+3)
            if(iflagnorm.gt.2) then 
      c=1.d0/2.d0/dsqrt(pi*40320.d0*(1.d0/(2.d0*dd1)**9+
     &        2.d0*dd3/(dd1+dd2)**9+dd3**2/(2.d0*dd2)**9))
            endif 

c
            do k=indtmax,indtm
            distp(k,1)=c*dexp(-dd1*r(k,k_1))
            distp(k,2)=c*dexp(-dd2*r(k,k_1))
            enddo

            do i=1,indtm
            z(indorbp,i)=(distp(i,1)+dd3*distp(i,2))*r(i,k_1)**3
            enddo            
c
            if(indt.ne.1) then
               rp1=r(0,k_1)**3
               rp2=r(0,k_1)**2
c
cc              the first derivative 
               fun=distp(0,1)*(3.d0*rp2-dd1*rp1)
     1       +dd3*distp(0,2)*(3.d0*rp2-dd2*rp1)
cc
c              the second derivative
               fun2=distp(0,1)*(6.d0*r(0,k_1)-6.d0*dd1*rp2+dd1**2*rp1)
     1         +dd3*distp(0,2)*(6.d0*r(0,k_1)-6.d0*dd2*rp2+dd2**2*rp1)
cc
               if(r(0,k_1).gt.1d-9) then
                 do i=1,3
                    z(indorbp,indt+i)=fun*rmu(i,0,k_1)/r(0,k_1)
                  enddo
c
                  z(indorbp,indt+4)=2.d0*fun/r(0,k_1)+fun2
               else
                  do i=1,4
                     z(indorbp,indt+i)=0.d0 
                 enddo
              endif
c
            endif

            indorb=indorbp

         endif
         indpar=indpar+3
         indshell=indshellp         

         case(14)     ! 1s single Z pseudo
! (1.d0 + dd1 r) * exp(-dd1 * r)   ! normalized

         k_1=kion(1)
         indshellp=indshell+1


         if(iocc(indshellp).eq.1) then 

            indorbp=indorb+1
            dd1=dd(indpar+1)
            
            do k=indtmax,indtm
               distp(k,1)=dexp(-dd1*r(k,k_1))
            enddo

            if(iflagnorm.gt.2) then 
            c=dsqrt(dd1**3.d0/7.d0/pi) 
            endif

            do i=1,indtm
               z(indorbp,i)=c*(1.d0+dd1*r(i,k_1))*distp(i,1)
            enddo            

            if(indt.ne.1) then
      
             fun=-distp(0,1)*dd1**2*r(0,k_1)
             fun2=-distp(0,1)*dd1**2*(1.d0-dd1*r(0,k_1))

               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
                     z(indorbp,indt+i)=c*fun*rmu(i,0,k_1)/r(0,k_1)
                  enddo
                  z(indorbp,indt+4)=c*2.d0*fun/r(0,k_1)+c*fun2
               else

                  do i=1,4
                     z(indorbp,indt+i)=0.d0 
                  enddo
               endif

            endif

            indorb=indorbp

         endif
         indpar=indpar+1
         indshell=indshellp
  


         case(15)     ! 1s single Z pseudo
! (r**2 + dd2*(1 + dd1*r))*exp(-dd1*r)  ! normalized

         k_1=kion(1)
         indshellp=indshell+1

         if(iocc(indshellp).eq.1) then 

            indorbp=indorb+1
            dd1=dd(indpar+1)
            dd2=dd(indpar+2)
            
            if(iflagnorm.gt.2) then 
            c=dsqrt(2.d0*dd1**7/pi/
     &      (45.d0+42.d0*dd1**2*dd2+14.d0*dd1**4*dd2**2))
            endif

            do k=indtmax,indtm
               distp(k,1)=c*dexp(-dd1*r(k,k_1))
            enddo


            do i=1,indtm
               z(indorbp,i)=(r(i,k_1)**2+dd2*(1.d0+dd1*r(i,k_1)))
     &                      *distp(i,1)
            enddo            

            if(indt.ne.1) then
      
             fun=distp(0,1)*r(0,k_1)*(2.d0-dd1**2*dd2-dd1*r(0,k_1))
             fun2=distp(0,1)*((1.d0-dd1*r(0,k_1))
     &            *(3.d0-dd1**2*dd2-dd1*r(0,k_1))-1.d0)

               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)/r(0,k_1)
                  enddo
                  z(indorbp,indt+4)=2.d0*fun/r(0,k_1)+fun2
               else

                  do i=1,4
                     z(indorbp,indt+i)=0.d0 
                  enddo
               endif

            endif

            indorb=indorbp

         endif
         indpar=indpar+2
         indshell=indshellp

      case(16)  ! 2s gaussian for pseudo
! R(r)=exp(-z*r**2) single zeta

      
         k_1=kion(1)   
         indshellp=indshell+1

         if(iocc(indshellp).eq.1) then 
                 
            indorbp=indorb+1

            dd(indpar+1)=abs(dd(indpar+1))
            dd1=dd(indpar+1)

            if(iflagnorm.gt.2) then
            c=(2.d0*dd1/pi)**(3.d0/4.d0)
            endif

            do k=indtmax,indtm
            distp(k,1)=c*dexp(-dd1*r(k,k_1)**2)
            enddo

            do i=1,indtm
            z(indorbp,i)=distp(i,1)
            enddo            

            if(indt.ne.1) then
!              the first derivative /r  
            fun=-2.d0*dd1*distp(0,1)

!              the second derivative
            fun2=fun*(1.d0-2.d0*dd1*r(0,k_1)**2)

                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

                  z(indorbp,indt+4)=2.d0*fun+fun2

            endif

            indorb=indorbp

         endif

         indpar=indpar+1
         indshell=indshellp


      case(17)  ! 2s gaussian for pseudo
! R(r)=r**2*exp(-z*r**2) single zeta

      
         k_1=kion(1)   
         indshellp=indshell+1

         if(iocc(indshellp).eq.1) then 
                 
            indorbp=indorb+1
            dd1=dd(indpar+1)

            if(iflagnorm.gt.2) then 
         c=4.d0*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)/dsqrt(15.d0)
            endif 
            
            do k=indtmax,indtm
            distp(k,1)=c*dexp(-dd1*r(k,k_1)**2)
            enddo

            do i=1,indtm
            z(indorbp,i)=distp(i,1)*r(i,k_1)**2
            enddo            

            if(indt.ne.1) then
               rp1=r(0,k_1)**2
!              the first derivative / r 
            fun=2.d0*distp(0,1)*(1.d0-dd1*rp1)
!              the second derivative 
            fun2=2.d0*distp(0,1)*(1.d0-5.d0*dd1*rp1+2.d0*dd1**2*rp1**2)
                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo
                  z(indorbp,indt+4)=2.d0*fun+fun2
            endif

            indorb=indorbp

         endif

         indpar=indpar+1
         indshell=indshellp

      case(18)  ! 2s gaussian for pseudo
! R(r)=r**4*exp(-z*r**2) single zeta

      
         k_1=kion(1)   
         indshellp=indshell+1

         if(iocc(indshellp).eq.1) then 
                 
            indorbp=indorb+1
            dd1=dd(indpar+1)

            if(iflagnorm.gt.2) then
            c=(2.d0*dd1**11/pi)**(1.d0/4.d0)*(512.d0/945.d0/pi)
            endif

            do k=indtmax,indtm
            distp(k,1)=c*dexp(-dd1*r(k,k_1)**2)
            enddo

            do i=1,indtm
            z(indorbp,i)=r(i,k_1)**4*distp(i,1)
            enddo            

            if(indt.ne.1) then
               rp1=r(0,k_1)**2

!              the first derivative 
            fun=distp(0,1)*rp1*(4.d0-2.d0*dd1*rp1)

!              the second derivative
            fun2=distp(0,1)*rp1*(12.d0-18.d0*dd1*rp1
     &           +4.d0*dd1**2*rp1**2)

                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

                  z(indorbp,indt+4)=2.d0*fun+fun2

            endif

            indorb=indorbp

         endif
         indpar=indpar+1
         indshell=indshellp

      case(19)  ! derivative of 16 with respect to z
! R(r)=c*exp(-z*r**2)*(3/4/z-r**2)


      
         k_1=kion(1)   
         indshellp=indshell+1

         if(iocc(indshellp).eq.1) then 
                 
            indorbp=indorb+1

            dd(indpar+1)=abs(dd(indpar+1))
            dd1=dd(indpar+1)

          if(iflagnorm.gt.2) then
            c=(2.d0*dd1/pi)**(3.d0/4.d0)
          endif

            do k=indtmax,indtm
            distp(k,1)=c*dexp(-dd1*r(k,k_1)**2)
            enddo

            do i=1,indtm
            z(indorbp,i)=distp(i,1)*(3.d0/4.d0/dd1-r(i,k_1)**2)
            enddo            

            if(indt.ne.1) then
!              the first derivative /r  
            fun=distp(0,1)*(2.d0*dd1*r(0,k_1)**2-7.d0/2.d0)

!              the second derivative
            fun2=distp(0,1)*(-4.d0*dd1**2*r(0,k_1)**4
     &      +13.d0*dd1*r(0,k_1)**2-7.d0/2.d0)

                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

                  z(indorbp,indt+4)=2.d0*fun+fun2

            endif

            indorb=indorbp

         endif

         indpar=indpar+1
         indshell=indshellp



      case(20)      ! 2p single zeta 
c     2p single Z with no  cusp condition

         k_1=kion(1)
         dd1=dd(indpar+1)
         if(iflagnorm.gt.2) then 
         c=dsqrt((2.d0*dd1)**5/8.d0/pi)/2.d0
         endif 


         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1))
         enddo
         
         indorbp=indorb
c         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
                  z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)
              enddo
            endif
         enddo

         if(indt.ne.1) then
            fun0=distp(0,1)
            fun=-dd1*distp(0,1)
            fun2=dd1**2*distp(0,1)

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun/r(0,k_1)
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
             z(indorbp,indt+4)=rmu(ic,0,k_1)*(4.d0*fun/r(0,k_1)+fun2)
                  endif
               enddo

            else
               
               indorbp=indorb 
c
               do ic=1,3 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,4
                        z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo

            endif  !endif for r(0)

         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp


      case(21)     ! 2p double zeta
c     2p without cusp condition

         k_1=kion(1)

         dd1=dd(indpar+1)
         dd2=dd(indpar+2)
         peff=dd(indpar+3)
         if(iflagnorm.gt.2) then 
         c=0.5d0/dsqrt(8.d0*pi*(1.d0/(2.d0*dd1)**5
     &   +2.d0*peff/(dd1+dd2)**5+peff**2/(2.d0*dd2)**5))
         endif
         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1))
         distp(k,2)=c*dexp(-dd2*r(k,k_1))
         enddo

         do i=indtmax,indtm
         distp(i,3)=distp(i,1)+peff*distp(i,2)
         enddo 
         
         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,3)
               enddo
            endif
         enddo


         if(indt.ne.1) then
       fun=(-dd1*distp(0,1)-dd2*peff*distp(0,2))/r(0,k_1)
       fun2=dd1**2*distp(0,1)+peff*dd2**2*distp(0,2)

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                  z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*fun
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+distp(0,3)
                     enddo
            z(indorbp,indt+4)=rmu(ic,0,k_1)*(4.d0*fun+fun2)
                  endif
               enddo

            else
               
               indorbp=indorb 

               do ic=1,3 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo

            endif  !endif for r(0,k_1)

         endif  !endif for indt

         indpar=indpar+3
         indshell=indshell+3
         indorb=indorbp


        case(22)       ! 3p single zeta
c      3p without cusp condition
c      r e^{-z1 r } 


         k_1=kion(1)
         dd1=dd(indpar+1)
         if(iflagnorm.gt.2) then 
         c=dsqrt((2.d0*dd1)**7/240.d0/pi)/2.d0
         endif 
c
         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1))
         distp(k,2)=r(k,k_1)*distp(k,1)
         enddo
c
         indorbp=indorb
c         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,2)
               enddo
            endif
         enddo
c
c
         if(indt.ne.1) then
            fun0=distp(0,2)
            fun=(1.d0-dd1*r(0,k_1))*distp(0,1)
            fun2=dd1*(dd1*r(0,k_1)-2.d0)*distp(0,1)
c
           if(r(0,k_1).gt.1d-9) then
               indorbp=indorb
c
               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun/r(0,k_1)
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                 z(indorbp,indt+4)=rmu(ic,0,k_1)*
     &           (4.d0*fun/r(0,k_1)+fun2)
c
                 endif
               enddo
c
            else
c               
               indorbp=indorb 

               do ic=1,3 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=0.d0
                     enddo
                     z(indorbp,indt+4)=0.d0
                  endif
               enddo
c
            endif  !endif for r(0)
c
         endif  !endif for indt
c
         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp


        case(23)       ! 3p double zeta
c      3p without cusp condition
c      r ( e^{-z2 r } + z1 e^{-z3 r } )         


         k_1=kion(1)
         dd1=dd(indpar+1)
         dd2=dd(indpar+2)
         dd3=dd(indpar+3)
         if(iflagnorm.gt.2) then 
         c=1.d0/2.d0/dsqrt(240.d0*pi*(1.d0/(2.d0*dd1)**7
     &   +2.d0*dd3/(dd1+dd2)**7+dd3**2/(2.d0*dd2)**7))
         endif 
c
         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1))
         distp(k,2)=c*dexp(-dd2*r(k,k_1))
         enddo
c
         do i=indtmax,indtm
         distp(i,3)=r(i,k_1)*(distp(i,1)+dd3*distp(i,2))
         enddo 
c         
         indorbp=indorb
c         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,3)
               enddo
            endif
         enddo
c
c
         if(indt.ne.1) then
            fun0=distp(0,3)
            fun=(1.d0-dd1*r(0,k_1))*distp(0,1)
     1+dd3*(1.d0-dd2*r(0,k_1))*distp(0,2)
            fun2=dd1*(dd1*r(0,k_1)-2.d0)*distp(0,1)
     1+dd3*dd2*(dd2*r(0,k_1)-2.d0)*distp(0,2)
c
           if(r(0,k_1).gt.1d-9) then
               indorbp=indorb
c
               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun/r(0,k_1)
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                 z(indorbp,indt+4)=rmu(ic,0,k_1)*
     &           (4.d0*fun/r(0,k_1)+fun2)
c
                 endif
               enddo
c
            else
c               
               indorbp=indorb 

               do ic=1,3 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=0.d0
                     enddo
                     z(indorbp,indt+4)=0.d0
                  endif
               enddo
c
            endif  !endif for r(0)
c
         endif  !endif for indt
c
         indpar=indpar+3
         indshell=indshell+3
         indorb=indorbp



      case(24)    ! 4p single zeta
cc     4p without cusp condition
cc      r^2   e^{-z1 r }           

         k_1=kion(1)
         dd1=dd(indpar+1)
         if(iflagnorm.gt.2) then 
         c=dsqrt((2.d0*dd1)**9/120960.d0/pi)/2.d0
         endif 

         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1))
         enddo

         do i=indtmax,indtm
         distp(i,3)=r(i,k_1)**2*distp(i,1)
         enddo 
         
         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,3)
               enddo
            endif
         enddo


         if(indt.ne.1) then
            fun0=distp(0,3)
            fun=(2.d0*r(0,k_1)-dd1*r(0,k_1)**2)*distp(0,1)
            fun2=((dd1*r(0,k_1))**2+2.d0-4.d0*dd1*r(0,k_1))*distp(0,1)
            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun/r(0,k_1)
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                z(indorbp,indt+4)=rmu(ic,0,k_1)*(4.d0*fun/r(0,k_1)+fun2)
                  endif
               enddo
            else
               indorbp=indorb 
               do ic=1,3 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=0.d0
                    enddo
                     z(indorbp,indt+4)=0.d0
                  endif
               enddo

            endif  !endif for r(0)

         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp


      case(25)      ! 4p double zeta
c     4p without cusp condition
c      r^2  ( e^{-z2 r } + z1 e^{-z3 r } )         

         k_1=kion(1)

         dd1=dd(indpar+1)
         dd2=dd(indpar+2)
         dd3=dd(indpar+3)
         if(iflagnorm.gt.2) then 
        c=1.d0/2.d0/dsqrt(120960.d0*pi*(1.d0/(2.d0*dd1)**9
     &   +2.d0*dd3/(dd1+dd2)**9+dd3**2/(2.d0*dd2)**9))
         endif 
         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1))
         distp(k,2)=c*dexp(-dd2*r(k,k_1))
         enddo

         do i=indtmax,indtm
         distp(i,3)=r(i,k_1)**2*(distp(i,1)+dd3*distp(i,2))
         enddo 
         
         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,3)
               enddo
            endif
         enddo


         if(indt.ne.1) then
            fun0=distp(0,3)
            fun=(2.d0*r(0,k_1)-dd1*r(0,k_1)**2)*distp(0,1)
     1+dd3*(2.d0*r(0,k_1)-dd2*r(0,k_1)**2)*distp(0,2)
            fun2=((dd1*r(0,k_1))**2+2.d0-4.d0*dd1*r(0,k_1))*distp(0,1)
     1+dd3*((dd2*r(0,k_1))**2+2.d0-4.d0*dd2*r(0,k_1))*distp(0,2)
            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb
               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun/r(0,k_1)
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                z(indorbp,indt+4)=rmu(ic,0,k_1)*(4.d0*fun/r(0,k_1)+fun2)
                  endif
               enddo
            else
               indorbp=indorb 
               do ic=1,3 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=0.d0
                     enddo
                     z(indorbp,indt+4)=0.d0
                  endif
               enddo

            endif  !endif for r(0)

         endif  !endif for indt

         indpar=indpar+3
         indshell=indshell+3
         indorb=indorbp


      case(26)     ! 2p triple zeta
!     2p without cusp condition

         k_1=kion(1)

         dd1=dd(indpar+1)
         dd2=dd(indpar+2)
         peff=dd(indpar+3)
         dd3=dd(indpar+4)
         peff2=dd(indpar+5)

         if(iflagnorm.gt.2) then 
         c=1.d0/2.d0/dsqrt(8.d0*pi*(1.d0/(2.d0*dd1)**5
     &   +2.d0*peff/(dd1+dd2)**5+peff**2/(2.d0*dd2)**5
     &   +2.d0*peff2/(dd1+dd3)**5+peff2**2/(2.d0*dd3)**5
     &   +2.d0*peff2*peff/(dd2+dd3)**5))
         endif

         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1))
         distp(k,2)=c*dexp(-dd2*r(k,k_1))
         distp(k,3)=c*dexp(-dd3*r(k,k_1))
         enddo

         do i=indtmax,indtm
         distp(i,4)=distp(i,1)+peff*distp(i,2)+peff2*distp(i,3)
         enddo 
         
         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,4)
               enddo
            endif
         enddo


         if(indt.ne.1) then
       fun=(-dd1*distp(0,1)-dd2*peff*distp(0,2)
     & -dd3*peff2*distp(0,3))/r(0,k_1)
       fun2=dd1**2*distp(0,1)+peff*dd2**2*distp(0,2)
     & +peff2*dd3**2*distp(0,3)

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                  z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*fun
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+distp(0,4)
                     enddo
            z(indorbp,indt+4)=rmu(ic,0,k_1)*(4.d0*fun+fun2)
                  endif
               enddo

            else
               
               indorbp=indorb 

               do ic=1,3 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo

            endif  !endif for r(0,k_1)

         endif  !endif for indt

         indpar=indpar+5
         indshell=indshell+3
         indorb=indorbp

      case(27)     ! 3p triple zeta
!     2p without cusp condition

         k_1=kion(1)

         dd1=dd(indpar+1)
         dd2=dd(indpar+2)
         peff=dd(indpar+3)
         dd3=dd(indpar+4)
         peff2=dd(indpar+5)

         if(iflagnorm.gt.2) then 
         c=1.d0/2.d0/dsqrt(240.d0*pi*(1.d0/(2.d0*dd1)**7
     &   +2.d0*peff/(dd1+dd2)**7+peff**2/(2.d0*dd2)**7
     &   +2.d0*peff2/(dd1+dd3)**7+peff2**2/(2.d0*dd3)**7
     &   +2.d0*peff2*peff/(dd2+dd3)**7))
         endif

         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1))
         distp(k,2)=c*dexp(-dd2*r(k,k_1))
         distp(k,3)=c*dexp(-dd3*r(k,k_1))
         enddo

         do i=indtmax,indtm
         distp(i,4)=r(i,k_1)*(distp(i,1)+peff*distp(i,2)
     &   +peff2*distp(i,3))
         enddo 
         
         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,4)
               enddo
            endif
         enddo


         if(indt.ne.1) then
         fun0=distp(0,4)
         fun=(1.d0-dd1*r(0,k_1))*distp(0,1)
     1+peff*(1.d0-dd2*r(0,k_1))*distp(0,2)
     1+peff2*(1.d0-dd3*r(0,k_1))*distp(0,3)
         fun2=dd1*(dd1*r(0,k_1)-2.d0)*distp(0,1)
     1+peff*dd2*(dd2*r(0,k_1)-2.d0)*distp(0,2)
     1+peff2*dd3*(dd3*r(0,k_1)-2.d0)*distp(0,3)

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                  z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &            fun/r(0,k_1)
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)*
     &               (4.d0*fun/r(0,k_1)+fun2)
                  endif
               enddo

            else
               
               indorbp=indorb 

               do ic=1,3 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo

            endif  !endif for r(0,k_1)

         endif  !endif for indt

         indpar=indpar+5
         indshell=indshell+3
         indorb=indorbp


       case(30) 
!     3d without cusp and one parmater

          k_1=kion(1)
          dd1=dd(indpar+1)
          if(iflagnorm.gt.2) then 
      c=
     & 1.d0/(2.d0**3*3.d0)*dsqrt(1.d0/pi)*(2.d0*dd1)**(7.d0/2.d0) 
          endif 

          do k=indtmax,indtm
          distp(k,1)=c*dexp(-dd1*r(k,k_1))
          enddo

          do i=indtmax,indtm
             distp(i,3)=distp(i,1)
             distp(i,4)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d ! lz=0
             distp(i,5)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d ! lz=+/-2
             distp(i,6)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d  ! lz=+/-2
             distp(i,7)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
             distp(i,8)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
          enddo 
         
          indorbp=indorb

          do ic=1,5
             if(iocc(indshell+ic).eq.1) then
                indorbp=indorbp+1
                do i=1,indtm
                   z(indorbp,i)=distp(i,3+ic)*distp(i,3)
                enddo
             endif
          enddo

          if(indt.ne.1) then
             fun0=distp(0,3)
             fun=-dd1*distp(0,1)
             fun2=dd1**2*distp(0,1)

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                     z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0,k_1)
     &                *fun/r(0,k_1)
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d                      
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
           z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0,k_1)+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

           else
c               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                      z(indorbp,indt+i)=0.d0
                   enddo
                 endif
              enddo
               
           endif  !endif for r(0)

        endif   !endif for indt
c
        indpar=indpar+1
        indshell=indshell+5
        indorb=indorbp

              case(31) 
!     3d without cusp condition double Z
 
         k_1=kion(1)          
         dd1=dd(indpar+1)
         dd2=dd(indpar+2)
         peff=dd(indpar+3)
         if(iflagnorm.gt.2) then 
      c=1/2.d0*dsqrt(5.d0/pi)
     $/dsqrt(1/dd1**7/128.d0+2*peff/(dd1+dd2)**7
     $+peff**2/dd2**7/128.d0)/dsqrt(720.d0)
         endif 

         do k=indtmax,indtm
            distp(k,1)=dexp(-dd1*r(k,k_1))
            distp(k,2)=dexp(-dd2*r(k,k_1))
         enddo

         do i=indtmax,indtm
            distp(i,3)=c*(distp(i,1)+peff*distp(i,2))
            distp(i,4)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d !lz=0 
            distp(i,5)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d !lz=+/-2
            distp(i,6)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d  ! lz=+/- 2
            distp(i,7)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
            distp(i,8)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
         enddo 
         
         indorbp=indorb

         do ic=1,5
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
                  z(indorbp,i)=distp(i,3+ic)*distp(i,3)
               enddo
            endif
         enddo

         if(indt.ne.1) then
            fun0=distp(0,3)
            fun=c*(-dd1*distp(0,1)-peff*dd2*distp(0,2))
            fun2=c*(dd1**2*distp(0,1)
     1           +peff*dd2**2*distp(0,2))

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
            z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0,k_1)*fun/r(0,k_1)
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d         
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d          
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
               z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0,k_1)+fun2)
                 endif  !endif for iocc
              enddo     ! enddo fot ic

           else
               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                       z(indorbp,indt+i)=0.d0
                    enddo
                 endif
              enddo
               
           endif  !endif for r(0)

        endif   !endif for indt

        indpar=indpar+3
        indshell=indshell+5
        indorb=indorbp


              case(32) 
!     3d without cusp condition triple Z
 
         k_1=kion(1)          
         dd1=dd(indpar+1)
         dd2=dd(indpar+2)
         peff=dd(indpar+3)
         dd3=dd(indpar+4)
         peff2=dd(indpar+5)

         if(iflagnorm.gt.2) then 
      c=1/2.d0*dsqrt(5.d0/pi)
     $/dsqrt(1/(2.d0*dd1)**7+2*peff/(dd1+dd2)**7
     $+peff**2/(2.d0*dd2)**7+2*peff2/(dd1+dd3)**7
     $+peff2**2/(2.d0*dd3)**7+2*peff*peff2/(dd2+dd3)**7)/dsqrt(720.d0)
         endif 

         do k=indtmax,indtm
            distp(k,1)=dexp(-dd1*r(k,k_1))
            distp(k,2)=dexp(-dd2*r(k,k_1))
            distp(k,3)=dexp(-dd3*r(k,k_1))
         enddo

         do i=indtmax,indtm
            distp(i,4)=c*(distp(i,1)+peff*distp(i,2)+peff2*distp(i,3))
            distp(i,5)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d !lz=0 
            distp(i,6)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d !lz=+/-2
            distp(i,7)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d  ! lz=+/- 2
            distp(i,8)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
            distp(i,9)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
         enddo 
         
         indorbp=indorb

         do ic=1,5
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
                  z(indorbp,i)=distp(i,4+ic)*distp(i,4)
               enddo
            endif
         enddo

         if(indt.ne.1) then
            fun0=distp(0,4)
            fun=c*(-dd1*distp(0,1)-peff*dd2*distp(0,2)
     &      -peff2*dd3*distp(0,3))
            fun2=c*(dd1**2*distp(0,1)+peff*dd2**2*distp(0,2)
     &      +peff2*dd3**2*distp(0,3))

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
            z(indorbp,indt+i)=distp(0,4+ic)*rmu(i,0,k_1)*fun/r(0,k_1)
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d         
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d          
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
               z(indorbp,indt+4)=distp(0,4+ic)*(6.d0*fun/r(0,k_1)+fun2)
                 endif  !endif for iocc
              enddo     ! enddo fot ic

           else
               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                       z(indorbp,indt+i)=0.d0
                    enddo
                 endif
              enddo
               
           endif  !endif for r(0)

        endif   !endif for indt

        indpar=indpar+5
        indshell=indshell+5
        indorb=indorbp

       case(33) 
!     4d without cusp and one parmater

          k_1=kion(1)
          dd1=dd(indpar+1)
          if(iflagnorm.gt.2) then 
!     c=
!    & 1/(2.d0**3*3.d0)*dsqrt(1.d0/pi)*(2.d0*dd1)**(7.d0/2.d0) 
      c=
     & 1.d0/(2.d0**3*3.d0)/dsqrt(56.d0*pi)*(2.d0*dd1)**(9.d0/2.d0) 
          endif 

          do k=indtmax,indtm
          distp(k,1)=c*dexp(-dd1*r(k,k_1))
          enddo

          do i=indtmax,indtm
             distp(i,3)=distp(i,1)*r(i,k_1)
             distp(i,4)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d ! lz=0
             distp(i,5)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d ! lz=+/-2
             distp(i,6)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d  ! lz=+/-2
             distp(i,7)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
             distp(i,8)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
          enddo 
         
          indorbp=indorb

          do ic=1,5
             if(iocc(indshell+ic).eq.1) then
                indorbp=indorbp+1
                do i=1,indtm
                   z(indorbp,i)=distp(i,3+ic)*distp(i,3)
                enddo
             endif
          enddo

          if(indt.ne.1) then
             fun0=distp(0,3)
             fun=-dd1*distp(0,3)+distp(0,1)
             fun2=dd1**2*distp(0,3)-2.d0*dd1*distp(0,1)

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                     z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0,k_1)
     &                *fun/r(0,k_1)
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d                      
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
           z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0,k_1)+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

           else
c               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                      z(indorbp,indt+i)=0.d0
                   enddo
                 endif
              enddo
               
           endif  !endif for r(0)

        endif   !endif for indt
c
        indpar=indpar+1
        indshell=indshell+5
        indorb=indorbp


          case(34)     ! 2s single  Z WITH CUSP zero
                                    ! normalized
! exp(-dd1*r) + dd1*r*exp(-dd1*r)

         k_1=kion(1)
         indshellp=indshell+1


         if(iocc(indshellp).eq.1) then 

            indorbp=indorb+1
            dd1=dd(indpar+1)
!           peff=dd1
            

            if(iflagnorm.gt.2) then 
            c=1.d0/dsqrt(1.d0/4.d0/dd1**3+12.d0*dd1/(2.d0*dd1)**4+
     &3.d0*dd1**2/4.d0/dd1**5)/dsqrt(4.d0*pi)
            endif 

            do k=indtmax,indtm
            distp(k,1)=c*dexp(-dd1*r(k,k_1))
            enddo
             
            do i=1,indtm
            z(indorbp,i)=distp(i,1)*(1.d0+r(i,k_1)*dd1)
            enddo            

            if(indt.ne.1) then
      
             fun=-dd1**2*distp(0,1)
             fun2=fun*(1.d0-dd1*r(0,k_1))

                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo
                  z(indorbp,indt+4)=(2.d0*fun+fun2)

            endif

            indorb=indorbp

         endif
         indpar=indpar+1
         indshell=indshellp

          case(35)     ! 2s single  Z WITH CUSP
                                    ! normalized
! exp(-dd1*r) + dd1* r * exp(-dd2*r)

         k_1=kion(1)
         indshellp=indshell+1


         if(iocc(indshellp).eq.1) then 

            indorbp=indorb+1
            dd1=dd(indpar+1)
            dd2=dd(indpar+2)
            peff=dd1
            
            do k=indtmax,indtm
               distp(k,1)=dexp(-dd1*r(k,k_1))
               distp(k,2)=dexp(-dd2*r(k,k_1))
            enddo

            if(iflagnorm.gt.2) then 
            c=1.d0/dsqrt(1/4.d0/dd1**3+12*peff/(dd1+dd2)**4+
     &3*peff**2/4/dd2**5)/dsqrt(4.0*pi)
            endif 
             
            do i=1,indtm
            z(indorbp,i)=c*(distp(i,1)+r(i,k_1)*distp(i,2)*peff)
            enddo            

            if(indt.ne.1) then
      
             fun=-dd1*distp(0,1)+peff*distp(0,2)*(1.d0-dd2*r(0,k_1))
             fun2=distp(0,1)*dd1**2
     &+peff*distp(0,2)*(dd2**2*r(0,k_1)-2.d0*dd2)

               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
                  z(indorbp,indt+i)=fun*c*rmu(i,0,k_1)/r(0,k_1)
                  enddo
                  z(indorbp,indt+4)=c*(2.d0*fun/r(0,k_1)+fun2)
               else

                  do i=1,4
                     z(indorbp,indt+i)=0.d0 
                  enddo
               endif

            endif

            indorb=indorbp

         endif
         indpar=indpar+2
         indshell=indshellp

      case(36)  ! single gaussian p orbitals

         k_1=kion(1)

         dd(indpar+1)=abs(dd(indpar+1))
         dd1=dd(indpar+1)


         if(iflagnorm.gt.2) then 
         c=dsqrt(2.d0)*pi**(-0.75d0)*(2.d0*dd1)**1.25d0
         endif 


         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1)**2)
         enddo
         
         indorbp=indorb
c         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
                  z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)
              enddo
            endif
         enddo

         if(indt.ne.1) then
            fun0=distp(0,1)
            fun=-2.d0*dd1*distp(0,1)
            fun2=fun*(1.d0-2.d0*dd1*r(0,k_1)**2)

               indorbp=indorb

               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
             z(indorbp,indt+4)=rmu(ic,0,k_1)*(4.d0*fun+fun2)
                  endif
               enddo


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

       case(37)
! d orbitals
! R(r)= exp(-alpha r^2)
! each gaussian term is normalized


         k_1=kion(1)
         indorbp=indorb
         indparp=indpar+1

         dd(indpar+1)=abs(dd(indpar+1))
         dd1=dd(indparp)

         if(iflagnorm.gt.2) then 
! overall normalization 
         c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)
         endif


         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1)**2)
         enddo


          do i=indtmax,indtm
      distp(i,2)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d ! lz=0
      distp(i,3)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d ! lz=+/-2
      distp(i,4)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d  ! lz=+/-2
      distp(i,5)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
      distp(i,6)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
          enddo 


         do ic=1,5
            if(iocc(indshell+ic).eq.1) then
              indorbp=indorbp+1
              do k=1,indtm
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)
              enddo
            endif
         enddo


         if(indt.ne.1) then
            
            dd1=dd(indparp)
            fun0=distp(0,1)
            fun=-2.d0*dd1*distp(0,1)
            fun2=fun*(1.d0-2.d0*dd1*r(0,k_1)**2)


               indorbp=indorb
               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0,k_1)
     &                *fun
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d     
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d         
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
      z(indorbp,indt+4)=distp(0,1+ic)*(6.d0*fun+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+5
         indorb=indorbp         

      case(47) 
! d orbitals cartesian !!!
! R(r)= exp(-alpha r^2)
! each gaussian term is normalized


         k_1=kion(1)
         indorbp=indorb
         indparp=indpar+1
         dd1=dd(indparp)

         if(iflagnorm.gt.2) then 
! overall normalization 
         c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)
         endif


         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1)**2)
         enddo


          do i=indtmax,indtm
      distp(i,2)=rmu(1,i,k_1)**2
      distp(i,3)=rmu(2,i,k_1)**2
      distp(i,4)=rmu(3,i,k_1)**2
      distp(i,5)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d  ! lz=+/-2
      distp(i,6)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
      distp(i,7)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
          enddo 


         do ic=1,6
            if(iocc(indshell+ic).eq.1) then
              indorbp=indorbp+1
              do k=1,indtm
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)
              enddo
            endif
         enddo


         if(indt.ne.1) then
            
            dd1=dd(indparp)
            fun0=distp(0,1)
            fun=-2.d0*dd1*distp(0,1)
            fun2=fun*(1.d0-2.d0*dd1*r(0,k_1)**2)

               indorbp=indorb
               do ic=1,6
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0,k_1)
     &                *fun
                       if(ic.le.3) then
                          if(i.eq.ic) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.6) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
                    if(ic.le.3) then
      z(indorbp,indt+4)=distp(0,1+ic)*(6.d0*fun+fun2)+2.d0*distp(0,1)
                    else
      z(indorbp,indt+4)=distp(0,1+ic)*(6.d0*fun+fun2)
                    endif
                 endif  !endif for iocc
             enddo     ! enddo fot ic


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+6
         indorb=indorbp   


      case(48)  
! f single gaussian orbital 
! R(r)= exp(-alpha r^2)
! normalized

         k_1=kion(1)
         indorbp=indorb
         indparp=indpar+1

         dd(indpar+1)=abs(dd(indpar+1))
         dd1=dd(indparp)

         if(iflagnorm.gt.2) then 
! overall normalization 
         c=8.d0/dsqrt(15.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(9.d0/4.d0)
         endif


         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1)**2)
         enddo


          do i=indtmax,indtm
      distp(i,2)=cost1f*rmu(3,i,k_1)
     &           *(5.d0*rmu(3,i,k_1)**2-3.d0*r(i,k_1)**2) ! lz=0
      distp(i,3)=cost2f*rmu(1,i,k_1)
     &           *(5.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)  ! lz=+/-1
      distp(i,4)=cost2f*rmu(2,i,k_1)
     &           *(5.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)  ! lz=+/-1
      distp(i,5)=cost3f*rmu(3,i,k_1)
     &           *(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2) ! lz=+/-2
      distp(i,6)=cost3f*2.d0*rmu(3,i,k_1)
     &           *rmu(1,i,k_1)*rmu(2,i,k_1) ! lz=+/-2 
      distp(i,7)=cost4f*rmu(1,i,k_1)
     &           *(rmu(1,i,k_1)**2-3.d0*rmu(2,i,k_1)**2)  ! lz=+/-3
      distp(i,8)=cost4f*rmu(2,i,k_1)
     &           *(3.d0*rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)  ! lz=+/-3
          enddo 


         do ic=1,7
            if(iocc(indshell+ic).eq.1) then
              indorbp=indorbp+1
              do k=1,indtm
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)
              enddo
            endif
         enddo


         if(indt.ne.1) then
            
            dd1=dd(indparp)
            fun0=distp(0,1)
            fun=-2.d0*dd1*distp(0,1)
            fun2=fun*(1.d0-2.d0*dd1*r(0,k_1)**2)


               indorbp=indorb
               do ic=1,7
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0,k_1)
     &                *fun
                       if(ic.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                 6.d0*cost1f*fun0*rmu(i,0,k_1)*rmu(3,0,k_1)
                          if(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &       cost1f*fun0*(15.d0*rmu(i,0,k_1)**2-3.d0*r(0,k_1)**2)
                          endif
                       elseif(ic.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                 2.d0*cost2f*fun0*rmu(i,0,k_1)*rmu(1,0,k_1)
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &             cost2f*fun0*(5.d0*rmu(3,0,k_1)**2-r(0,k_1)**2)
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                10.d0*cost2f*fun0*rmu(i,0,k_1)*rmu(1,0,k_1)
                          endif
                       elseif(ic.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                 2.d0*cost2f*fun0*rmu(i,0,k_1)*rmu(2,0,k_1)
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &             cost2f*fun0*(5.d0*rmu(3,0,k_1)**2-r(0,k_1)**2)
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                10.d0*cost2f*fun0*rmu(i,0,k_1)*rmu(2,0,k_1)
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                 2.d0*cost3f*fun0*rmu(1,0,k_1)*rmu(3,0,k_1)
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                 2.d0*cost3f*fun0*rmu(2,0,k_1)*rmu(3,0,k_1)
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &              cost3f*fun0*(rmu(1,0,k_1)**2-rmu(2,0,k_1)**2)
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                 2.d0*cost3f*fun0*rmu(2,0,k_1)*rmu(3,0,k_1)
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                 2.d0*cost3f*fun0*rmu(1,0,k_1)*rmu(3,0,k_1)
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                 2.d0*cost3f*fun0*rmu(1,0,k_1)*rmu(2,0,k_1)
                          endif
                       elseif(ic.eq.6) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &         3.d0*cost4f*fun0*(rmu(1,0,k_1)**2-rmu(2,0,k_1)**2)
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                 6.d0*cost4f*fun0*rmu(1,0,k_1)*rmu(2,0,k_1)
                          endif
                       else
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                 6.d0*cost4f*fun0*rmu(1,0,k_1)*rmu(2,0,k_1)
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &         3.d0*cost4f*fun0*(rmu(1,0,k_1)**2-rmu(2,0,k_1)**2)
                          endif
                       endif    !endif for ic 
                    enddo  !enddo for i
      z(indorbp,indt+4)=distp(0,1+ic)*(8.d0*fun+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+7
         indorb=indorbp         


       case(49) ! derivative of 48 with respect to z
! f orbitals
! R(r)= c*exp(-z r^2)*(9/4/z-r^2)


         k_1=kion(1)
         indorbp=indorb
         indparp=indpar+1
         dd(indpar+1)=abs(dd(indpar+1))
         dd1=dd(indparp)

         if(iflagnorm.gt.2) then 
! overall normalization 
         c=8.d0/dsqrt(15.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(9.d0/4.d0)
         endif


         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1)**2)
         enddo


          do i=indtmax,indtm
      distp(i,2)=cost1f*rmu(3,i,k_1)
     &           *(5.d0*rmu(3,i,k_1)**2-3.d0*r(i,k_1)**2) ! lz=0
      distp(i,3)=cost2f*rmu(1,i,k_1)
     &           *(5.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)  ! lz=+/-1
      distp(i,4)=cost2f*rmu(2,i,k_1)
     &           *(5.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)  ! lz=+/-1
      distp(i,5)=cost3f*rmu(3,i,k_1)
     &           *(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2) ! lz=+/-2
      distp(i,6)=cost3f*2.d0*rmu(3,i,k_1)
     &           *rmu(1,i,k_1)*rmu(2,i,k_1) ! lz=+/-2 
      distp(i,7)=cost4f*rmu(1,i,k_1)
     &           *(rmu(1,i,k_1)**2-3.d0*rmu(2,i,k_1)**2)  ! lz=+/-3
      distp(i,8)=cost4f*rmu(2,i,k_1)
     &           *(3.d0*rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)  ! lz=+/-3
          enddo 


         do ic=1,7
            if(iocc(indshell+ic).eq.1) then
              indorbp=indorbp+1
              do k=1,indtm
              z(indorbp,k)=distp(k,1)*(9.d0/4.d0/dd1-r(k,k_1)**2)*
     &        distp(k,1+ic)
              enddo
            endif
         enddo


         if(indt.ne.1) then
            
            dd1=dd(indparp)
            fun0=distp(0,1)*(9.d0/4.d0/dd1-r(0,k_1)**2)
            fun=distp(0,1)*(2.d0*dd1*r(0,k_1)**2-13.d0/2.d0)
            fun2=distp(0,1)*(-4.d0*dd1**2*r(0,k_1)**4
     &      +19.d0*dd1*r(0,k_1)**2-13.d0/2.d0)


               indorbp=indorb
               do ic=1,7
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0,k_1)
     &                *fun
                       if(ic.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                 6.d0*cost1f*fun0*rmu(i,0,k_1)*rmu(3,0,k_1)
                          if(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &       cost1f*fun0*(15.d0*rmu(i,0,k_1)**2-3.d0*r(0,k_1)**2)
                          endif
                       elseif(ic.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                 2.d0*cost2f*fun0*rmu(i,0,k_1)*rmu(1,0,k_1)
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &             cost2f*fun0*(5.d0*rmu(3,0,k_1)**2-r(0,k_1)**2)
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                10.d0*cost2f*fun0*rmu(i,0,k_1)*rmu(1,0,k_1)
                          endif
                       elseif(ic.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                 2.d0*cost2f*fun0*rmu(i,0,k_1)*rmu(2,0,k_1)
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &             cost2f*fun0*(5.d0*rmu(3,0,k_1)**2-r(0,k_1)**2)
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                10.d0*cost2f*fun0*rmu(i,0,k_1)*rmu(2,0,k_1)
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                 2.d0*cost3f*fun0*rmu(1,0,k_1)*rmu(3,0,k_1)
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                 2.d0*cost3f*fun0*rmu(2,0,k_1)*rmu(3,0,k_1)
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &              cost3f*fun0*(rmu(1,0,k_1)**2-rmu(2,0,k_1)**2)
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                 2.d0*cost3f*fun0*rmu(2,0,k_1)*rmu(3,0,k_1)
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                 2.d0*cost3f*fun0*rmu(1,0,k_1)*rmu(3,0,k_1)
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                 2.d0*cost3f*fun0*rmu(1,0,k_1)*rmu(2,0,k_1)
                          endif
                       elseif(ic.eq.6) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &         3.d0*cost4f*fun0*(rmu(1,0,k_1)**2-rmu(2,0,k_1)**2)
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                 6.d0*cost4f*fun0*rmu(1,0,k_1)*rmu(2,0,k_1)
                          endif
                       else
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                 6.d0*cost4f*fun0*rmu(1,0,k_1)*rmu(2,0,k_1)
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &         3.d0*cost4f*fun0*(rmu(1,0,k_1)**2-rmu(2,0,k_1)**2)
                          endif
                       endif    !endif for ic 
                    enddo  !enddo for i
      z(indorbp,indt+4)=distp(0,1+ic)*(8.d0*fun+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+7
         indorb=indorbp   



      case(38)      ! 3s -derivative of 34 with respect to dd1 
c R(r)=r**2*exp(-z1*r)

      
         k_1=kion(1)   
         indshellp=indshell+1

         if(iocc(indshellp).eq.1) then 
                 
            indorbp=indorb+1
            dd1=dd(indpar+1)
            if(iflagnorm.gt.2) then 
!           c=dsqrt((2*dd1)**7/720.d0/pi)/2.d0
            c=1.d0/dsqrt(1.d0/4.d0/dd1**3+12.d0*dd1/(2.d0*dd1)**4+
     &3.d0*dd1**2/4.d0/dd1**5)/dsqrt(4.d0*pi)


            endif 
           
            c0=-c*dd1

            c1=1.5d0*c/dd1


            
            do k=indtmax,indtm
            distp(k,1)=dexp(-dd1*r(k,k_1))
            enddo

            do i=1,indtm
            z(indorbp,i)=(c0*r(i,k_1)**2+c1*(1.d0+dd1*r(i,k_1)))
     1*distp(i,1)
            enddo            

            c1=c1*dd1**2

            if(indt.ne.1) then
!            fun=-dd1**2*distp(0,1)
!            fun2=fun*(1.d0-dd1*r(0,k_1))
               fun=(c0*(2.d0-dd1*r(0,k_1))-c1)*distp(0,1)
               fun2=(c0*(2.d0-4*dd1*r(0,k_1)+(dd1*r(0,k_1))**2) 
     1+c1*(dd1*r(0,k_1)-1.d0))*distp(0,1)
                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo
                  z(indorbp,indt+4)=2.d0*fun+fun2
            endif

            indorb=indorbp

         endif
         indpar=indpar+1
         indshell=indshellp         

      case(39)      ! 4s single zeta derivative of 10
c R(r)=r**3*exp(-z1*r)
c         
         indshellp=indshell+1
         k_1=kion(1)


        if(iocc(indshellp).eq.1) then 

            indorbp=indorb+1
            dd1=dd(indpar+1)
            if(iflagnorm.gt.2) then 
!           c=dsqrt((2*dd1)**9/40320.d0/pi)/2.d0
            c=dsqrt((2*dd1)**7/720.d0/pi)/2.d0
!           c=-c
            endif 

            c0=-c
            c1=3.5d0*c/dd1

            do k=indtmax,indtm
            distp(k,1)=dexp(-dd1*r(k,k_1))
            enddo

            do i=1,indtm
            z(indorbp,i)=(c0*r(i,k_1)**3+c1*r(i,k_1)**2)*distp(i,1)
            enddo            

            if(indt.ne.1) then
               rp1=r(0,k_1)**3
               rp2=r(0,k_1)**2

!              fun=(2.d0-dd1*r(0,k_1))*distp(0,1)
!              fun2=(2.d0-4*dd1*r(0,k_1)+(dd1*r(0,k_1))**2)*distp(0,1)
c
cc              the first derivative/r 
               fun=distp(0,1)*(c0*(3.d0*r(0,k_1)-dd1*rp2)
     1+c1*(2.d0-dd1*r(0,k_1)))

cc

cc              the second derivative
               fun2=distp(0,1)*
     1(c0*(6.d0*r(0,k_1)-6.d0*dd1*rp2+dd1**2*rp1)
     1+c1*(2.d0-4*dd1*r(0,k_1)+(dd1*r(0,k_1))**2))
cc
                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

                  z(indorbp,indt+4)=2.d0*fun+fun2

            endif
c
            indorb=indorbp
c
         endif
         indpar=indpar+1
         indshell=indshellp         
c
        case(40)       ! 3p single zeta
c      3p without cusp condition derivative of 20
c      r e^{-z1 r } 


         k_1=kion(1)
         dd1=dd(indpar+1)
         if(iflagnorm.gt.2) then 
         c=dsqrt((2.d0*dd1)**5/8.d0/pi)/2.d0

         endif 

         c0=-c
         c1=2.5d0*c/dd1

c
         do k=indtmax,indtm
         distp(k,1)=dexp(-dd1*r(k,k_1))
         distp(k,2)=r(k,k_1)*distp(k,1)
         enddo
c
         indorbp=indorb
c         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*(c0*distp(i,2)+c1*distp(i,1))
               enddo
            endif
         enddo
c
c
         if(indt.ne.1) then
            fun0=c0*distp(0,2)+c1*distp(0,1)
            fun=(c0*(1.d0-dd1*r(0,k_1))-c1*dd1)*distp(0,1)
            fun2=(c0*dd1*(dd1*r(0,k_1)-2.d0)+c1*dd1**2)*distp(0,1)
c
           if(r(0,k_1).gt.1d-9) then
               indorbp=indorb
c
               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun/r(0,k_1)
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                 z(indorbp,indt+4)=rmu(ic,0,k_1)*
     &           (4.d0*fun/r(0,k_1)+fun2)
c
                 endif
               enddo
c
            else
c               
               indorbp=indorb 

               do ic=1,3 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=0.d0
                     enddo
                     z(indorbp,indt+4)=0.d0
                  endif
               enddo
c
            endif  !endif for r(0)
c
         endif  !endif for indt
c
         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp


      case(41)    ! 4p single zeta
cc     4p without cusp condition derivative of 22 
cc      r^2   e^{-z1 r }           

         k_1=kion(1)
         dd1=dd(indpar+1)
         if(iflagnorm.gt.2) then 
         c=dsqrt((2.d0*dd1)**7/240.d0/pi)/2.d0
!        c=dsqrt((2.d0*dd1)**9/120960.d0/pi)/2.d0
         endif 

         c0=-c

         c1=3.5d0*c/dd1

         do k=indtmax,indtm
         distp(k,1)=dexp(-dd1*r(k,k_1))
         enddo

         do i=indtmax,indtm
         distp(i,3)=r(i,k_1)**2*distp(i,1)
         enddo 
         
         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
      z(indorbp,i)=rmu(ic,i,k_1)*(c0*distp(i,3)+c1*r(i,k_1)*distp(i,1))
               enddo
            endif
         enddo


         if(indt.ne.1) then
!           fun=(1.d0-dd1*r(0,k_1))*distp(0,1)
!           fun2=dd1*(dd1*r(0,k_1)-2.d0)*distp(0,1)

            fun0=c0*distp(0,3)+c1*r(0,k_1)*distp(0,1)
            fun=(c0*(2.d0-dd1*r(0,k_1))*r(0,k_1)
     1+c1*(1.d0-dd1*r(0,k_1)))*distp(0,1)
            fun2=(c0*((dd1*r(0,k_1))**2+2.d0-4.d0*dd1*r(0,k_1))
     1+c1*dd1*(dd1*r(0,k_1)-2.d0))*distp(0,1)

           if(r(0,k_1).gt.1d-9) then
               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun/r(0,k_1)
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
          z(indorbp,indt+4)=rmu(ic,0,k_1)*(4.d0*fun/r(0,k_1)+fun2)
                  endif
               enddo
           else

               indorbp=indorb 

               do ic=1,3 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=0.d0
                     enddo
                     z(indorbp,indt+4)=0.d0
                  endif
               enddo


           endif 





         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

       case(42) 
!     4d without cusp and one parmater derivative of 30

          k_1=kion(1)
          dd1=dd(indpar+1)
          if(iflagnorm.gt.2) then 
      c=
     & 1.d0/(2.d0**3*3.d0)*dsqrt(1.d0/pi)*(2.d0*dd1)**(7.d0/2.d0) 
!     c=
!    & 1.d0/(2.d0**3*3.d0)/dsqrt(56.d0*pi)*(2.d0*dd1)**(9.d0/2.d0) 
          endif 

          c0=-c
          c1=3.5d0*c/dd1

          do k=indtmax,indtm
          distp(k,1)=dexp(-dd1*r(k,k_1))
          enddo

          do i=indtmax,indtm
             distp(i,3)=distp(i,1)*(c0*r(i,k_1)+c1)
             distp(i,4)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d ! lz=0
             distp(i,5)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d ! lz=+/-2
             distp(i,6)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d  ! lz=+/-2
             distp(i,7)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
             distp(i,8)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
          enddo 
         
          indorbp=indorb

          do ic=1,5
             if(iocc(indshell+ic).eq.1) then
                indorbp=indorbp+1
                do i=1,indtm
                   z(indorbp,i)=distp(i,3+ic)*distp(i,3)
                enddo
             endif
          enddo

          if(indt.ne.1) then
             fun0=distp(0,3)
             fun=-dd1*distp(0,3)+c0*distp(0,1)
             fun2=dd1**2*distp(0,3)-2.d0*dd1*c0*distp(0,1)

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                     z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0,k_1)
     &                *fun/r(0,k_1)
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d                      
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
           z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0,k_1)+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

           else
c               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                      z(indorbp,indt+i)=0.d0
                   enddo
                 endif
              enddo
               
           endif  !endif for r(0)

        endif   !endif for indt
c
        indpar=indpar+1
        indshell=indshell+5
        indorb=indorbp

       case(43) 
!     4d without cusp and one parmater  derivative of 33 

          k_1=kion(1)
          dd1=dd(indpar+1)
          if(iflagnorm.gt.2) then 
      c=
     & 1.d0/(2.d0**3*3.d0)/dsqrt(56.d0*pi)*(2.d0*dd1)**(9.d0/2.d0) 
          endif 

          c0=-c
          c1=4.5d0*c/dd1


          do k=indtmax,indtm
          distp(k,1)=dexp(-dd1*r(k,k_1))
          enddo

          do i=indtmax,indtm
             distp(i,3)=distp(i,1)*(c0*r(i,k_1)**2+c1*r(i,k_1))
             distp(i,4)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d ! lz=0
             distp(i,5)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d ! lz=+/-2
             distp(i,6)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d  ! lz=+/-2
             distp(i,7)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
             distp(i,8)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
          enddo 
         
          indorbp=indorb

          do ic=1,5
             if(iocc(indshell+ic).eq.1) then
                indorbp=indorbp+1
                do i=1,indtm
                   z(indorbp,i)=distp(i,3+ic)*distp(i,3)
                enddo
             endif
          enddo

          if(indt.ne.1) then
             fun0=distp(0,3)
             fun=-dd1*distp(0,3)+distp(0,1)*(2.d0*c0*r(0,k_1)+c1) 
             fun2=dd1**2*distp(0,3)+distp(0,1)*
     1(-2.d0*dd1*(2.d0*c0*r(0,k_1)+c1)+2.d0*c0)
            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                     z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0,k_1)
     &                *fun/r(0,k_1)
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d                      
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
           z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0,k_1)+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

           else
c               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                      z(indorbp,indt+i)=0.d0
                   enddo
                 endif
              enddo
               
           endif  !endif for r(0)

        endif   !endif for indt
c
        indpar=indpar+1
        indshell=indshell+5
        indorb=indorbp


       case(44) ! derivative of 36 with respect zeta
! R(r)=x*exp(-z*r**2)*(5/4/z-r**2)

         k_1=kion(1)
         dd(indpar+1)=abs(dd(indpar+1))
         dd1=dd(indpar+1)
         if(iflagnorm.gt.2) then 
         c=dsqrt(2.d0)*pi**(-0.75d0)*(2.d0*dd1)**1.25d0
         endif 


         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1)**2)
         enddo
         
         indorbp=indorb
c         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
                   z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)*
     &             (5.d0/4.d0/dd1-r(i,k_1)**2)
               enddo
            endif
         enddo

         if(indt.ne.1) then
            fun0=distp(0,1)*(5.d0/4.d0/dd1-r(0,k_1)**2)
            fun=distp(0,1)*(2.d0*dd1*r(0,k_1)**2-9.d0/2.d0)
            fun2=distp(0,1)*(-4.d0*dd1**2*r(0,k_1)**4
     &      +15.d0*dd1*r(0,k_1)**2-9.d0/2.d0)

               indorbp=indorb

               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
             z(indorbp,indt+4)=rmu(ic,0,k_1)*(4.d0*fun+fun2)
                  endif
               enddo


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp


       case(45) ! derivative of 37 with respect to z
! d orbitals
! R(r)= c*exp(-z r^2)*(7/4/z-r^2)


         k_1=kion(1)
         indorbp=indorb
         indparp=indpar+1
         dd(indpar+1)=abs(dd(indpar+1))
         dd1=dd(indparp)

         if(iflagnorm.gt.2) then 
! overall normalization 
         c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)
         endif


         do k=indtmax,indtm
         distp(k,1)=c*dexp(-dd1*r(k,k_1)**2)
         enddo


          do i=indtmax,indtm
      distp(i,2)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d ! lz=0
      distp(i,3)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d ! lz=+/-2
      distp(i,4)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d  ! lz=+/-2
      distp(i,5)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
      distp(i,6)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
          enddo 


         do ic=1,5
            if(iocc(indshell+ic).eq.1) then
              indorbp=indorbp+1
              do k=1,indtm
              z(indorbp,k)=distp(k,1)*(7.d0/4.d0/dd1-r(k,k_1)**2)*
     &        distp(k,1+ic)
              enddo
            endif
         enddo


         if(indt.ne.1) then
            
            dd1=dd(indparp)
            fun0=distp(0,1)*(7.d0/4.d0/dd1-r(0,k_1)**2)
            fun=distp(0,1)*(2.d0*dd1*r(0,k_1)**2-11.d0/2.d0)
            fun2=distp(0,1)*(-4.d0*dd1**2*r(0,k_1)**4
     &      +17.d0*dd1*r(0,k_1)**2-11.d0/2.d0)


               indorbp=indorb
               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0,k_1)
     &                *fun
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d                      
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
      z(indorbp,indt+4)=distp(0,1+ic)*(6.d0*fun+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+5
         indorb=indorbp   

      case(46)  ! derivative of 17 with respect to z
! R(r)=c*r**2*exp(-z*r**2)*(7/4/z-r**2)

         k_1=kion(1)   
         indshellp=indshell+1

         if(iocc(indshellp).eq.1) then 
                 
            indorbp=indorb+1
            dd1=dd(indpar+1)

            if(iflagnorm.gt.2) then 
         c=4.d0*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)/dsqrt(15.d0)
            endif 
            
            do k=indtmax,indtm
            distp(k,1)=c*dexp(-dd1*r(k,k_1)**2)
            enddo

            do i=1,indtm
            z(indorbp,i)=distp(i,1)*(7.d0/4.d0/dd1*r(i,k_1)**2
     &                   -r(i,k_1)**4)
            enddo            

            if(indt.ne.1) then
               rp1=r(0,k_1)**2
!              the first derivative / r 
            fun=distp(0,1)*(7.d0-15.d0*dd1*rp1
     &                   +4.d0*(dd1*rp1)**2)/2.d0/dd1
!              the second derivative
            fun2=distp(0,1)*(7.d0-59*dd1*rp1+50*(dd1*rp1)**2
     &                   -8*(dd1*rp1)**3)/2.d0/dd1 
                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo
                  z(indorbp,indt+4)=2.d0*fun+fun2
            endif

            indorb=indorbp

         endif

         indpar=indpar+1
         indshell=indshellp



      

! **************** GAUSSIAN BASIS *****************************

       case(300:399)
! s orbitals 
! R(r)= \sum_i a_i exp(-alpha_i r^2)
! each gaussian term is normalized

         ngauss=iopt-300

         k_1=kion(1)   
         indshellp=indshell+1

         if(iocc(indshellp).eq.1) then 
                 
            indorbp=indorb+1


         if(iflagnorm.gt.2) then 
! overall normalization 
         c=(2.d0/pi)**(3.d0/4.d0) 
         peff=0.d0
         do i=1,ngauss
            do k=1,i-1
               peff=peff+dd(indpar+ngauss+i)*dd(indpar+ngauss+k)
     &         *(dd(indpar+i)*dd(indpar+k))**(3.d0/4.d0)
     &         /(dd(indpar+i)+dd(indpar+k))**(3.d0/2.d0)
            enddo
         enddo
         peff=peff*4.d0*sqrt(2.d0)
         do i=1,ngauss
            peff=peff+dd(indpar+ngauss+i)**2
         enddo
!         c=c/sqrt(peff)
         endif

            do k=indtmax,indtm
            do i=1,ngauss
            dd1=dd(indpar+i)
            distp(k,i)=dd1**(3.d0/4.d0)*dexp(-dd1*r(k,k_1)**2)
            enddo
            enddo

            do k=1,indtm
            z(indorbp,k)=0.d0
            do i=1,ngauss
            z(indorbp,k)=z(indorbp,k)+distp(k,i)*dd(indpar+ngauss+i)
            enddo
            z(indorbp,k)=c*z(indorbp,k)
            enddo            

            if(indt.ne.1) then
               rp1=r(0,k_1)**2

!              the first and second derivative 
            fun=0.d0
            fun2=0.d0
            do i=1,ngauss
            dd1=dd(indpar+i)
            fun=fun-2.d0*dd1*dd(indpar+ngauss+i)*distp(0,i)*r(0,k_1)
            fun2=fun2-2.d0*dd1*dd(indpar+ngauss+i)*
     &           distp(0,i)*(1.d0-2.d0*dd1*rp1)
            enddo
            fun=c*fun
            fun2=c*fun2

               if(r(0,k_1).gt.1d-9) then
                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)/r(0,k_1)
                  enddo

                  z(indorbp,indt+4)=2.d0*fun/r(0,k_1)+fun2
               else
                  do i=1,4
                     z(indorbp,indt+i)=0.d0 
                  enddo
               endif

            endif

            indorb=indorbp

         endif
         indpar=indpar+2*ngauss
         indshell=indshellp



       case(400:499)
! p orbitals
! R(r)= \sum_i a_i x exp(-alpha_i r^2)
! each gaussian term is normalized

         ngauss=iopt-400

         k_1=kion(1)
         indorbp=indorb

         if(iflagnorm.gt.2) then 
! overall normalization 
         c=2.d0*(2.d0/pi)**(3.d0/4.d0)
         peff=0.d0
         do i=1,ngauss
            do k=1,i-1
               peff=peff+dd(indpar+ngauss+i)*dd(indpar+ngauss+k)
     &         *(dd(indpar+i)*dd(indpar+k))**(5.d0/4.d0)
     &         /(dd(indpar+i)+dd(indpar+k))**(5.d0/2.d0)
            enddo
         enddo
         peff=peff*8.d0*sqrt(2.d0)
         do i=1,ngauss
            peff=peff+dd(indpar+ngauss+i)**2
         enddo
!         c=c/sqrt(peff)
         endif


         do k=indtmax,indtm
         do i=1,ngauss
         dd1=dd(indpar+i)
         distp(k,i)=dd1**(5.d0/4.d0)*dexp(-dd1*r(k,k_1)**2)
         enddo
         enddo

         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do k=1,indtm
                  z(indorbp,k)=0.d0
                  do i=1,ngauss
                     z(indorbp,k)=z(indorbp,k)+
     &               distp(k,i)*dd(indpar+ngauss+i)
                  enddo
                  z(indorbp,k)=c*rmu(ic,k,k_1)*z(indorbp,k)
              enddo
            endif
         enddo

         if(indt.ne.1) then
            rp1=r(0,k_1)**2
            fun0=0.d0
            fun=0.d0
            fun2=0.d0
            do i=1,ngauss
            dd1=dd(indpar+i)
            fun0=fun0+dd(indpar+ngauss+i)*distp(0,i)
            fun=fun-2.d0*dd1*dd(indpar+ngauss+i)*distp(0,i)*r(0,k_1)
            fun2=fun2-2.d0*dd1*dd(indpar+ngauss+i)*
     &           distp(0,i)*(1.d0-2.d0*dd1*rp1)
            enddo
            fun0=c*fun0
            fun=c*fun
            fun2=c*fun2

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                     do i=1,3
                        z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun/r(0,k_1)
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
             z(indorbp,indt+4)=rmu(ic,0,k_1)*(4.d0*fun/r(0,k_1)+fun2)
                  endif
               enddo

            else
               
               indorbp=indorb 
c
               do ic=1,3 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,4
                        z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo

            endif  !endif for r(0)

         endif  !endif for indt

         indpar=indpar+2*ngauss
         indshell=indshell+3
         indorb=indorbp         

  
       case(500:599)
! d orbitals
! R(r)= \sum_i a_i xy exp(-alpha_i r^2)
! each gaussian term is normalized

         ngauss=iopt-500

         k_1=kion(1)
         indorbp=indorb

         if(iflagnorm.gt.2) then 
! overall normalization 
         c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)
         peff=0.d0
         do i=1,ngauss
            do k=1,i-1
               peff=peff+dd(indpar+ngauss+i)*dd(indpar+ngauss+k)
     &         *(dd(indpar+i)*dd(indpar+k))**(7.d0/4.d0)
     &         /(dd(indpar+i)+dd(indpar+k))**(7.d0/2.d0)
            enddo
         enddo
         peff=peff*16.d0*sqrt(2.d0)
         do i=1,ngauss
            peff=peff+dd(indpar+ngauss+i)**2
         enddo
!         c=c/sqrt(peff)
         endif


         do k=indtmax,indtm
         do i=1,ngauss
         dd1=dd(indpar+i)
         distp(k,i)=dd1**(7.d0/4.d0)*dexp(-dd1*r(k,k_1)**2)
         enddo
         enddo


          do i=indtmax,indtm
      distp(i,ngauss+1)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d ! lz=0
      distp(i,ngauss+2)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d ! lz=+/-2
      distp(i,ngauss+3)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d  ! lz=+/-2
      distp(i,ngauss+4)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
      distp(i,ngauss+5)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d  ! lz=+/-1
          enddo 


         do ic=1,5
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do k=1,indtm
                  z(indorbp,k)=0.d0
                  do i=1,ngauss
                     z(indorbp,k)=z(indorbp,k)+
     &               distp(k,i)*dd(indpar+ngauss+i)
                  enddo
                  z(indorbp,k)=c*distp(k,ngauss+ic)*z(indorbp,k)
              enddo
            endif
         enddo


         if(indt.ne.1) then
            rp1=r(0,k_1)**2
            fun0=0.d0
            fun=0.d0
            fun2=0.d0
            do i=1,ngauss
            dd1=dd(indpar+i)
            fun0=fun0+dd(indpar+ngauss+i)*distp(0,i)
            fun=fun-2.d0*dd1*dd(indpar+ngauss+i)*distp(0,i)*r(0,k_1)
            fun2=fun2-2.d0*dd1*dd(indpar+ngauss+i)*
     &           distp(0,i)*(1.d0-2.d0*dd1*rp1)
            enddo
            fun0=c*fun0
            fun=c*fun
            fun2=c*fun2


            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                z(indorbp,indt+i)=distp(0,ngauss+ic)*rmu(i,0,k_1)
     &                *fun/r(0,k_1)
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d                      
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
      z(indorbp,indt+4)=distp(0,ngauss+ic)*(6.d0*fun/r(0,k_1)+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

           else
               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                      z(indorbp,indt+i)=0.d0
                   enddo
                 endif
              enddo
               
           endif  !endif for r(0)


         endif  !endif for indt

         indpar=indpar+2*ngauss
         indshell=indshell+5
         indorb=indorbp         





! ******************* END GAUSSIAN BASIS ************************

! ** ** ** ** ** ** **  JASTROW ORBITALS ** ** ** ** ** ** ** ** *
       case(100)  
c     2s single gaussian  
c     exp(-dd2*r^2)
         k_1=kion(1)

         dd2=dd(indpar+1)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=dexp(-dd2*r(k,k_1)*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)
               enddo
            endif


         if(indt.ne.1) then
               fun=-dd2*distp(0,1)*2.d0

                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*dd2*(-3.d0+2.d0*dd2*r(0,k_1)**2)*
     1     distp(0,1)


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshellp
         indorb=indorbp



       case(101)  
c     2s without cusp condition
c     dd1*( dd3 +exp(-dd2*r^2))
         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=dexp(-dd2*r(k,k_1)*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)+dd3
               enddo
            endif


         if(indt.ne.1) then
               fun=-dd2*distp(0,1)*2.d0

                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*dd2*(-3.d0+2.d0*dd2*r(0,k_1)**2)*
     1     distp(0,1)


         endif  !endif for indt

         indpar=indpar+2
         indshell=indshellp
         indorb=indorbp


            case(102)  
c     2s double gaussian with constant
c     (dd3+ exp (-dd2 r^2)+dd4*exp(-dd5*r^2)) 

         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)
         dd4=dd(indpar+3)
         dd5=dd(indpar+4)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=dexp(-dd2*r(k,k_1)*r(k,k_1))
               distp(k,2)=dexp(-dd5*r(k,k_1)*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=(distp(i,1)+dd3+dd4*distp(i,2))
c              write(6,*) ' function inside = ',z(indorbp,i)
               enddo
            endif


         if(indt.ne.1) then
               fun=-2.d0*(dd2*distp(0,1)+dd5*dd4*distp(0,2))
               fun2=r(0,k_1)**2

c              write(6,*) ' fun inside = ',fun,fun2 

                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*(dd2*(-3.d0+2.d0*dd2*fun2)*
     1     distp(0,1)+dd5*dd4*(-3.d0+2.d0*dd5*fun2)*distp(0,2))

c           write(6,*) ' lap 106 =',z(indorbp,indt+4)


c           stop 

         endif  !endif for indt

         indpar=indpar+4
         indshell=indshellp
         indorb=indorbp



       case(104)  
c     2p  double gaussian 
c       dd1 * x_mu  (exp(-dd2 r^2)+dd3 * exp(-dd4*r^2))
         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)
         dd4=dd(indpar+3)


         do k=indtmax,indtm
         distp(k,1)=dexp(-dd2*r(k,k_1)**2)
         distp(k,2)=dexp(-dd4*r(k,k_1)**2)
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
        z(indorbp,i)=rmu(ic,i,k_1)*(distp(i,1)+dd3*distp(i,2))
               enddo
            endif
         enddo


         if(indt.ne.1) then
            fun0=(distp(0,1)+dd3*distp(0,2))
            fun=2.d0*(-dd2*distp(0,1)
     1      -dd4*dd3*distp(0,2))
      fun2=2.d0*(dd2*(-1.d0+2.d0*dd2*r(0,k_1)**2)*distp(0,1)    
     1+dd4*dd3*(-1.d0+2.d0*dd4*r(0,k_1)**2)*distp(0,2)) 


               indorbp=indorb

               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo

         endif  !endif for indt

         indpar=indpar+3
         indshell=indshell+3
         indorb=indorbp
          
       case(103)  
c     2p single gaussian  

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=dexp(-dd2*r(k,k_1)**2)
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)
               enddo
            endif
         enddo


         if(indt.ne.1) then
            fun0=distp(0,1)
            fun=-dd2*distp(0,1)*2.d0
            fun2=2.d0*dd2*(-1.d0+2.d0*dd2*r(0,k_1)**2)*
     1     distp(0,1)    

               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo

         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

            case(105)  
c     2s double gaussian without  constant
c     (exp (-dd2 r^2)+dd4*exp(-dd5*r^2)) 

         k_1=kion(1)

c        dd1=1.d0 
         dd2=dd(indpar+1)
c        dd3=dd(indpar+2)
c        dd4=dd(indpar+3)
c        dd5=dd(indpar+4)
         dd4=dd(indpar+2)
         dd5=dd(indpar+3)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=dexp(-dd2*r(k,k_1)*r(k,k_1))
               distp(k,2)=dexp(-dd5*r(k,k_1)*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)+dd4*distp(i,2)
c              write(6,*) ' function inside = ',z(indorbp,i)
               enddo
            endif


         if(indt.ne.1) then
               fun=-2.d0*(dd2*distp(0,1)+dd5*dd4*distp(0,2))
               fun2=r(0,k_1)**2

c              write(6,*) ' fun inside = ',fun,fun2 

                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*(dd2*(-3.d0+2.d0*dd2*fun2)*
     1     distp(0,1)+dd5*dd4*(-3.d0+2.d0*dd5*fun2)*distp(0,2))




         endif  !endif for indt

         indpar=indpar+3
         indshell=indshellp
         indorb=indorbp

       case(106)  
c     2s without cusp condition
c     dd1*( dd3 +1/(1+dd2*r^2))
         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=1.d0/(1.d0+dd2*r(k,k_1)*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)+dd3
               enddo
            endif


         if(indt.ne.1) then
               fun=-dd2*distp(0,1)**2*2.d0
               fun2=fun*distp(0,1)*(1.-3.d0*dd2*r(0,k_1)**2)

                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=fun2+2.d0*fun


         endif  !endif for indt

         indpar=indpar+2
         indshell=indshellp
         indorb=indorbp

       case(107)  
c     2p single  lorentian  parent of 103 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=1.d0/(1.d0+dd2*r(k,k_1)**2)
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)
               enddo
            endif
         enddo


         if(indt.ne.1) then
            fun0=distp(0,1)
            fun=-dd2*distp(0,1)**2*2.d0
            fun2=fun*distp(0,1)*(1.d0-3.d0*dd2*r(0,k_1)**2)

               indorbp=indorb

               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

            case(108)  
c     2s double lorentian with constant  parent of 102 
c     (dd3+ L(dd2 r^2)+dd4*L(dd5*r^2)) ;  L(x)=1/1+x^2

         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)
         dd4=dd(indpar+3)
         dd5=dd(indpar+4)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=1.d0/(1.d0+dd2*r(k,k_1)*r(k,k_1))
               distp(k,2)=1.d0/(1.d0+dd5*r(k,k_1)*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=(distp(i,1)+dd3+dd4*distp(i,2))
c              write(6,*) ' function inside = ',z(indorbp,i)
               enddo
            endif


         if(indt.ne.1) then
          fun=-2.d0*(dd2*distp(0,1)**2+dd5*dd4*distp(0,2)**2)
          fun2=2.d0*dd2*distp(0,1)**3*(-1.d0+3.d0*dd2*r(0,k_1)**2)
     1    +2.d0*dd5*dd4*distp(0,2)**3*(-1.d0+3.d0*dd5*r(0,k_1)**2)

c              write(6,*) ' fun inside = ',fun,fun2 

                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=fun2+2.d0*fun

c           write(6,*) ' lap 106 =',z(indorbp,indt+4)

         endif  !endif for indt

         indpar=indpar+4
         indshell=indshellp
         indorb=indorbp


       case(109)  
c     2p  double  Lorentian  
c       dd1 * x_mu  (L(dd2 r^2)+dd3 * L(dd4*r^2)) ; L(x)=1/(1+x^2)

         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)
         dd4=dd(indpar+3)

         do k=indtmax,indtm
         distp(k,1)=1.d0/(1.d0+dd2*r(k,k_1)**2)
         distp(k,2)=1.d0/(1.d0+dd4*r(k,k_1)**2)
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
        z(indorbp,i)=rmu(ic,i,k_1)*(distp(i,1)+dd3*distp(i,2))
               enddo
            endif
         enddo


         if(indt.ne.1) then
            fun0=distp(0,1)+dd3*distp(0,2)

            fun=2.d0*(-dd2*distp(0,1)**2-dd4*dd3*distp(0,2)**2)
c     fun2=2.d0*(dd2*(-1.d0+2.d0*dd2*r(0,k_1)**2)*distp(0,1)    
c    1+dd4*dd3*(-1.d0+2.d0*dd4*r(0,k_1)**2)*distp(0,2)) 

        fun2=2*dd2*distp(0,1)**3*(-1.d0+3.d0*dd2*r(0,k_1)**2)
     1+2*dd3*dd4*distp(0,2)**3*(-1.d0+3.d0*dd4*r(0,k_1)**2)

               indorbp=indorb

               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo

         endif  !endif for indt

         indpar=indpar+3
         indshell=indshell+3
         indorb=indorbp
          
       case(110)  
c     2s without cusp condition
c     dd1*( dd3 +1/(1+dd2*r^3))
         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=1.d0/(1.d0+dd2*r(k,k_1)**3)
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)+dd3
               enddo
            endif


         if(indt.ne.1) then
               fun=-dd2*distp(0,1)**2*3.d0*r(0,k_1)
               fun2=fun*distp(0,1)*(2.d0-4.d0*dd2*r(0,k_1)**3)

                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=fun2+2.d0*fun


         endif  !endif for indt

         indpar=indpar+2
         indshell=indshellp
         indorb=indorbp

       case(111)  
c     2p single   r_mu/(1+b r^3)   parent of 103 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=1.d0/(1.d0+dd2*r(k,k_1)**3)
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)
               enddo
            endif
         enddo


         if(indt.ne.1) then
            fun0=distp(0,1)
            fun=-dd2*distp(0,1)**2*3.d0*r(0,k_1)
            fun2=fun*distp(0,1)*(2.d0-4.d0*dd2*r(0,k_1)**3)

               indorbp=indorb

               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp
       case(112)  
c     2p single   r_mu/(1+b r)^3   parent of 103 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=1.d0/(1.d0+dd2*r(k,k_1))**3
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)
               enddo
            endif
         enddo



         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then 
            fun0=distp(0,1)
            fun=-3.d0*dd2*distp(0,1)/(r(0,k_1)*(1.d0+dd2*r(0,k_1)))
            fun2=12.d0*dd2**2/(1.+dd2*r(0,k_1))**5

               indorbp=indorb

               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo
            else

               indorbp=indorb

               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo

            endif  !endif for r(0,k_1)

         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

       case(113)  
c     2s without cusp condition
c     dd1*( dd3 +r^2/(1+dd2*r)^4)
         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=r(k,k_1)**2/(1.d0+dd2*r(k,k_1))**4
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)+dd3
               enddo
            endif


         if(indt.ne.1) then
               fun= (2.d0-2.d0*dd2*r(0,k_1))/(1+dd2*r(0,k_1))**5
          fun2=2.d0*(1.d0-6.d0*dd2*r(0,k_1)+3.d0*(dd2*r(0,k_1))**2)
     1/(1+dd2*r(0,k_1))**6
                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=fun2+2.d0*fun


         endif  !endif for indt

         indpar=indpar+2
         indshell=indshellp
         indorb=indorbp

       case(114)  
c     2s without cusp condition
c     dd1*( dd3 +r^2/(1+dd2*r)^3)
         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=r(k,k_1)**2/(1.d0+dd2*r(k,k_1))**3
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)+dd3
               enddo
            endif


         if(indt.ne.1) then
               fun= (2.d0-dd2*r(0,k_1))/(1+dd2*r(0,k_1))**4
               fun2=2.d0*(1.d0-4.d0*dd2*r(0,k_1)+(dd2*r(0,k_1))**2)
     1/(1+dd2*r(0,k_1))**5
                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=fun2+2.d0*fun


         endif  !endif for indt

         indpar=indpar+2
         indshell=indshellp
         indorb=indorbp

         case(115)  
c     2s double lorentian with constant  parent of 102 
c     (dd3+ r^2/(1+dd2*r)^3+dd4*r^3/(1+dd5*r)^4; 

         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)
         dd4=dd(indpar+3)
         dd5=dd(indpar+4)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=r(k,k_1)**2/(1.d0+dd2*r(k,k_1))**3
               distp(k,2)=r(k,k_1)**3/(1.d0+dd5*r(k,k_1))**4
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=(distp(i,1)+dd3+dd4*distp(i,2))
c              write(6,*) ' function inside = ',z(indorbp,i)
               enddo
            endif


         if(indt.ne.1) then

               fun= (2.d0-dd2*r(0,k_1))/(1+dd2*r(0,k_1))**4
     1   -dd4*r(0,k_1)*(-3.d0+dd5*r(0,k_1))/(1.d0+dd5*r(0,k_1))**5
               fun2=2.d0*(1.d0-4.d0*dd2*r(0,k_1)+(dd2*r(0,k_1))**2)
     1/(1+dd2*r(0,k_1))**5
     1+dd4*2.d0*r(0,k_1)*(3.d0-6.d0*dd5*r(0,k_1)+(dd5*r(0,k_1))**2)
     1/(1.d0+dd5*r(0,k_1))**6


c              write(6,*) ' fun inside = ',fun,fun2 

                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=fun2+2.d0*fun

c           write(6,*) ' lap 106 =',z(indorbp,indt+4)

         endif  !endif for indt

         indpar=indpar+4
         indshell=indshellp
         indorb=indorbp

       case(116)  
c     2p  double  Lorentian  
c       dd1 * x_mu  (L^3(dd2 r)+dd3 r * L(dd4*r)^4) ; L(x)=1/(1+x)

         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)
         dd4=dd(indpar+3)

         do k=indtmax,indtm
         distp(k,1)=1.d0/(1.d0+dd2*r(k,k_1))**3
         distp(k,2)=r(k,k_1)/(1.d0+dd4*r(k,k_1))**4
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
        z(indorbp,i)=rmu(ic,i,k_1)*(distp(i,1)+dd3*distp(i,2))
               enddo
            endif
         enddo


         if(indt.ne.1) then

            if(r(0,k_1).gt.1.d-9) then  

            fun0=distp(0,1)+dd3*distp(0,2)
       fun=-3.d0*dd2*distp(0,1)/(r(0,k_1)*(1.d0+dd2*r(0,k_1)))
     1+dd3*distp(0,2)/r(0,k_1)**2*(1.d0-3*dd4*r(0,k_1))
     1/(1.d0+dd4*r(0,k_1))
            fun2=12.d0*dd2**2/(1.+dd2*r(0,k_1))**5
     1+dd3*4.d0*dd4*(-2.d0+3.d0*dd4*r(0,k_1))/(1.+dd4*r(0,k_1))**6

c           fun0=distp(0,1)+dd3*distp(0,2)
c           fun=2.d0*(-dd2*distp(0,1)**2-dd4*dd3*distp(0,2)**2)

c       fun2=2*dd2*distp(0,1)**3*(-1.d0+3.d0*dd2*r(0,k_1)**2)
c    1+2*dd3*dd4*distp(0,2)**3*(-1.d0+3.d0*dd4*r(0,k_1)**2)

               indorbp=indorb

               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo

               else

               indorbp=indorb

               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo

            endif  !endif for r(0,k_1)


         endif  !endif for indt

         indpar=indpar+3
         indshell=indshell+3
         indorb=indorbp

         case(117)  
c     2s double lorentian with constant  parent of 102 
c     (dd3+r^3/(1+dd5*r)^4; 

         k_1=kion(1)

         dd3=dd(indpar+1)
         dd5=dd(indpar+2)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=r(k,k_1)**3/(1.d0+dd5*r(k,k_1))**4
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=dd3+distp(i,1)
c              write(6,*) ' function inside = ',z(indorbp,i)
               enddo
            endif


         if(indt.ne.1) then

               fun= 
     1   -r(0,k_1)*(-3.d0+dd5*r(0,k_1))/(1.d0+dd5*r(0,k_1))**5
               fun2=
     1+2.d0*r(0,k_1)*(3.d0-6.d0*dd5*r(0,k_1)+(dd5*r(0,k_1))**2)
     1/(1.d0+dd5*r(0,k_1))**6


c              write(6,*) ' fun inside = ',fun,fun2 

                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=fun2+2.d0*fun

c           write(6,*) ' lap 106 =',z(indorbp,indt+4)

         endif  !endif for indt

         indpar=indpar+2
         indshell=indshellp
         indorb=indorbp

         case(118)  
c     2s double lorentian with constant  parent of 102 
c     (dd1+  1/ (1 + Exp[  dd2 (r^2 - r_0^2) ] )   | dd3=r_0
c      Fermi distribution with r^2 
         k_1=kion(1)

         dd1=dd(indpar+1)
         dd2=dd(indpar+2)
         dd3=-dd2*dd(indpar+3)**2

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               arg=dd2*r(k,k_1)**2+dd3
               if(arg.gt.200) then 
               distp(k,1)=dexp(200.d0)
               else
               distp(k,1)=dexp(arg)
               endif 
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=dd1+1.d0/(1.d0+distp(i,1))
c              write(6,*) ' function inside = ',z(indorbp,i)
               enddo
            endif


         if(indt.ne.1) then

               fun= -2.d0*dd2*distp(0,1)/(1.d0+distp(0,1))**2
         fun2=-2.d0*dd2*(-distp(0,1)*(-1.d0-2.d0*dd2*r(0,k_1)**2)
     1+distp(0,1)**2*(1.d0-2.d0*dd2*r(0,k_1)**2))/(1.d0+distp(0,1))**3


c              write(6,*) ' fun inside = ',fun,fun2 

                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=fun2+2.d0*fun

c           write(6,*) ' lap 106 =',z(indorbp,indt+4)

         endif  !endif for indt

         indpar=indpar+3
         indshell=indshellp
         indorb=indorbp

       case(119)  
c     2p single   r_mu/(1+b r^2)^(3/2)   parent of 103 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=1.d0/(1.d0+dd2*r(k,k_1)**2)**1.5d0
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)
               enddo
            endif
         enddo



         if(indt.ne.1) then

            fun0=distp(0,1)
            fun=-3.d0*dd2*distp(0,1)/(1.d0+dd2*r(0,k_1)**2)
            fun2=3.d0*dd2*(-1.d0+4.d0*dd2*r(0,k_1)**2)
     1/(1.d0+dd2*r(0,k_1)**2)**3.5d0

               indorbp=indorb

               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo

         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

       case(120)  
c     2p  double  cubic   
c       dd1 * x_mu  (L^3(dd2 r)+dd3 L(dd4*r)^3) ; L(x)=1/(1+x)

         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)
         dd4=dd(indpar+3)

         do k=indtmax,indtm
         distp(k,1)=1.d0/(1.d0+dd2*r(k,k_1))**3
         distp(k,2)=1.d0/(1.d0+dd4*r(k,k_1))**3
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
        z(indorbp,i)=rmu(ic,i,k_1)*(distp(i,1)+dd3*distp(i,2))
               enddo
            endif
         enddo


         if(indt.ne.1) then

            if(r(0,k_1).gt.1.d-9) then  

            fun0=distp(0,1)+dd3*distp(0,2)
       fun=-3.d0*dd2*distp(0,1)/(r(0,k_1)*(1.d0+dd2*r(0,k_1)))
     1 -3.d0*dd4*dd3*distp(0,2)/(r(0,k_1)*(1.d0+dd4*r(0,k_1)))
            fun2=12.d0*dd2**2/(1.+dd2*r(0,k_1))**5
     1      +12.d0*dd3*dd4**2/(1.+dd4*r(0,k_1))**5

c           fun0=distp(0,1)+dd3*distp(0,2)
c           fun=2.d0*(-dd2*distp(0,1)**2-dd4*dd3*distp(0,2)**2)

c       fun2=2*dd2*distp(0,1)**3*(-1.d0+3.d0*dd2*r(0,k_1)**2)
c    1+2*dd3*dd4*distp(0,2)**3*(-1.d0+3.d0*dd4*r(0,k_1)**2)

               indorbp=indorb

               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo

               else

               indorbp=indorb

               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo

            endif  !endif for r(0,k_1)


         endif  !endif for indt

         indpar=indpar+3
         indshell=indshell+3
         indorb=indorbp

       case(121)  
c     2p single exponential   

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=dexp(-dd2*r(k,k_1))
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)
               enddo
            endif
         enddo


         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then 

            fun0=distp(0,1)
            fun=-dd2*distp(0,1)/r(0,k_1)
            fun2=dd2**2*distp(0,1)    

               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo

           else
               indorbp=indorb

               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo

           endif 


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

       case(122)  
c     2s with cusp condition
c     dd1*( dd3 +exp(-dd2*r)*(1+dd2*r))
         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=dexp(-dd2*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)*(1.d0+dd2*r(i,k_1))+dd3
               enddo
            endif


         if(indt.ne.1) then
               fun=-dd2**2*distp(0,1)
               fun2=fun*(1.d0-dd2*r(0,k_1))


                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*fun+fun2


         endif  !endif for indt

         indpar=indpar+2
         indshell=indshellp
         indorb=indorbp

       case(123)  
c     2p  double  exp  
c       dd1 * x_mu  (exp(-dd2 r)+dd3 * exp(-dd4*r))
         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)
         dd4=dd(indpar+3)


         do k=indtmax,indtm
         distp(k,1)=dexp(-dd2*r(k,k_1))
         distp(k,2)=dexp(-dd4*r(k,k_1))
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
        z(indorbp,i)=rmu(ic,i,k_1)*(distp(i,1)+dd3*distp(i,2))
               enddo
            endif
         enddo


         if(indt.ne.1) then
            if(r(0,k_1).gt.1d-9) then 

            fun0=distp(0,1)+dd3*distp(0,2)
            fun=-(dd2*distp(0,1)+dd3*dd4*distp(0,2))/r(0,k_1)
            fun2=dd2**2*distp(0,1)+dd3*dd4**2*distp(0,2)    


               indorbp=indorb

               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo

         else
               indorbp=indorb

               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo

         endif 


         endif  !endif for indt

         indpar=indpar+3
         indshell=indshell+3
         indorb=indorbp

            case(124)  
c     2s double exp  with constant and cusp cond. 
c     (dd3+ exp (-dd2 r)*(1+dd2*r)+dd4*exp(-dd5*r)*(1+dd5*r)) 

         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)
         dd4=dd(indpar+3)
         dd5=dd(indpar+4)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,3)=dexp(-dd2*r(k,k_1))
               distp(k,4)=dexp(-dd5*r(k,k_1))
               distp(k,1)=distp(k,3)*(1.d0+dd2*r(k,k_1))
               distp(k,2)=distp(k,4)*(1.d0+dd5*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)+dd3+dd4*distp(i,2)
c              write(6,*) ' function inside = ',z(indorbp,i)
               enddo
            endif


         if(indt.ne.1) then

               fun=-dd2**2*distp(0,3)-dd5**2*dd4*distp(0,4)
               fun2=-dd2**2*distp(0,3)*(1.d0-dd2*r(0,k_1))
     1        -dd4*dd5**2*distp(0,4)*(1.d0-dd5*r(0,k_1))



                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

         z(indorbp,indt+4)=2.d0*fun+fun2


         endif  !endif for indt

         indpar=indpar+4
         indshell=indshellp
         indorb=indorbp

       case(125)  
c     2s with cusp condition
c     dd1*( dd3 +exp(-dd2*r))  ! with no cusp condition 
         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=dexp(-dd2*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)+dd3
               enddo
            endif


         if(indt.ne.1) then
              if(r(0,k_1).gt.1d-9) then 
               fun=-dd2*distp(0,1)/r(0,k_1)
               fun2=dd2**2*distp(0,1)

                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*fun+fun2

               else

               do i=1,4
               z(indorbp,indt+i)=0.d0
               enddo
               endif 




         endif  !endif for indt

         indpar=indpar+2
         indshell=indshellp
         indorb=indorbp

            case(126)  
c     2s double exp  with constant
c     (dd3+ exp (-dd2 r)+dd4*exp(-dd5*r)) 

         k_1=kion(1)

         dd2=dd(indpar+1)
         dd3=dd(indpar+2)
         dd4=dd(indpar+3)
         dd5=dd(indpar+4)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=dexp(-dd2*r(k,k_1))
               distp(k,2)=dexp(-dd5*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)+dd3+dd4*distp(i,2)
c              write(6,*) ' function inside = ',z(indorbp,i)
               enddo
            endif


         if(indt.ne.1) then

             if(r(0,k_1).gt.1d-9) then 

               fun=-(dd2*distp(0,1)+dd5*dd4*distp(0,2))/r(0,k_1)
               fun2=dd2**2*distp(0,1)+dd4*dd5**2*distp(0,2)



                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

         z(indorbp,indt+4)=2.d0*fun+fun2


                else

                do i=1,4
                z(indorbp,indt+i)=0.d0
                enddo

                endif

         endif  !endif for indt

         indpar=indpar+4
         indshell=indshellp
         indorb=indorbp

       case(127) 
!     3d without cusp and one parmater

          k_1=kion(1)
          dd1=dd(indpar+1)

          do k=indtmax,indtm
          distp(k,1)=dexp(-dd1*r(k,k_1))
          enddo

          do i=indtmax,indtm
             distp(i,3)=distp(i,1)
             distp(i,4)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d
             distp(i,5)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d
             distp(i,6)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d
             distp(i,7)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d
             distp(i,8)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d
          enddo 
         
          indorbp=indorb

          do ic=1,5
             if(iocc(indshell+ic).eq.1) then
                indorbp=indorbp+1
                do i=1,indtm
                   z(indorbp,i)=distp(i,3+ic)*distp(i,3)
                enddo
             endif
          enddo

          if(indt.ne.1) then
             fun0=distp(0,3)
             fun=-dd1*distp(0,1)
             fun2=dd1**2*distp(0,1)

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                     z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0,k_1)
     &                *fun/r(0,k_1)
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d                      
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
           z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0,k_1)+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

           else
c               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                      z(indorbp,indt+i)=0.d0
                   enddo
                 endif
              enddo
               
           endif  !endif for r(0)

        endif   !endif for indt
c
        indpar=indpar+1
        indshell=indshell+5
        indorb=indorbp

       case(128)  
c     2s with cusp condition
c     ( r^2*exp(-dd2*r))  ! with no cusp condition 
         k_1=kion(1)

         dd2=dd(indpar+1)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=dexp(-dd2*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)*r(i,k_1)**2
               enddo
            endif


         if(indt.ne.1) then
               fun=(2.d0-dd2*r(0,k_1))*distp(0,1)
               fun2=(2.d0-4*dd2*r(0,k_1)+(dd2*r(0,k_1))**2)*distp(0,1)
                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*fun+fun2


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshellp
         indorb=indorbp

       case(129)  
c     2p single exponential  r e^{-z r}  ! parent of 121 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=dexp(-dd2*r(k,k_1))
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)*r(i,k_1)
               enddo
            endif
         enddo


         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then 
            fun0=distp(0,1)*r(0,k_1)
            fun=distp(0,1)*(1.d0-dd2*r(0,k_1))/r(0,k_1)
            fun2=dd2*(dd2*r(0,k_1)-2.d0)*distp(0,1)
               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo
           else
               indorbp=indorb
               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo
           endif 
         endif  !endif for indt
         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

       case(130)  
c     2p single exponential  r^2  e^{-z r}  ! parent of 121 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=dexp(-dd2*r(k,k_1))
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)*r(i,k_1)**2
               enddo
            endif
         enddo


         if(indt.ne.1) then

            fun0=distp(0,1)*r(0,k_1)**2
            fun=distp(0,1)*(2.d0-dd2*r(0,k_1))
            fun2=(2.d0-4.d0*dd2*r(0,k_1)+(dd2*r(0,k_1))**2)*distp(0,1)
               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo
         endif  !endif for indt
         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

       case(131)  
!     2s without cusp condition
!     dd1*(r^2*exp(-dd2*r^2))
         k_1=kion(1)

         dd2=dd(indpar+1)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=dexp(-dd2*r(k,k_1)*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)*r(i,k_1)**2
               enddo
            endif


         if(indt.ne.1) then
               fun0=dd2*r(0,k_1)**2
               fun=2.d0*distp(0,1)*(1.d0-fun0)
               fun2=2.d0*distp(0,1)*(1.d0-5.d0*fun0+2.d0*fun0**2)
  
                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*fun+fun2


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshellp
         indorb=indorbp

       case(132)  
c     2s with cusp condition
c     ( r^3*exp(-dd2*r))  ! with no cusp condition 
         k_1=kion(1)

         dd2=dd(indpar+1)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=dexp(-dd2*r(k,k_1))*r(k,k_1)
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)*r(i,k_1)**2
               enddo
            endif


         if(indt.ne.1) then
               fun=(3.d0-dd2*r(0,k_1))*distp(0,1)
               fun2=(6.d0-6*dd2*r(0,k_1)+(dd2*r(0,k_1))**2)*distp(0,1)
                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*fun+fun2


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshellp
         indorb=indorbp


       case(133) 
!     4d  one parmater

          k_1=kion(1)
          dd1=dd(indpar+1)

          do k=indtmax,indtm
          distp(k,1)=dexp(-dd1*r(k,k_1))
          enddo

          do i=indtmax,indtm
             distp(i,3)=distp(i,1)*r(i,k_1)
             distp(i,4)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d
             distp(i,5)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d
             distp(i,6)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d
             distp(i,7)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d
             distp(i,8)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d
          enddo 
         
          indorbp=indorb

          do ic=1,5
             if(iocc(indshell+ic).eq.1) then
                indorbp=indorbp+1
                do i=1,indtm
                z(indorbp,i)=distp(i,3+ic)*distp(i,3)
                enddo
             endif
          enddo

          if(indt.ne.1) then
             fun0=distp(0,3)
             fun=(1.d0-dd1*r(0,k_1))*distp(0,1)
             fun2=dd1*(dd1*r(0,k_1)-2.d0)*distp(0,1)

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                     z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0,k_1)
     &                *fun/r(0,k_1)
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d                      
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
           z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0,k_1)+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

           else
c               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                      z(indorbp,indt+i)=0.d0
                   enddo
                 endif
              enddo
               
           endif  !endif for r(0)

        endif   !endif for indt
c
        indpar=indpar+1
        indshell=indshell+5
        indorb=indorbp

            case(134) 
c     2p single exponential  r^3 e^{-z r}  ! 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=dexp(-dd2*r(k,k_1))
         enddo

         indorbp=indorb

         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)*r(i,k_1)**3
               enddo
            endif
         enddo

         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then
            fun0=distp(0,1)*r(0,k_1)**3
            fun=distp(0,1)*(3.d0-dd2*r(0,k_1))*r(0,k_1)
c           fun= derivative of fun0 respect to r divided dy r
            fun2=distp(0,1)*(dd2**2*r(0,k_1)**3-6*dd2*r(0,k_1)**2
     &      +6*r(0,k_1))
c           fun2= second derivative of fun0 respect to r
               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo
           else
               indorbp=indorb
               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo
           endif
         endif  !endif for indt
         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp


            case(135) 
c     2p single exponential  r^4 e^{-z r}  ! 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=dexp(-dd2*r(k,k_1))
         enddo

         indorbp=indorb

         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)*r(i,k_1)**4
               enddo
            endif
         enddo

         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then
            fun0=distp(0,1)*r(0,k_1)**4
            fun=distp(0,1)*(4.d0-dd2*r(0,k_1))*r(0,k_1)**2
            fun2=distp(0,1)*(12*r(0,k_1)**2-8*dd2*r(0,k_1)**3
     &      +dd2**2*r(0,k_1)**4)
               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo
           else
               indorbp=indorb
               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo
           endif
         endif  !endif for indt
         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp


            case(136) 
c     2p single exponential  r^5 e^{-z r}  ! 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=dexp(-dd2*r(k,k_1))
         enddo

         indorbp=indorb

         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)*r(i,k_1)**5
               enddo
            endif
         enddo

         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then
            fun0=distp(0,1)*r(0,k_1)**5
            fun=distp(0,1)*(5.d0-dd2*r(0,k_1))*r(0,k_1)**3
            fun2=distp(0,1)*(20*r(0,k_1)**3-10*dd2*r(0,k_1)**4
     &      +dd2**2*r(0,k_1)**5)
               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo
           else
               indorbp=indorb
               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo
           endif
         endif  !endif for indt
         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

       case(137)  
c     2s with cusp condition
c     dd1*(exp(-dd2*r)*(1+dd2*r))
         k_1=kion(1)

         dd2=dd(indpar+1)

!           if(iflagnorm.gt.2) then 
!           c=1.d0/dsqrt(1/4.d0/dd2**3+12*dd2/(2.d0*dd2)**4+
!    &3*dd2**2/4/dd2**5)/dsqrt(4.0*pi)
!           endif 

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=dexp(-dd2*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)*(1.d0+dd2*r(i,k_1))
               enddo
            endif


         if(indt.ne.1) then
               fun=-dd2**2*distp(0,1)
               fun2=fun*(1.d0-dd2*r(0,k_1))


                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*fun+fun2


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshellp
         indorb=indorbp

       case(138)  
c     2s with cusp condition
c     ( -dd2*r^2*exp(-dd2*r))  ! with no cusp condition der of 137 
         k_1=kion(1)

         dd2=dd(indpar+1)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=-dd2*dexp(-dd2*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)*r(i,k_1)**2
               enddo
            endif


         if(indt.ne.1) then
               fun=(2.d0-dd2*r(0,k_1))*distp(0,1)
               fun2=(2.d0-4*dd2*r(0,k_1)+(dd2*r(0,k_1))**2)*distp(0,1)
                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*fun+fun2


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshellp
         indorb=indorbp

       case(139)  
c     2s with cusp condition
c     ( r^3*exp(-dd2*r))  !  der of 128  
         k_1=kion(1)

         dd2=dd(indpar+1)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=-dexp(-dd2*r(k,k_1))*r(k,k_1)
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=distp(i,1)*r(i,k_1)**2
               enddo
            endif


         if(indt.ne.1) then
               fun=(3.d0-dd2*r(0,k_1))*distp(0,1)
               fun2=(6.d0-6*dd2*r(0,k_1)+(dd2*r(0,k_1))**2)*distp(0,1)
                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*fun+fun2


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshellp
         indorb=indorbp


       case(140)  
c     2p single exponential  -r e^{-z r}  ! der  of 121 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=-dexp(-dd2*r(k,k_1))
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)*r(i,k_1)
               enddo
            endif
         enddo


         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then 
            fun0=distp(0,1)*r(0,k_1)
            fun=distp(0,1)*(1.d0-dd2*r(0,k_1))/r(0,k_1)
            fun2=dd2*(dd2*r(0,k_1)-2.d0)*distp(0,1)
               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo
           else
               indorbp=indorb
               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo
           endif 
         endif  !endif for indt
         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

       case(141)  
c     2p single exponential  r^2  e^{-z r}  ! parent of 121 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=-dexp(-dd2*r(k,k_1))
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)*r(i,k_1)**2
               enddo
            endif
         enddo


         if(indt.ne.1) then

            fun0=distp(0,1)*r(0,k_1)**2
            fun=distp(0,1)*(2.d0-dd2*r(0,k_1))
            fun2=(2.d0-4.d0*dd2*r(0,k_1)+(dd2*r(0,k_1))**2)*distp(0,1)
               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo
         endif  !endif for indt
         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

       case(142)  ! der of 127 
!     4d  one parmater

          k_1=kion(1)
          dd1=dd(indpar+1)

          do k=indtmax,indtm
          distp(k,1)=dexp(-dd1*r(k,k_1))
          enddo

          do i=indtmax,indtm
             distp(i,3)=distp(i,1)*r(i,k_1)
             distp(i,4)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d
             distp(i,5)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d
             distp(i,6)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d
             distp(i,7)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d
             distp(i,8)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d
          enddo 
         
          indorbp=indorb

          do ic=1,5
             if(iocc(indshell+ic).eq.1) then
                indorbp=indorbp+1
                do i=1,indtm
                z(indorbp,i)=-distp(i,3+ic)*distp(i,3)
                enddo
             endif
          enddo

          if(indt.ne.1) then
             fun0=-distp(0,3)
             fun=-(1.d0-dd1*r(0,k_1))*distp(0,1)
             fun2=-dd1*(dd1*r(0,k_1)-2.d0)*distp(0,1)

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                     z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0,k_1)
     &                *fun/r(0,k_1)
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d                      
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
           z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0,k_1)+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

           else
c               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                      z(indorbp,indt+i)=0.d0
                   enddo
                 endif
              enddo
               
           endif  !endif for r(0)

        endif   !endif for indt
c
        indpar=indpar+1
        indshell=indshell+5
        indorb=indorbp

       case(143) 
!     4d  one parmater der of 133 

          k_1=kion(1)
          dd1=dd(indpar+1)

          do k=indtmax,indtm
          distp(k,1)=dexp(-dd1*r(k,k_1))
          enddo

          do i=indtmax,indtm
             distp(i,3)=distp(i,1)*r(i,k_1)**2
             distp(i,4)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d
             distp(i,5)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d
             distp(i,6)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d
             distp(i,7)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d
             distp(i,8)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d
          enddo 
         
          indorbp=indorb

          do ic=1,5
             if(iocc(indshell+ic).eq.1) then
                indorbp=indorbp+1
                do i=1,indtm
                z(indorbp,i)=-distp(i,3+ic)*distp(i,3)
                enddo
             endif
          enddo

          if(indt.ne.1) then
             fun0=-distp(0,3)
      fun=-(-2.d0+dd1*r(0,k_1))*distp(0,1)
      fun2=((dd1*r(0,k_1))**2 -4.d0*r(0,k_1)*dd1+2.d0)*distp(0,1)

               indorbp=indorb

               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                     z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0,k_1)
     &                *fun
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d                      
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
           z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

        endif   !endif for indt
c
        indpar=indpar+1
        indshell=indshell+5
        indorb=indorbp

            case(144) 
c     2p single exponential  -r^3 e^{-z r}  ! derivative of  130 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=-dexp(-dd2*r(k,k_1))
         enddo

         indorbp=indorb

         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)*r(i,k_1)**3
               enddo
            endif
         enddo

         if(indt.ne.1) then

            fun0=distp(0,1)*r(0,k_1)**3
            fun=distp(0,1)*(3.d0-dd2*r(0,k_1))*r(0,k_1)
c           fun= derivative of fun0 respect to r divided dy r
            fun2=distp(0,1)*(dd2**2*r(0,k_1)**3-6*dd2*r(0,k_1)**2
     &      +6*r(0,k_1))
c           fun2= second derivative of fun0 respect to r
               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo
         endif  !endif for indt
         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

       case(145)  
!     2s without cusp condition  !derivative 100
!     -(r^2*exp(-dd2*r^2))
         k_1=kion(1)

         dd2=dd(indpar+1)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=dexp(-dd2*r(k,k_1)*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=-distp(i,1)*r(i,k_1)**2
               enddo
            endif


         if(indt.ne.1) then
               fun0=dd2*r(0,k_1)**2
               fun=-2.d0*distp(0,1)*(1.d0-fun0)
               fun2=-2.d0*distp(0,1)*(1.d0-5.d0*fun0+2.d0*fun0**2)
  
                  do i=1,3
                  z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*fun+fun2


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshellp
         indorb=indorbp

       case(146)  
!     2p single exponential  -r^2  e^{-z r^2}  ! derivative  of 103 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=dexp(-dd2*r(k,k_1)*r(k,k_1))
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=-rmu(ic,i,k_1)*distp(i,1)*r(i,k_1)*r(i,k_1)
               enddo
            endif
         enddo


         if(indt.ne.1) then
            rp2=dd2*r(0,k_1)*r(0,k_1)
            fun0=-distp(0,1)*r(0,k_1)*r(0,k_1)
            fun=distp(0,1)*(-2.d0+2.d0*rp2)
            fun2=(-2.d0+10.d0*rp2-4.d0*rp2*rp2)*distp(0,1)
               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo
         endif  !endif for indt
         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

       case(147)
! 3d single gaussian 

          k_1=kion(1)
          dd1=dd(indpar+1)

          do k=indtmax,indtm
          distp(k,1)=dexp(-dd1*r(k,k_1)**2)
          enddo

          do i=indtmax,indtm
             distp(i,3)=distp(i,1)
             distp(i,4)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d
             distp(i,5)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d
             distp(i,6)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d
             distp(i,7)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d
             distp(i,8)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d
          enddo 
         
          indorbp=indorb

          do ic=1,5
             if(iocc(indshell+ic).eq.1) then
                indorbp=indorbp+1
                do i=1,indtm
                   z(indorbp,i)=distp(i,3+ic)*distp(i,3)
                enddo
             endif
          enddo

          if(indt.ne.1) then
             fun0=distp(0,3)
             fun=-2.d0*dd1*r(0,k_1)*distp(0,1)
             fun2=((2.d0*dd1*r(0,k_1))**2-2.d0*dd1)*distp(0,1)

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                     z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0,k_1)
     &                *fun/r(0,k_1)
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d           
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
           z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0,k_1)+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

           else
c               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                      z(indorbp,indt+i)=0.d0
                   enddo
                 endif
              enddo
               
           endif  !endif for r(0)

        endif   !endif for indt
c
        indpar=indpar+1
        indshell=indshell+5
        indorb=indorbp

       case(148)
!  derivative of 147 with respect to dd1

          k_1=kion(1)
          dd1=dd(indpar+1)

          do k=indtmax,indtm
          distp(k,1)=dexp(-dd1*r(k,k_1)**2)
          enddo

          do i=indtmax,indtm
             distp(i,3)=-r(i,k_1)**2*distp(i,1)
             distp(i,4)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d
             distp(i,5)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d
             distp(i,6)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d
             distp(i,7)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d
             distp(i,8)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d
          enddo 
         
          indorbp=indorb

          do ic=1,5
             if(iocc(indshell+ic).eq.1) then
                indorbp=indorbp+1
                do i=1,indtm
                   z(indorbp,i)=distp(i,3+ic)*distp(i,3)
                enddo
             endif
          enddo

          if(indt.ne.1) then
             fun0=distp(0,3)
             fun=2.d0*(dd1*r(0,k_1)**2-1.d0)*r(0,k_1)*distp(0,1)
             fun2=-2.d0*(2.d0*dd1**2*r(0,k_1)**4+1.d0
     &       -5.d0*dd1*r(0,k_1)**2)*distp(0,1)

            if(r(0,k_1).gt.1d-9) then
               indorbp=indorb

               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                     z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0,k_1)
     &                *fun/r(0,k_1)
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d                      
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
           z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0,k_1)+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

           else
c               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                      z(indorbp,indt+i)=0.d0
                   enddo
                 endif
              enddo
               
           endif  !endif for r(0)

        endif   !endif for indt
c
        indpar=indpar+1
        indshell=indshell+5
        indorb=indorbp

       case(149)
! derivative of 131 with respect z_1
! - r^4 exp(-z_1 r^2)

         k_1=kion(1)

         dd2=dd(indpar+1)

            indorbp=indorb+1
            indshellp=indshell+1 
               do k=indtmax,indtm
               distp(k,1)=dexp(-dd2*r(k,k_1)*r(k,k_1))
               enddo

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=-distp(i,1)*r(i,k_1)**4
               enddo
            endif


         if(indt.ne.1) then
         fun0=dd2*r(0,k_1)**2
         fun=-2.d0*r(0,k_1)**2*distp(0,1)*(2.d0-fun0)
         fun2=-2.d0*r(0,k_1)**2*distp(0,1)*(6.d0-9.d0*fun0+2.d0*fun0**2)
  
                  do i=1,3
                     z(indorbp,indt+i)=fun*rmu(i,0,k_1)
                  enddo

        z(indorbp,indt+4)=2.d0*fun+fun2


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshellp
         indorb=indorbp

       case(150)
c     2p single exponential  r e^{-z r^2}  

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=dexp(-dd2*r(k,k_1)**2)
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)*r(i,k_1)
               enddo
            endif
         enddo


         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then 
            fun0=distp(0,1)*r(0,k_1)
            cost=2.d0*dd2*r(0,k_1)**2
            fun=distp(0,1)*(1.d0-cost)/r(0,k_1)
            fun2=2.d0*dd2*fun0*(cost-3.d0)
               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo
           else
               indorbp=indorb
               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo
           endif 
         endif  !endif for indt
         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp

       case(151)
c     2p single exponential  -r^3  e^{-z r^2}  ! parent of 150 

         k_1=kion(1)

         dd2=dd(indpar+1)

         do k=indtmax,indtm
         distp(k,1)=dexp(-dd2*r(k,k_1)**2)
         enddo

         indorbp=indorb
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
               z(indorbp,i)=-rmu(ic,i,k_1)*distp(i,1)*r(i,k_1)**3
               enddo
            endif
         enddo


         if(indt.ne.1) then

          fun0=-distp(0,1)*r(0,k_1)**3
          cost=dd2*r(0,k_1)**2
          fun=distp(0,1)*(-3.d0+2.d0*cost)*r(0,k_1)
          fun2=-2.d0*distp(0,1)*r(0,k_1)*(3.d0-7.d0*cost+2.d0*cost**2)
               indorbp=indorb
               do ic=1,3
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,3
                       z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*
     &                       fun
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
                     enddo
                     z(indorbp,indt+4)=rmu(ic,0,k_1)
     $                *(4.d0*fun+fun2)
                  endif
               enddo
         endif  !endif for indt
         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp
          

       case(199) 
! derivative of 200 LA COSTANTE

            indorbp=indorb+1
            indshellp=indshell+1

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=0.d0
               enddo
            endif

         if(indt.ne.1) then
                  do i=1,3
                     z(indorbp,indt+i)=0.d0
                  enddo

        z(indorbp,indt+4)=0
         endif  !endif for indt

         indshell=indshellp
         indorb=indorbp


       case(200) 
!     THE  COSTANT

            indorbp=indorb+1
            indshellp=indshell+1

            if(iocc(indshellp).eq.1) then
               do i=1,indtm
               z(indorbp,i)=1.d0
               enddo
            endif

         if(indt.ne.1) then
                  do i=1,3
                     z(indorbp,indt+i)=0
                  enddo

        z(indorbp,indt+4)=0
         endif  !endif for indt

         indshell=indshellp
         indorb=indorbp


      case(1000:1099) 
!     s gaussian  r**(2*npower)*exp(-alpha*r**2)
 
         npower=iopt-1000
         k_1=kion(1)
         indorbp=indorb+1
         indshellp=indshell+1

         dd(indpar+1)=abs(dd(indpar+1))
         dd2=dd(indpar+1)
         do k=indtmax,indtm
            distp(k,1)=r(k,k_1)**(2*npower)*dexp(-dd2*r(k,k_1)**2)
         enddo
         
         if(iocc(indshellp).eq.1) then
            do i=1,indtm
               z(indorbp,i)=distp(i,1)
            enddo
         endif


         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then

            rp1=r(0,k_1)**2
            fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1
            fun2=(npower*(2.d0*npower-1.d0)-
     1      (1.d0+4.d0*npower)*dd2*rp1+2.d0*(dd2*rp1)**2)*
     1      distp(0,1)*2.d0/rp1

            if(iocc(indshellp).eq.1) then 
               do i=1,3
                  z(indorbp,indt+i)=rmu(i,0,k_1)*fun
               enddo
               z(indorbp,indt+4)=2.d0*fun+fun2
            endif
            
            else
               
            do i=1,4
               z(indorbp,indt+i)=0.d0
            enddo

            endif

         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+1
         indorb=indorbp

      case(2000:2099) 
!     s gaussian  -r**(2*(npower+1))*exp(-alpha*r**2) derivative of 1000
 
         npower=iopt+1-2000
         k_1=kion(1)
         indorbp=indorb+1
         indshellp=indshell+1


         dd(indpar+1)=abs(dd(indpar+1))
         dd2=dd(indpar+1)
         do k=indtmax,indtm
            distp(k,1)=-r(k,k_1)**(2*npower)*dexp(-dd2*r(k,k_1)**2)
         enddo
         
         if(iocc(indshellp).eq.1) then
            do i=1,indtm
               z(indorbp,i)=distp(i,1)
            enddo
         endif


         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then

            rp1=r(0,k_1)**2
            fun0=distp(0,1)
            fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1
            fun2=(npower*(2.d0*npower-1.d0)-
     1      (1.d0+4.d0*npower)*dd2*rp1+2.d0*(dd2*rp1)**2)*
     1      distp(0,1)*2.d0/rp1

            if(iocc(indshellp).eq.1) then 
               do i=1,3
                  z(indorbp,indt+i)=rmu(i,0,k_1)*fun
               enddo
               z(indorbp,indt+4)=2.d0*fun+fun2
            endif

            else
               
            do i=1,4
               z(indorbp,indt+i)=0.d0
            enddo

            endif

         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+1
         indorb=indorbp


      case(1100:1199)
!     p gaussian  r**(2*npower)*exp(-alpha*r**2)

         npower=iopt-1100
         k_1=kion(1)
         indorbp=indorb

         dd(indpar+1)=abs(dd(indpar+1))
         dd2=dd(indpar+1)
         do k=indtmax,indtm
            distp(k,1)=r(k,k_1)**(2*npower)*dexp(-dd2*r(k,k_1)**2)
         enddo
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
                  z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)
               enddo
            endif
         enddo


         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then

            rp1=r(0,k_1)**2
            fun0=distp(0,1)
            fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1
            fun2=(npower*(2.d0*npower-1.d0)-
     1      (1.d0+4.d0*npower)*dd2*rp1+2.d0*(dd2*rp1)**2)*
     1      distp(0,1)*2.d0/rp1

            indorbp=indorb
            do ic=1,3
               if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1
               do i=1,3
               z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*fun
               if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
               enddo
               z(indorbp,indt+4)=rmu(ic,0,k_1)*(4.d0*fun+fun2)
               endif
            enddo

           else
               indorbp=indorb
               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo
           endif 


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp         


      case(2100:2199)
!     p gaussian  -r**(2*(npower+1))*exp(-alpha*r**2) derivative 1100

         npower=iopt+1-2100
         k_1=kion(1)
         indorbp=indorb

         dd(indpar+1)=abs(dd(indpar+1))
         dd2=dd(indpar+1)
         do k=indtmax,indtm
            distp(k,1)=-r(k,k_1)**(2*npower)*dexp(-dd2*r(k,k_1)**2)
         enddo
         
         do ic=1,3
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
                  z(indorbp,i)=rmu(ic,i,k_1)*distp(i,1)
               enddo
            endif
         enddo


         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then

            rp1=r(0,k_1)**2
            fun0=distp(0,1)
            fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1
            fun2=(npower*(2.d0*npower-1.d0)-
     1      (1.d0+4.d0*npower)*dd2*rp1+2.d0*(dd2*rp1)**2)*
     1      distp(0,1)*2.d0/rp1

            indorbp=indorb
            do ic=1,3
               if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1
               do i=1,3
               z(indorbp,indt+i)=rmu(ic,0,k_1)*rmu(i,0,k_1)*fun
               if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
               enddo
               z(indorbp,indt+4)=rmu(ic,0,k_1)*(4.d0*fun+fun2)
               endif
            enddo

           else
               indorbp=indorb
               do ic=1,3
                  if(iocc(indshell+ic).eq.1) then
                     indorbp=indorbp+1
                     do i=1,4
                     z(indorbp,indt+i)=0.d0
                     enddo
                  endif
               enddo
           endif 


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+3
         indorb=indorbp


      case(1200:1299)
!     d gaussian  r**(2*npower)*exp(-alpha*r**2)

         npower=iopt-1200
         k_1=kion(1)
         indorbp=indorb

         dd(indpar+1)=abs(dd(indpar+1))
         dd2=dd(indpar+1)
         do k=indtmax,indtm
            distp(k,1)=r(k,k_1)**(2*npower)*dexp(-dd2*r(k,k_1)**2)
         enddo

         do i=indtmax,indtm
            distp(i,2)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d
            distp(i,3)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d
            distp(i,4)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d
            distp(i,5)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d
            distp(i,6)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d
         enddo
         
         do ic=1,5
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
                  z(indorbp,i)=distp(i,1+ic)*distp(i,1)
               enddo
            endif
         enddo

         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then

            rp1=r(0,k_1)**2
            fun0=distp(0,1)
            fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1
            fun2=(npower*(2.d0*npower-1.d0)-
     1      (1.d0+4.d0*npower)*dd2*rp1+2.d0*(dd2*rp1)**2)*
     1      distp(0,1)*2.d0/rp1


               indorbp=indorb
               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0,k_1)
     &                *fun
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d                 
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
      z(indorbp,indt+4)=distp(0,1+ic)*(6.d0*fun+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

           else
               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                      z(indorbp,indt+i)=0.d0
                   enddo
                 endif
              enddo
               
           endif  !endif for r(0)


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+5
         indorb=indorbp


      case(2200:2299)
!     d gaussian  -r**(2*(npower+1))*exp(-alpha*r**2) derivative 1200

         npower=iopt+1-2200
         k_1=kion(1)
         indorbp=indorb

         dd(indpar+1)=abs(dd(indpar+1))
         dd2=dd(indpar+1)
         do k=indtmax,indtm
            distp(k,1)=-r(k,k_1)**(2*npower)*dexp(-dd2*r(k,k_1)**2)
         enddo

         do i=indtmax,indtm
            distp(i,2)=(3.d0*rmu(3,i,k_1)**2-r(i,k_1)**2)*cost1d
            distp(i,3)=(rmu(1,i,k_1)**2-rmu(2,i,k_1)**2)*cost2d
            distp(i,4)=rmu(1,i,k_1)*rmu(2,i,k_1)*cost3d
            distp(i,5)=rmu(2,i,k_1)*rmu(3,i,k_1)*cost3d
            distp(i,6)=rmu(1,i,k_1)*rmu(3,i,k_1)*cost3d
         enddo
         
         do ic=1,5
            if(iocc(indshell+ic).eq.1) then
               indorbp=indorbp+1
               do i=1,indtm
                  z(indorbp,i)=distp(i,1+ic)*distp(i,1)
               enddo
            endif
         enddo

         if(indt.ne.1) then

            if(r(0,k_1).gt.1d-9) then

            rp1=r(0,k_1)**2
            fun0=distp(0,1)
            fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1
            fun2=(npower*(2.d0*npower-1.d0)-
     1      (1.d0+4.d0*npower)*dd2*rp1+2.d0*(dd2*rp1)**2)*
     1      distp(0,1)*2.d0/rp1


               indorbp=indorb
               do ic=1,5
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1
                     do i=1,3
                z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0,k_1)
     &                *fun
                       if(ic.eq.1) then
                          if(i.ne.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost1d
                          else
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            4.d0*rmu(i,0,k_1)*fun0*cost1d
                          endif
                       elseif(ic.eq.2) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d            
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)-
     &                            2.d0*rmu(i,0,k_1)*fun0*cost2d                 
                          endif
                       elseif(ic.eq.3) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          elseif(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.4) then
                          if(i.eq.2) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(2,0,k_1)*fun0*cost3d
                          endif
                       elseif(ic.eq.5) then
                          if(i.eq.1) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(3,0,k_1)*fun0*cost3d
                          elseif(i.eq.3) then
                             z(indorbp,indt+i)=z(indorbp,indt+i)+
     &                            rmu(1,0,k_1)*fun0*cost3d
                          endif !endif for i
                       endif    !endif for ic 
                    enddo  !enddo for i
      z(indorbp,indt+4)=distp(0,1+ic)*(6.d0*fun+fun2)
                 endif  !endif for iocc
             enddo     ! enddo fot ic

           else
               
              indorbp=indorb 
              
              do ic=1,5
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1
                    do i=1,4
                      z(indorbp,indt+i)=0.d0
                   enddo
                 endif
              enddo
               
           endif  !endif for r(0)


         endif  !endif for indt

         indpar=indpar+1
         indshell=indshell+5
         indorb=indorbp


      case default
      write(6,*) 'WARNING makefun: orbital',iopt,'not found'
      stop


      end select 
! ** ** ** ** ** ** ** END OF JASTROW ORBITALS ** ** ** ** ** ** ** ** *
 
      return
      end


