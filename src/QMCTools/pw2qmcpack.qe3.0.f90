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
    
    



!
! Copyright (C) 2004 PWSCF group 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
!----------------------------------------------------------------------- 
PROGRAM pw2qmcpack
  !----------------------------------------------------------------------- 

  ! This subroutine writes the file "prefix".pwscf.xml and "prefix".pwscf.h5
  ! containing the  plane wave coefficients and other stuff needed by the 
  ! QMC code QMCPACK. 

#include "f_defs.h"

  USE io_files,  ONLY : nd_nmbr, prefix, outdir, tmp_dir, trimcheck
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  INTEGER :: ios

  NAMELIST / inputpp / prefix, outdir

  CALL start_postproc(nd_nmbr)
  ! 
  !   set default values for variables in namelist 
  ! 
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'

  IF ( ionode )  THEN 
     !
     READ (5, inputpp, err=200, iostat=ios)
200  CALL errore('pw2qmcpack', 'reading inputpp namelist', ABS(ios))
     tmp_dir = trimcheck (outdir)
     !
  END IF
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( prefix, ionode_id ) 
  CALL mp_bcast(tmp_dir, ionode_id ) 
  !
  CALL read_file
  CALL openfil_pp
  !
  CALL compute_qmcpack
  !
  CALL stop_pp
  STOP

END PROGRAM pw2qmcpack


SUBROUTINE compute_qmcpack

  USE kinds, ONLY: DP
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
  USE cell_base, ONLY: omega, alat, tpiba2, at, bg
  USE char, ONLY: title
  USE constants, ONLY: tpi
  USE ener, ONLY: ewld, ehart, etxc, vtxc, etot, etxcc
  USE gvect, ONLY: ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
       nrxx, g, gg, ecutwfc, gcutm, nl, igtongl
  USE klist , ONLY: nks, nelec, nelup, neldw, xk
  USE lsda_mod, ONLY: lsda, nspin
  USE scf, ONLY: rho, rho_core, rhog, rhog_core
  USE vlocal, ONLY: vloc, vnew, strf
  USE wvfct, ONLY: npw, npwx, nbnd, gamma_only, igk, g2kin, wg, et
  USE uspp, ONLY: nkb, vkb, dvan
  USE uspp_param, ONLY: nh
  USE becmod,   ONLY: becp 
  USE io_global, ONLY: stdout
  USE io_files, ONLY: nd_nmbr, nwordwfc, iunwfc, iun => iunsat, tmp_dir, prefix
  USE wavefunctions_module, ONLY : evc, psic
  use gsmooth,         only: nls, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
  use iotk_module
  use iotk_xtox_interf
  IMPLICIT NONE
  INTEGER :: ig, ibnd, ik, io, na, j, ispin, nbndup, nbnddown, &
       nk, ngtot, ig7, ikk, nt, ijkb0, ikb, ih, jh, jkb, at_num, &
       nelec_tot, nelec_up, nelec_down, ii, igx, igy, igz
  INTEGER, ALLOCATABLE :: INDEX(:), igtog(:)
  LOGICAL :: exst, found
  REAL(DP) :: ek, eloc, enl, charge, etotefield
  REAL(DP) :: bg_qmc(3,3), tau_r(3), g_red(3)
  COMPLEX(DP), ALLOCATABLE :: phase(:),aux(:), hpsi(:,:), eigpacked(:)
  COMPLEX(DP), ALLOCATABLE :: psitr(:)
  REAL(DP), ALLOCATABLE ::  g_cart(:,:),psireal(:)
  INTEGER :: ios, ierr, h5len,ig_c,save_complex, nup,ndown
  INTEGER, EXTERNAL :: atomic_number, is_complex
  INTEGER, ALLOCATABLE::  g_qmc(:,:)
  REAL (DP), EXTERNAL :: ewald

  CHARACTER(256)          :: tmp,h5name,eigname
  CHARACTER(iotk_attlenx) :: attr

  CALL init_us_1
  CALL newd
  !####io = 77

  !####WRITE (6,'(/,5x,''Writing file pwfn.data for program CASINO'')')

  !####CALL seqopn( 77, 'pwfn.data', 'formatted',exst)  

  !ALLOCATE (hpsi(npwx, nbnd))
  ALLOCATE (aux(nrxx))
  ALLOCATE (becp (nkb,nbnd))
  ! four times npwx should be enough
  ALLOCATE (INDEX (4*npwx) )
  ALLOCATE (igtog (4*npwx) )


  !hpsi (:,:) = (0.d0, 0.d0)
  INDEX(:) = 0
  igtog(:) = 0

  IF( lsda ) THEN
     nbndup = nbnd
     nbnddown = nbnd
     nk = nks/2
     !     nspin = 2
  ELSE
     nbndup = nbnd
     nbnddown = 0
     nk = nks
     !     nspin = 1
  ENDIF

  !  if(nks > 1) rewind(iunigk)
  !  do ik=1,nks
  !     if(nks > 1) read(iunigk) npw, igk
  !     
  !  if(nks > 1) rewind(iunigk)
  ek  = 0.d0
  eloc= 0.d0
  enl = 0.d0
  !
  DO ispin = 1, nspin 
     !
     !     calculate the local contribution to the total energy
     !
     !      bring rho to G-space
     !
     aux(:) = CMPLX ( rho(:,ispin), 0.d0)
     CALL cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
     !
     DO nt=1,ntyp
        DO ig = 1, ngm
           eloc = eloc + vloc(igtongl(ig),nt) * strf(ig,nt) &
                * CONJG(aux(nl(ig)))
        ENDDO
     ENDDO

     DO ik = 1, nk
        ikk = ik + nk*(ispin-1)
        CALL gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        CALL davcio (evc, nwordwfc, iunwfc, ikk, - 1)
        CALL init_us_2 (npw, igk, xk (1, ikk), vkb)
        CALL ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)

        DO ig =1, npw
           IF( igk(ig) > 4*npwx ) & 
                CALL errore ('pw2qmcpack','increase allocation of index', ig)
           INDEX( igk(ig) ) = 1
        ENDDO
     ENDDO
  ENDDO

#ifdef __PARA
  CALL reduce(1,eloc)
  CALL reduce(1,ek)
  CALL poolreduce(1,ek)
  CALL poolreduce(1,enl)
#endif
  eloc = eloc * omega 
  ek = ek * tpiba2

  ngtot = 0
  DO ig = 1, 4*npwx
     IF( INDEX(ig) == 1 ) THEN
        ngtot = ngtot + 1
        igtog(ngtot) = ig
     ENDIF
  ENDDO

  ALLOCATE (g_qmc(3,ngtot))
  ALLOCATE (g_cart(3,ngtot))

  ! get the number of electrons
  nelec_tot= NINT(nelec)
  nup=NINT(nelup)
  ndown=NINT(neldw)

  if(nup .eq. 0) then
    ndown=nelec_tot/2
    nup=nelec_tot-ndown
  endif

  h5name = TRIM( prefix ) // '.pwscf.h5'
  eigname = "eigenstates_"//trim(iotk_itoa(nr1s))//'_'//trim(iotk_itoa(nr2s))//'_'//trim(iotk_itoa(nr3s))

  bg_qmc(:,:)=bg(:,:)/alat

  !! create a file for particle set
  tmp = TRIM( tmp_dir ) // TRIM( prefix )// '.ptcl.xml'
  CALL iotk_open_write(iun, FILE=TRIM(tmp), ROOT="qmcsystem", IERR=ierr )

  CALL iotk_write_attr (attr,"name","global",first=.true.)
  CALL iotk_write_begin(iun, "simulationcell",ATTR=attr)
  CALL iotk_write_attr (attr,"name","lattice",first=.true.)
  CALL iotk_write_attr (attr,"units","bohr")
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  WRITE(iun,100) alat*at(1,1), alat*at(2,1), alat*at(3,1)
  WRITE(iun,100) alat*at(1,2), alat*at(2,2), alat*at(3,2)
  WRITE(iun,100) alat*at(1,3), alat*at(2,3), alat*at(3,3)
  CALL iotk_write_end(iun, "parameter")
  CALL iotk_write_attr (attr,"name","reciprocal",first=.true.)
  CALL iotk_write_attr (attr,"units","2pi/bohr")
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  WRITE(iun,100) bg_qmc(1,1), bg_qmc(2,1), bg_qmc(3,1)
  WRITE(iun,100) bg_qmc(1,2), bg_qmc(2,2), bg_qmc(3,2)
  WRITE(iun,100) bg_qmc(1,3), bg_qmc(2,3), bg_qmc(3,3)
  CALL iotk_write_end(iun, "parameter")

  CALL iotk_write_attr (attr,"name","bconds",first=.true.)
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  WRITE(iun,'(a)') 'p p p'
  CALL iotk_write_end(iun, "parameter")

  CALL iotk_write_attr (attr,"name","LR_dim_cutoff",first=.true.)
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  WRITE(iun,'(a)') '10'
  CALL iotk_write_end(iun, "parameter")
  CALL iotk_write_end(iun, "simulationcell")

  ! <particleset name="ions">
  CALL iotk_write_attr (attr,"name","ion0",first=.true.)
  CALL iotk_write_attr (attr,"size",nat)
  CALL iotk_write_begin(iun, "particleset",ATTR=attr)

  ! ionic species --> group
  DO na=1,ntyp
  CALL iotk_write_attr (attr,"name",TRIM(atm(na)),first=.true.)
  CALL iotk_write_begin(iun, "group",ATTR=attr)

  CALL iotk_write_attr (attr,"name","charge",first=.true.)
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  write(iun,*) zv(na)
  CALL iotk_write_end(iun, "parameter")

  CALL iotk_write_end(iun, "group")
  ENDDO

  ! <attrib name="ionid"/>
  CALL iotk_write_attr (attr,"name","ionid",first=.true.)
  CALL iotk_write_attr (attr,"datatype","stringArray")
  CALL iotk_write_begin(iun, "attrib",ATTR=attr)
  write(iun,'(a)') (TRIM(atm(ityp(na))),na=1,nat)
  CALL iotk_write_end(iun, "attrib")

  ! <attrib name="position"/>
  CALL iotk_write_attr (attr,"name","position",first=.true.)
  CALL iotk_write_attr (attr,"datatype","posArray")
  CALL iotk_write_attr (attr,"condition","0")
  CALL iotk_write_begin(iun, "attrib",ATTR=attr)
  ! write in cartesian coordinates in bohr
  ! problem with xyz ordering inrelation to real-space grid
  DO na = 1, nat
  !tau_r(1)=alat*(tau(1,na)*bg_qmc(1,1)+tau(2,na)*bg_qmc(1,2)+tau(3,na)*bg_qmc(1,3))
  !tau_r(2)=alat*(tau(1,na)*bg_qmc(2,1)+tau(2,na)*bg_qmc(2,2)+tau(3,na)*bg_qmc(2,3))
  !tau_r(3)=alat*(tau(1,na)*bg_qmc(3,1)+tau(2,na)*bg_qmc(3,2)+tau(3,na)*bg_qmc(3,3))
  tau_r(1)=alat*tau(1,na)
  tau_r(2)=alat*tau(2,na)
  tau_r(3)=alat*tau(3,na)
  WRITE(iun,100) (tau_r(j),j=1,3)
  ENDDO
  !write(iun,100) tau
  CALL iotk_write_end(iun, "attrib")
  CALL iotk_write_end(iun, "particleset")
  ! </particleset>

  ! <particleset name="e">
  CALL iotk_write_attr (attr,"name","e",first=.true.)
  CALL iotk_write_attr (attr,"random","yes")
  CALL iotk_write_begin(iun, "particleset",ATTR=attr)

  ! <group name="u" size="" >
  CALL iotk_write_attr (attr,"name","u",first=.true.)
  CALL iotk_write_attr (attr,"size",nup)
  CALL iotk_write_begin(iun, "group",ATTR=attr)
  CALL iotk_write_attr (attr,"name","charge",first=.true.)
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  write(iun,*) -1
  CALL iotk_write_end(iun, "parameter")
  CALL iotk_write_end(iun, "group")

  ! <group name="d" size="" >
  CALL iotk_write_attr (attr,"name","d",first=.true.)
  CALL iotk_write_attr (attr,"size",ndown)
  CALL iotk_write_begin(iun, "group",ATTR=attr)
  CALL iotk_write_attr (attr,"name","charge",first=.true.)
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  write(iun,*) -1
  CALL iotk_write_end(iun, "parameter")
  CALL iotk_write_end(iun, "group")
  CALL iotk_write_end(iun, "particleset")
  CALL iotk_close_write(iun)

  !! close the file
  DO ik = 0, nk-1
  ! create a xml input file for each k-point
     IF(nk .gt. 1) THEN
       tmp = TRIM( tmp_dir ) // TRIM( prefix ) //TRIM(iotk_index(ik))// '.wfs.xml'
     ELSE
       tmp = TRIM( tmp_dir ) // TRIM( prefix )// '.wfs.xml'
     ENDIF
     CALL iotk_open_write(iun, FILE=TRIM(tmp), ROOT="qmcsystem", IERR=ierr )
     ! <wavefunction name="psi0">
     CALL iotk_write_attr (attr,"name","psi0",first=.true.)
     CALL iotk_write_attr (attr,"target","e")
     CALL iotk_write_begin(iun, "wavefunction",ATTR=attr)
       write(iun,'(a)') '<!-- Uncomment this out to use plane-wave basis functions'
       CALL iotk_write_attr (attr,"type","PW",first=.true.)
       CALL iotk_write_attr (attr,"href",TRIM(h5name))
       CALL iotk_write_attr (attr,"version","1.10")
       CALL iotk_write_begin(iun, "determinantset",ATTR=attr)
       write(iun,'(a)') '--> '
       CALL iotk_write_attr (attr,"type","bspline",first=.true.)
       CALL iotk_write_attr (attr,"href",TRIM(h5name))
       CALL iotk_write_attr (attr,"version","0.10")
       CALL iotk_write_begin(iun, "determinantset",ATTR=attr)
          CALL iotk_write_attr (attr,"ecut",ecutwfc/2,first=.true.)
          ! basisset to overwrite cutoff to a smaller value
          CALL iotk_write_begin(iun, "basisset",ATTR=attr)
             ! add grid to use spline on FFT grid
             CALL iotk_write_attr (attr,"dir","0",first=.true.)
             CALL iotk_write_attr (attr,"npts",nr1s)
             CALL iotk_write_attr (attr,"closed","no")
             CALL iotk_write_empty(iun, "grid",ATTR=attr)
             CALL iotk_write_attr (attr,"dir","1",first=.true.)
             CALL iotk_write_attr (attr,"npts",nr2s)
             CALL iotk_write_attr (attr,"closed","no")
             CALL iotk_write_empty(iun, "grid",ATTR=attr)
             CALL iotk_write_attr (attr,"dir","2",first=.true.)
             CALL iotk_write_attr (attr,"npts",nr3s)
             CALL iotk_write_attr (attr,"closed","no")
             CALL iotk_write_empty(iun, "grid",ATTR=attr)
          CALL iotk_write_end(iun, "basisset")
          
          !CALL iotk_write_attr (attr,"href",TRIM(h5name),first=.true.)
          !CALL iotk_write_empty(iun, "coefficients",ATTR=attr)
  
          ! write the index of the twist angle
          CALL iotk_write_attr (attr,"name","twistIndex",first=.true.)
          CALL iotk_write_begin(iun, "h5tag",ATTR=attr)
          write(iun,*) ik
          CALL iotk_write_end(iun, "h5tag")

          CALL iotk_write_attr (attr,"name","twistAngle",first=.true.)
          CALL iotk_write_begin(iun, "h5tag",ATTR=attr)
          g_red(1)=at(1,1)*xk(1,ik+1)+at(2,1)*xk(2,ik+1)+at(3,1)*xk(3,ik+1)
          g_red(2)=at(1,2)*xk(1,ik+1)+at(2,2)*xk(2,ik+1)+at(3,2)*xk(3,ik+1)
          g_red(3)=at(1,3)*xk(1,ik+1)+at(2,3)*xk(2,ik+1)+at(3,3)*xk(3,ik+1)
          !write(iun,100) xk(1,ik+1),xk(2,ik+1),xk(3,ik+1)
          write(iun,100) g_red(1),g_red(2),g_red(3)
          CALL iotk_write_end(iun, "h5tag")
          !write(iun,'(a)') '<!-- Uncomment this out for bspline wavefunctions '
          CALL iotk_write_attr (attr,"name","eigenstates",first=.true.)
          CALL iotk_write_begin(iun, "h5tag",ATTR=attr)
          write(iun,'(a)') TRIM(eigname)
          CALL iotk_write_end(iun, "h5tag")
          !write(iun,'(a)') '--> '

  
          CALL iotk_write_begin(iun, "slaterdeterminant")
             ! build determinant for up electrons
             CALL iotk_write_attr (attr,"id","updet",first=.true.)
             CALL iotk_write_attr (attr,"size",nup)
             CALL iotk_write_begin(iun, "determinant",ATTR=attr)
                CALL iotk_write_attr (attr,"mode","ground",first=.true.)
                CALL iotk_write_attr (attr,"spindataset",0)
                CALL iotk_write_begin(iun, "occupation",ATTR=attr)
                CALL iotk_write_end(iun, "occupation")
             CALL iotk_write_end(iun, "determinant")
  
             ! build determinant for down electrons
             CALL iotk_write_attr (attr,"id","downdet",first=.true.)
             CALL iotk_write_attr (attr,"size",ndown)
             IF( lsda ) CALL iotk_write_attr (attr,"ref","updet")
             CALL iotk_write_begin(iun, "determinant",ATTR=attr)
               CALL iotk_write_attr (attr,"mode","ground",first=.true.)
               IF( lsda ) THEN
                 CALL iotk_write_attr (attr,"spindataset",1)
               ELSE
                 CALL iotk_write_attr (attr,"spindataset",0)
               ENDIF
               CALL iotk_write_begin(iun, "occupation",ATTR=attr)
               CALL iotk_write_end(iun, "occupation")
             CALL iotk_write_end(iun, "determinant")
          CALL iotk_write_end(iun, "slaterdeterminant")
  
       CALL iotk_write_end(iun, "determinantset")
     CALL iotk_write_end(iun, "wavefunction")
  
     CALL iotk_close_write(iun)
  ENDDO

  tmp = TRIM( tmp_dir )//TRIM( prefix ) //'.pwscf.h5'
  h5len = LEN_TRIM(tmp)

  DO ig=1, ngtot
    ig_c =igtog(ig)
    g_cart(1,ig)=tpi/alat*g(1,ig_c)
    g_cart(2,ig)=tpi/alat*g(2,ig_c)
    g_cart(3,ig)=tpi/alat*g(3,ig_c)
    g_qmc(1,ig)=at(1,1)*g(1,ig_c)+at(2,1)*g(2,ig_c)+at(3,1)*g(3,ig_c)
    g_qmc(2,ig)=at(1,2)*g(1,ig_c)+at(2,2)*g(2,ig_c)+at(3,2)*g(3,ig_c)
    g_qmc(3,ig)=at(1,3)*g(1,ig_c)+at(2,3)*g(2,ig_c)+at(3,3)*g(3,ig_c)
  ENDDO


  ! generate hdf5 file containing all the parameters
  ! print out 2*occupied bands
  CALL pwhdf_open_file(tmp,h5len)
  !CALL pwhdf_write_basis(g,igtog,ngtot)
  CALL pwhdf_write_basis(g_qmc,g_cart,ngtot)
  CALL pwhdf_write_parameters(nelec_tot,nspin,nbnd,nk,ecutwfc/2,alat,at)
  !
100 FORMAT (3(1x,f20.15))

  ALLOCATE (eigpacked(ngtot))
  
  ! start a main section to save eigen values and vector
  CALL pwhdf_open_eigg
  DO ik = 1, nk
     g_red(1)=at(1,1)*xk(1,ik)+at(2,1)*xk(2,ik)+at(3,1)*xk(3,ik)
     g_red(2)=at(1,2)*xk(1,ik)+at(2,2)*xk(2,ik)+at(3,2)*xk(3,ik)
     g_red(3)=at(1,3)*xk(1,ik)+at(2,3)*xk(2,ik)+at(3,3)*xk(3,ik)
     CALL pwhdf_open_twist(ik,g_red(1),nbnd,nspin)
     !CALL pwhdf_open_twist(ik,xk(1,ik),2*nbnd,nspin)
     DO ispin = 1, nspin 
        ikk = ik + nk*(ispin-1)
        IF( nks > 1 ) THEN
           CALL gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
           CALL davcio(evc,nwordwfc,iunwfc,ikk,-1)
        ENDIF
        DO ibnd = 1, nbnd
           DO ig=1, ngtot
              ! now for all G vectors find the PW coefficient for this k-point
              found = .FALSE.
              DO ig7 = 1, npw
                 IF( igk(ig7) == igtog(ig) )THEN
                    eigpacked(ig)=evc(ig7,ibnd)
                    found = .TRUE.
                    GOTO 17
                 ENDIF
              ENDDO
              ! if can't find the coefficient this is zero
17            IF( .NOT. found ) eigpacked(ig)=(0.d0,0.d0) 
           ENDDO
           CALL pwhdf_write_band(ibnd,ispin,et(ibnd,ikk)/2,eigpacked,ngtot)
        ENDDO
     ENDDO
     CALL pwhdf_close_twist
  ENDDO
  CALL pwhdf_close_eigg

  !ALLOCATE (phase(nrxxs))

  ALLOCATE(psireal(nrx1s*nrx2s*nrx3s))
  ALLOCATE(psitr(nrx1s*nrx2s*nrx3s))

  ! open real-space wavefunction on FFT grid
  CALL pwhdf_open_eigr(nr1s,nr2s,nr3s)
  DO ik = 1, nk
     !! evaluate the phase
     !phase(:) = (0.d0,0.d0)
     !if ( ig_(ik,ib)>0) phase( nls(ig_(ik,ib)) ) = (1.d0,0.d0)
     !call cft3s (phase, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
     g_red(1)=at(1,1)*xk(1,ik)+at(2,1)*xk(2,ik)+at(3,1)*xk(3,ik)
     g_red(2)=at(1,2)*xk(1,ik)+at(2,2)*xk(2,ik)+at(3,2)*xk(3,ik)
     g_red(3)=at(1,3)*xk(1,ik)+at(2,3)*xk(2,ik)+at(3,3)*xk(3,ik)

     IF(g_red(1)*g_red(1)+g_red(2)*g_red(2)+g_red(3)*g_red(3)<1e-12) THEN
       save_complex=0
     ELSE
       save_complex=1
     ENDIF

     ! open a twsit
     CALL pwhdf_open_twist(ik,g_red(1),nbnd,nspin)
     DO ispin = 1, nspin 
        ikk = ik + nk*(ispin-1)
        IF( nks > 1 ) THEN
           CALL gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
           CALL davcio(evc,nwordwfc,iunwfc,ikk,-1)
        ENDIF
        DO ibnd = 1, nbnd !!transform G to R
           psic(:)=(0.d0,0.d0)
           psic(nls(igk(1:npw)))=evc(1:npw,ibnd)
           call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)

           !! THIS IS ORIGINAL: changes are made to convert real and reorder for C arrays
           !!CALL pwhdf_write_wfr(ibnd,ispin,et(ibnd,ikk)/2,psic,1)
           IF(save_complex .eq. 1) THEN
             !psic(:)=psic(:)/omega
             !CALL pwhdf_write_wfr(ibnd,ispin,et(ibnd,ikk)/2,psic,save_complex)
             ii=1
             DO igx=1,nr1s
             DO igy=0,nr2s-1
             DO igz=0,nr3s-1
             psitr(ii)=psic(igx+nr1s*(igy+igz*nr2s))/omega
             ii=ii+1
             ENDDO
             ENDDO
             ENDDO
             CALL pwhdf_write_wfr(ibnd,ispin,et(ibnd,ikk)/2,psitr,save_complex)
           ELSE
             !DO ig=1,nr1s*nr2s*nr3s
             !  psireal(ig)=real(psic(ig))
             !ENDDO
             !psireal(1:nr1s*nr2s*nr3s)=real(psic(1:nr1s*nr2s*nr3s))/omega
             ii=1
             DO igx=1,nr1s
             DO igy=0,nr2s-1
             DO igz=0,nr3s-1
             psireal(ii)=real(psic(igx+nr1s*(igy+igz*nr2s)))/omega
             ii=ii+1
             ENDDO
             ENDDO
             ENDDO
             CALL pwhdf_write_wfr(ibnd,ispin,et(ibnd,ikk)/2,psireal,save_complex)
           ENDIF
           !! conversion and output complete for each band
        ENDDO
     ENDDO
     CALL pwhdf_close_twist
  ENDDO
  CALL pwhdf_close_eigr

  ! write charge density
  ! ignore spin for the time being
  !CALL pwhdf_write_rho(rho,rhog(1,1),ngm)
  ii=1
  DO igx=1,nr1s
  DO igy=0,nr2s-1
  DO igz=0,nr3s-1
  psireal(ii)=rho(igx+nr1s*(igy+igz*nr2s),1)
  ii=ii+1
  ENDDO
  ENDDO
  ENDDO
  CALL pwhdf_write_rho(psireal,rhog(1,1),ngm)

  ! close the file
  CALL pwhdf_close_file

  DEALLOCATE (igtog)
  DEALLOCATE (index)
  DEALLOCATE (becp)
  DEALLOCATE (aux)
  !DEALLOCATE (hpsi)
  DEALLOCATE (eigpacked)
  DEALLOCATE (g_qmc)
  DEALLOCATE (g_cart)
  DEALLOCATE (psireal)
  DEALLOCATE (psitr)
  !DEALLOCATE (phase)

END SUBROUTINE compute_qmcpack


