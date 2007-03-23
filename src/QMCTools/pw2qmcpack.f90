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
  USE klist , ONLY: nks, nelec, xk
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
  USE iotk_module
  IMPLICIT NONE
  INTEGER :: ig, ibnd, ik, io, na, j, ispin, nbndup, nbnddown, &
       nk, ngtot, ig7, ikk, nt, ijkb0, ikb, ih, jh, jkb, at_num, &
       nelec_tot, nelec_up, nelec_down
  INTEGER, ALLOCATABLE :: INDEX(:), igtog(:)
  LOGICAL :: exst, found
  REAL(DP) :: ek, eloc, enl, charge, etotefield
  COMPLEX(DP), ALLOCATABLE :: aux(:), hpsi(:,:), eigpacked(:), aux_ext(:,:,:)
  INTEGER :: ios, ierr, h5len
  INTEGER, EXTERNAL :: atomic_number
  REAL (DP), EXTERNAL :: ewald

  CHARACTER(256)          :: tmp,h5name
  CHARACTER(iotk_attlenx) :: attr

  CALL init_us_1
  CALL newd
  !####io = 77

  !####WRITE (6,'(/,5x,''Writing file pwfn.data for program CASINO'')')

  !####CALL seqopn( 77, 'pwfn.data', 'formatted',exst)  

  ALLOCATE (hpsi(npwx, nbnd))
  ALLOCATE (aux(nrxx))
  ALLOCATE (becp (nkb,nbnd))
  ! four times npwx should be enough
  ALLOCATE (INDEX (4*npwx) )
  ALLOCATE (igtog (4*npwx) )

  hpsi (:,:) = (0.d0, 0.d0)
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

  ! get the number of electrons
  nelec_tot= NINT(nelec)

  h5name = TRIM( prefix ) // '.pwscf.h5'

  ! create a xml input file for each k-point
  DO ik = 0, nk-1
     IF(nk .gt. 1) THEN
       tmp = TRIM( tmp_dir ) // TRIM( prefix ) //TRIM(iotk_index(ik))// '.pwscf.xml'
     ELSE
       tmp = TRIM( tmp_dir ) // TRIM( prefix )// '.pwscf.xml'
     ENDIF
     CALL iotk_open_write(iun, FILE=TRIM(tmp), ROOT="qmcsystem", IERR=ierr )
  
     CALL iotk_write_attr (attr,"name","global",first=.true.)
     CALL iotk_write_begin(iun, "simulationcell",ATTR=attr)
        CALL iotk_write_attr (attr,"name","lattice",first=.true.)
        CALL iotk_write_begin(iun, "parameter",ATTR=attr)
        WRITE(iun,100) alat*at(1,1), alat*at(2,1), alat*at(3,1)
        WRITE(iun,100) alat*at(1,2), alat*at(2,2), alat*at(3,2)
        WRITE(iun,100) alat*at(1,3), alat*at(2,3), alat*at(3,3)
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
        CALL iotk_write_attr (attr,"condition","1")
        CALL iotk_write_begin(iun, "attrib",ATTR=attr)
        write(iun,100) tau
        CALL iotk_write_end(iun, "attrib")
     CALL iotk_write_end(iun, "particleset")
     ! </particleset>
  
     ! <particleset name="e">
     CALL iotk_write_attr (attr,"name","e",first=.true.)
     CALL iotk_write_attr (attr,"random","yes")
     CALL iotk_write_begin(iun, "particleset",ATTR=attr)
  
        ! <group name="u" size="" >
        CALL iotk_write_attr (attr,"name","u",first=.true.)
        CALL iotk_write_attr (attr,"size",nbndup)
        CALL iotk_write_begin(iun, "group",ATTR=attr)
        CALL iotk_write_attr (attr,"name","charge",first=.true.)
        CALL iotk_write_begin(iun, "parameter",ATTR=attr)
        write(iun,*) -1
        CALL iotk_write_end(iun, "parameter")
        CALL iotk_write_end(iun, "group")
  
        ! <group name="d" size="" >
        CALL iotk_write_attr (attr,"name","d",first=.true.)
        CALL iotk_write_attr (attr,"size",nbndup)
        CALL iotk_write_begin(iun, "group",ATTR=attr)
        CALL iotk_write_attr (attr,"name","charge",first=.true.)
        CALL iotk_write_begin(iun, "parameter",ATTR=attr)
        write(iun,*) -1
        CALL iotk_write_end(iun, "parameter")
        CALL iotk_write_end(iun, "group")
     CALL iotk_write_end(iun, "particleset")
     ! </particleset>
  
     ! <wavefunction name="psi0">
     CALL iotk_write_attr (attr,"name","psi0",first=.true.)
     CALL iotk_write_attr (attr,"type","PW")
     CALL iotk_write_attr (attr,"target","e")
     CALL iotk_write_begin(iun, "wavefunction",ATTR=attr)
       CALL iotk_write_attr (attr,"version","0.10",first=.true.)
       CALL iotk_write_begin(iun, "determinantset",ATTR=attr)
          CALL iotk_write_attr (attr,"ecut",ecutwfc/2,first=.true.)
          CALL iotk_write_begin(iun, "basisset",ATTR=attr)
          CALL iotk_write_end(iun, "basisset")
          
          CALL iotk_write_attr (attr,"href",TRIM(h5name),first=.true.)
          CALL iotk_write_begin(iun, "coefficients",ATTR=attr)
          CALL iotk_write_end(iun, "coefficients")
  
          ! write the index of the twist angle
          CALL iotk_write_attr (attr,"name","twistIndex",first=.true.)
          CALL iotk_write_begin(iun, "h5tag",ATTR=attr)
          write(iun,*) ik
          CALL iotk_write_end(iun, "h5tag")

          CALL iotk_write_attr (attr,"name","twistAngle",first=.true.)
          CALL iotk_write_begin(iun, "h5tag",ATTR=attr)
          write(iun,100) xk(1,ik+1),xk(2,ik+1),xk(3,ik+1)
          CALL iotk_write_end(iun, "h5tag")
  
          CALL iotk_write_begin(iun, "slaterdeterminant")
             ! build determinant for up electrons
             CALL iotk_write_attr (attr,"id","updet",first=.true.)
             CALL iotk_write_attr (attr,"size",nbndup)
             CALL iotk_write_begin(iun, "determinant",ATTR=attr)
                CALL iotk_write_attr (attr,"mode","ground",first=.true.)
                CALL iotk_write_attr (attr,"spindataset",0)
                CALL iotk_write_begin(iun, "occupation",ATTR=attr)
                CALL iotk_write_end(iun, "occupation")
             CALL iotk_write_end(iun, "determinant")
  
             ! build determinant for down electrons
             CALL iotk_write_attr (attr,"id","downdet",first=.true.)
             CALL iotk_write_attr (attr,"size",nbndup)
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

  ! generate hdf5 file containing all the parameters
  ! print out 2*occupied bands
  CALL pwhdf_open_file(tmp,h5len)
  CALL pwhdf_write_basis(g,igtog,ngtot)
  CALL pwhdf_write_parameters(nelec_tot,nspin,2*nbnd,nk,ecutwfc/2,alat,at)
  !
100 FORMAT (3(1x,f20.15))

  ALLOCATE (eigpacked(ngtot))
  
  ! start a main section to save eigen values and vector
  CALL pwhdf_open_eigg
  DO ik = 1, nk
     CALL pwhdf_open_twist(ik,xk(1,ik),2*nbnd,nspin)
     DO ispin = 1, nspin 
        ikk = ik + nk*(ispin-1)
        IF( nks > 1 ) THEN
           CALL gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
           CALL davcio(evc,nwordwfc,iunwfc,ikk,-1)
        ENDIF
        DO ibnd = 1, 2*nbnd
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

  !ALLOCATE (aux_ext(nr1s+1,nr2s+1,nr3s+1))

  ! open real-space wavefunction on FFT grid
  CALL pwhdf_open_eigr(nr1s,nr2s,nr3s)
  DO ik = 1, nk
     CALL pwhdf_open_twist(ik,xk(1,ik),2*nbnd,nspin)
     DO ispin = 1, nspin 
        ikk = ik + nk*(ispin-1)
        IF( nks > 1 ) THEN
           CALL gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
           CALL davcio(evc,nwordwfc,iunwfc,ikk,-1)
        ENDIF
        DO ibnd = 1, 2*nbnd
           !!transform G to R
           psic(:)=(0.d0,0.d0)
           psic(nls(igk(1:npw)))=evc(1:npw,ibnd)
           call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)

           !! is this correct????
           !aux_ext(1:nr1s,  1:nr2s,1:nr3s)=psic(1:nr1s,1:nr2s,1:nr3s)
           !aux_ext(knr1s+1, 1:nr2s,1:nr3s)=psic(1,     1:nr2s,1:nr3s)
           !aux_ext(1:nr1s,  nr2s+1,1:nr3s)=psic(1:nr1s,1,     1:nr3s)
           !aux_ext(1:nr1s,  1:nr2s, nr3s+1)=psic(1:nr1s,1:nr2s,1)
           !aux_ext(nr1s+1,nr2s+1,1:nr3s)=psic(1,     1,     1:nr3s)
           !aux_ext(nr1s+1,1:nr2s,nr3s+1)=psic(1,     1:nr2s,1)
           !aux_ext(1:nr1s,nr2s+1,nr3s+1)=psic(1:nr1s,1,     1)

           CALL pwhdf_write_wfr(ibnd,ispin,et(ibnd,ikk)/2,psic)
        ENDDO
     ENDDO
     CALL pwhdf_close_twist
  ENDDO
  CALL pwhdf_close_eigr

  ! ignore spin for the time being
  CALL pwhdf_write_rho(rho,rhog(1,1),ngm)

  ! close the file
  CALL pwhdf_close_file

  DEALLOCATE (igtog)
  DEALLOCATE (index)
  DEALLOCATE (becp)
  DEALLOCATE (aux)
  DEALLOCATE (hpsi)
  DEALLOCATE (eigpacked)
  !DEALLOCATE (aux_ext)

END SUBROUTINE compute_qmcpack


