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
  ! containing the  plane wave coefficients and other stuff needed by QMCPACK. 

  USE io_files,  ONLY : nd_nmbr, prefix, tmp_dir
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  USE mp_global,  ONLY : mp_startup, npool, nimage, nproc_pool, nproc_file, nproc_pool_file
  USE control_flags, ONLY : gamma_only
  USE mp_world,   ONLY : world_comm, nproc
  USE environment,ONLY : environment_start, environment_end
  USE KINDS, ONLY : DP
  !
  IMPLICIT NONE
  INTEGER :: ios
  LOGICAL :: write_psir, expand_kp, debug
  REAL(DP) :: t1, t2, dt
  ! directory for temporary files
  CHARACTER(len=256) :: outdir
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck

  NAMELIST / inputpp / prefix, outdir, write_psir, expand_kp, debug
#ifdef __MPI
  CALL mp_startup ( )
#endif

  CALL environment_start ( 'pw2qmcpack' )
#if defined(__HDF5) || defined(__HDF5_C)
  IF ( nimage > 1) &
     CALL errore('pw2qmcpack', ' image parallelization not (yet) implemented', 1)

  !   CALL start_postproc(nd_nmbr)
  ! 
  !   set default values for variables in namelist 
  ! 
  prefix = 'pwscf'
  write_psir = .false.
  expand_kp = .false.
  debug = .false.
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  ios = 0
  IF ( ionode )  THEN 
     !
     CALL input_from_file ( )
     !READ (5, inputpp, err=200, iostat=ios)
     READ (5, inputpp, iostat=ios)
     tmp_dir = trimcheck (outdir)
     !
  END IF
  CALL mp_bcast( ios, ionode_id, world_comm ) 
  IF ( ios/=0 ) CALL errore('pw2qmcpack', 'reading inputpp namelist', ABS(ios))
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast(prefix, ionode_id, world_comm ) 
  CALL mp_bcast(tmp_dir, ionode_id, world_comm ) 
  CALL mp_bcast(write_psir, ionode_id, world_comm ) 
  CALL mp_bcast(expand_kp, ionode_id, world_comm ) 
  CALL mp_bcast(debug, ionode_id, world_comm ) 
  !
  CALL start_clock ( 'read_file' )
  CALL read_file
  CALL stop_clock ( 'read_file' )
  IF ( gamma_only ) &
    CALL errore('pw2qmcpack', 'Using gamma trick results a reduced G space that is not supported by QMCPACK &
                             & though pw2qmcpack itself still can convert the WF to an h5 file in this case (experts only). &
                             & Please run pw.x with k point 0.0 0.0 0.0 instead of gamma.',1)
  !
  CALL openfil_pp
  !
  CALL start_clock ( 'compute_qmcpack' )
  CALL compute_qmcpack(write_psir, expand_kp, debug)
  CALL stop_clock ( 'compute_qmcpack' )
  !
  IF ( ionode ) THEN 
    WRITE( 6, * )
    !
    CALL print_clock( 'read_file_lite' )
    CALL print_clock( 'compute_qmcpack' )
    !
    WRITE( 6, '(/5x,"Called by read_file_lite:")' )
    CALL print_clock ( 'read_pseudo' )
    CALL print_clock ( 'read_rho' )
    CALL print_clock ( 'fft_rho' )
    CALL print_clock ( 'read_wave' )
    !
    WRITE( 6, '(/5x,"Called by compute_qmcpack:")' )
    CALL print_clock ( 'big_loop' )
    CALL print_clock ( 'write_h5' )
    CALL print_clock ( 'glue_h5' )
  ENDIF
#else
#error HDF5 flag neither enabled during configure nor added manually in make.inc
  CALL errore('pw2qmcpack', ' HDF5 flag not enabled during configure',1)
#endif
  CALL environment_end ( 'pw2qmcpack' )
  CALL stop_pp
  STOP

  

END PROGRAM pw2qmcpack

SUBROUTINE check_norm(ng, eigvec, collect_intra_pool, jks, ibnd, tag)
  USE mp_pools, ONLY: me_pool, intra_pool_comm
  USE mp,       ONLY: mp_sum
  USE KINDS,    ONLY: DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: eigvec(ng)
  INTEGER, INTENT(IN) :: ng, jks, ibnd
  LOGICAL, INTENT(IN) :: collect_intra_pool
  CHARACTER(len=*), INTENT(IN) :: tag
  !
  INTEGER :: ig
  REAL(DP) :: total_norm
  !
  ! check normalization before writing h5
  total_norm = 0.d0
  !$omp parallel do reduction(+:total_norm)
  DO ig=1, ng
    total_norm = total_norm + eigvec(ig)*CONJG(eigvec(ig))
  ENDDO
  ! collect within a pool
  IF(collect_intra_pool) CALL mp_sum ( total_norm , intra_pool_comm )
  ! compute the norm
  total_norm = SQRT(total_norm)
  ! check abort if necessary
  IF (me_pool==0 .AND. ABS(total_norm-1.d0)>1.d-6) THEN
    WRITE(*,"(A,I3,A,I3,3A,F20.15)") "The wrong norm of k-point ", jks, " band ", ibnd, " , ", &
                                      tag, ", is ", total_norm
    IF (.NOT. collect_intra_pool) THEN
      WRITE(*,"(3A)") "The orbitals went wrong before being written to h5 file. ", &
                      "Please first add debug=.true. in the pw2qmcpack input file to check ", &
                      "if the orbitals can be read from QE files correctly."
    ELSE
      WRITE(*,"(2A)") "The orbitals read from QE files are incorrect. ", &
                      "Perhaps your QE orbital files are corrupted."
    ENDIF
    STOP
  ENDIF
  !
END SUBROUTINE

SUBROUTINE compute_qmcpack(write_psir, expand_kp, debug)

  USE kinds, ONLY: DP
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
  USE cell_base, ONLY: omega, alat, tpiba2, at, bg
  USE constants, ONLY: tpi
  USE run_info,  ONLY: title
  USE gvect, ONLY: ngm, ngm_g, g, ig_l2g
  USE klist , ONLY: nks, nelec, nelup, neldw, wk, xk, nkstot
  USE lsda_mod, ONLY: lsda, nspin, isk
  USE scf, ONLY: rho, rho_core, rhog_core, vnew
  USE wvfct, ONLY: npw, npwx, nbnd, g2kin, wg, et
  USE klist, ONLY: igk_k
  USE gvecw, ONLY : ecutwfc
  USE control_flags, ONLY: gamma_only
  USE becmod, ONLY: becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_global, ONLY: stdout, ionode,  ionode_id
  USE mp_world, ONLY: world_comm, mpime
  USE io_files, ONLY: nd_nmbr, nwordwfc, iunwfc, iun => iunsat, tmp_dir, prefix
  USE wavefunctions, ONLY : evc, psic
  use iotk_module
  use iotk_xtox_interf
  USE mp_global,            ONLY: inter_pool_comm, intra_pool_comm, nproc_pool, kunit
  USE mp_global,            ONLY: npool, my_pool_id, intra_image_comm
  USE mp_pools,             ONLY: me_pool
  USE mp,                   ONLY: mp_sum, mp_max, mp_bcast, mp_barrier
  use scatter_mod,          ONLY : gather_grid, scatter_grid 
  use fft_base,             ONLY : dffts
  use fft_interfaces,       ONLY : invfft, fwfft
  USE dfunct, ONLY : newd
  USE symm_base,            ONLY : nsym, s, ftau

  IMPLICIT NONE
  LOGICAL :: write_psir, expand_kp, debug
  LOGICAL :: pool_ionode
  INTEGER :: ig, ibnd, ik, io, na, j, ispin, nbndup, nbnddown, &
       nk, ngtot, ig7, ikk, iks, kpcnt, jks, nt, ijkb0, ikb, ih, jh, jkb, at_num, &
       nelec_tot, nelec_up, nelec_down, ii, igx, igy, igz, n_rgrid(3), &
       nkqs, nr1s,nr2s,nr3s
  INTEGER, ALLOCATABLE :: indx(:), igtog(:), igtomin(:), g_global_to_local(:)
  LOGICAL :: exst, found
  REAL(DP) :: ek, eloc, enl, charge, etotefield
  REAL(DP) :: bg_qmc(3,3), g_red(3), lattice_real(3,3)
  COMPLEX(DP), ALLOCATABLE :: phase(:),eigpacked(:)
  COMPLEX(DP), ALLOCATABLE :: psitr(:)
  REAL(DP), ALLOCATABLE ::  tau_r(:,:),psireal(:),eigval(:) !,g_cart(:,:)
  INTEGER :: ios, ierr, h5len,oldh5,ig_c,save_complex, nup,ndown
  INTEGER, EXTERNAL :: atomic_number, is_complex
  !REAL(DP), ALLOCATABLE :: g_qmc(:,:)
  INTEGER, ALLOCATABLE :: gint_den(:,:), gint_qmc(:,:)
  COMPLEX(DP), ALLOCATABLE :: den_g_global(:,:)
  REAL (DP), EXTERNAL :: ewald
  COMPLEX(DP), ALLOCATABLE, TARGET :: tmp_psic(:)
  COMPLEX(DP), DIMENSION(:), POINTER :: psiCptr
  REAL(DP), DIMENSION(:), POINTER :: psiRptr
! **********************************************************************
  INTEGER :: npw_sym  
  INTEGER, ALLOCATABLE, TARGET :: igk_sym(:)
  REAL(DP), ALLOCATABLE :: g2kin_sym(:)
! **********************************************************************
  INTEGER :: nkfull,max_nk,max_sym,isym,nxxs
  INTEGER , ALLOCATABLE :: num_irrep(:) 
  INTEGER, ALLOCATABLE :: xkfull_index(:,:) ! maps to sym_list and xk_full_list  
  INTEGER, ALLOCATABLE :: sym_list(:)
  REAL(DP),    ALLOCATABLE :: xk_full_list(:,:)
  REAL(DP) :: t1, t2, dt
  integer, allocatable :: rir(:)  
  COMPLEX(DP), ALLOCATABLE :: tmp_evc(:)

  CHARACTER(256)          :: tmp,h5name,eigname,tmp_combo
  CHARACTER(iotk_attlenx) :: attr
  
  INTEGER :: rest, nbase, basekindex, nktot
  REAL(DP) :: xk_cryst(3)

  INTEGER :: npwx_tot, igk_g;
  
  NULLIFY(psiRptr)
  NULLIFY(psiCptr)

  ! Ye Luo
  ! define the pool level ionode
  ! an image ionode must be pool ionode
  if(me_pool==0) then 
    pool_ionode=.true.
  else
    pool_ionode=.false.
  endif
  
  ! MAMorales:
  ! removed USPP functions

  ! Ye Luo:
  ! sum up npwx to npwx_tot inside a pool and maximize it among pools.
  npwx_tot = npwx
  CALL mp_sum ( npwx_tot, intra_pool_comm )
  CALL mp_max ( npwx_tot, inter_pool_comm )
  !write(*,*) mpime, ionode_id, npwx_tot, npw
  
  ! this limits independent definition of ecutrho to < 4*ecutwf
  ! four times npwx should be enough
  ALLOCATE (indx (4*npwx_tot) )
  ALLOCATE (igtog (4*npwx_tot) )
  ALLOCATE (g_global_to_local(ngm_g) )
  ALLOCATE (igtomin(4*npwx_tot) )
  ALLOCATE (tmp_evc(npwx_tot) )

  indx(:) = 0
  igtog(:) = 0
  igtomin(:) = 0

  rest = ( nkstot - kunit * ( nkstot / kunit / npool ) * npool ) / kunit
  nbase = nks * my_pool_id
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest * kunit
  !write(*,*) "debug",mpime, nks, nbase
  
  IF( lsda ) THEN
!      IF( expand_kp ) &
!        CALL errore ('pw2qmcpack','expand_kp not implemented with nspin>1`', 1)     
     nbndup = nbnd
     nbnddown = nbnd
     nk = nks/2
     nktot = nkstot/2
     !     nspin = 2
  ELSE
     nbndup = nbnd
     nbnddown = 0
     nk = nks
     nktot = nkstot
     !     nspin = 1
  ENDIF
  
! !    sanity check for lsda logic to follow 
!   if (ionode) then
!     DO ik = 1, nktot
!       iks=ik+nktot
!       xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik) - ( at(1,:)*xk(1,iks) + at(2,:)*xk(2,iks) + at(3,:)*xk(3,iks))
!       if (abs(xk_cryst(1))+abs(xk_cryst(2))+abs(xk_cryst(3)) .gt. 1e-12) then
!         print *,"not paired %i %i",ik,iks
!       endif
!     ENDDO
!   endif
   
   
  !
  
  ! for now, I'm assuming that symmetry rotations do not affect npw, 
  ! meaning that rotations don't displace elements outside the cutoff 
  nr1s = dffts%nr1
  nr2s = dffts%nr2
  nr3s = dffts%nr3
  nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x
  allocate (igk_sym( npwx ), g2kin_sym ( npwx ) )

  if (ionode) then
    if(expand_kp) then
      max_sym = min(48, 2 * nsym)
      max_nk = nktot * max_sym 
      ALLOCATE(num_irrep(nktot),xkfull_index(nktot,max_sym),sym_list(max_nk))
      ALLOCATE(xk_full_list(3,max_nk))
      ALLOCATE(rir(nxxs))
      call generate_symmetry_equivalent_list() 
      if(ionode) print *,'Total number of k-points after expansion:',nkfull
    else
      ALLOCATE(num_irrep(nktot),xkfull_index(nktot,1),sym_list(nktot))
      ALLOCATE(xk_full_list(3,nktot))
      nkfull = nktot
      do ik = 1, nktot
        xk_full_list(:,ik) = xk(:,ik) 
        num_irrep(ik) = 1
        sym_list(ik) = 1 
        xkfull_index(ik,1) = ik  
      enddo 
    endif
  else
    if(expand_kp) then
      max_sym = min(48, 2 * nsym)
      max_nk = nktot * max_sym 
      ALLOCATE(num_irrep(nktot),xkfull_index(nktot,max_sym),sym_list(max_nk))
      ALLOCATE(xk_full_list(3,max_nk))
      ALLOCATE(rir(nxxs))
    else
      ALLOCATE(num_irrep(nktot),xkfull_index(nktot,1),sym_list(nktot))
      ALLOCATE(xk_full_list(3,nktot))
      nkfull = nktot
    endif
  endif

  CALL mp_bcast(xkfull_index, ionode_id, world_comm ) 
  CALL mp_bcast(xk_full_list, ionode_id, world_comm ) 
  CALL mp_bcast(sym_list, ionode_id, world_comm ) 
  CALL mp_bcast(num_irrep, ionode_id, world_comm ) 
  CALL mp_bcast(nkfull, ionode_id, world_comm ) 
  
  ! IF ( nbase > 0 ) THEN
  !    num_irrep(1:nks) = num_irrep(nbase+1:nbase+nks)
  !    xk_full_list(:,1:nks) = xk_full_list(:,nbase+1:nbase+nks)
  ! END IF  

  DO ik = 1, nks
     basekindex = ik + nbase
     CALL gk_sort (xk (1, basekindex), ngm, g, ecutwfc / tpiba2, npw, igk_k(:,ik), g2kin)

     DO ig =1, npw
        ! mapping to the global g vectors.
        igk_g = ig_l2g(igk_k(ig,ik))
        IF( igk_g > 4*npwx_tot ) & 
             CALL errore ('pw2qmcpack','increase allocation of index', ig)
        indx( igk_g ) = 1
     ENDDO
  ENDDO
  call mp_max(indx, world_comm)

  ngtot = 0
  ! igtomin maps indices from the full set of G-vectors to the
  ! minimal set which includes the G-spheres of all k-points
  DO ig = 1, 4*npwx_tot
    IF( indx(ig) > 0 ) THEN
        ngtot = ngtot + 1
        igtog(ngtot) = ig
        igtomin(ig) = ngtot
    ENDIF
  ENDDO

  if (debug) then
    write(6,"(A)")  "     pw2qmcpack found"
    write(6,"(A,A8,2A8,4A10,1A4)") "     MPI rank", "pool id", "npwx", "npw", "npwx_tot", "ngtot", "ngm", "ngm_g", "nks"
    write(*,"(A,I9,I8,2I8,4I10,1I4)") "    ", mpime, me_pool, npwx, npw, npwx_tot, ngtot, ngm, ngm_g, nks
  endif

  ALLOCATE (gint_qmc(3,ngtot))
  ALLOCATE (gint_den(3,ngm_g))
  ALLOCATE (den_g_global(ngm_g,nspin))
  !ALLOCATE (g_qmc(3,ngtot))
  !ALLOCATE (g_cart(3,ngtot))
  ALLOCATE (tau_r(3,nat))

  ! get the number of electrons
  nelec_tot= NINT(nelec)
  nup=NINT(nelup)
  ndown=NINT(neldw)

  if(nup .eq. 0) then
    ndown=nelec_tot/2
    nup=nelec_tot-ndown
  endif

  bg_qmc(:,:)=bg(:,:)/alat

  if((npool>1) .and. (my_pool_id>0)) then
    h5name = TRIM( prefix ) // '.pwscf.h5' // "_part"//trim(iotk_itoa(my_pool_id))
  else
    h5name = TRIM( prefix ) // '.pwscf.h5'
  endif
  eigname = "eigenstates_"//trim(iotk_itoa(nr1s))//'_'//trim(iotk_itoa(nr2s))//'_'//trim(iotk_itoa(nr3s))

  tmp = TRIM( tmp_dir )//TRIM( h5name ) 
  h5len = LEN_TRIM(tmp)
  
#if defined(__HDF5) || defined(__HDF5_C)
  ! writing to xml and hdf5
  ! open hdf5 file 
  oldh5=0
  if(pool_ionode) CALL esh5_open_file(tmp,h5len,oldh5)


  if(ionode) then
  !! create a file for particle set
  tmp = TRIM( tmp_dir ) // TRIM( prefix )// '.ptcl.xml'
  CALL iotk_open_write(iun, FILE=TRIM(tmp), ROOT="qmcsystem", IERR=ierr )

  CALL iotk_write_attr (attr,"name","global",first=.true.)
  CALL iotk_write_begin(iun, "simulationcell",ATTR=attr)
  CALL iotk_write_attr (attr,"name","lattice",first=.true.)
  CALL iotk_write_attr (attr,"units","bohr")
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
   
  lattice_real=alat*at
  WRITE(iun,100) lattice_real(1,1), lattice_real(2,1), lattice_real(3,1)
  WRITE(iun,100) lattice_real(1,2), lattice_real(2,2), lattice_real(3,2)
  WRITE(iun,100) lattice_real(1,3), lattice_real(2,3), lattice_real(3,3)

  CALL esh5_write_supercell(lattice_real)

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
  WRITE(iun,'(a)') '15'
  CALL iotk_write_end(iun, "parameter")
  CALL iotk_write_end(iun, "simulationcell")

  ! <particleset name="ions">
  CALL iotk_write_attr (attr,"name","ion0",first=.true.)
  CALL iotk_write_attr (attr,"size",nat)
  CALL iotk_write_begin(iun, "particleset",ATTR=attr)

  CALL esh5_open_atoms(nat,ntyp)

  ! ionic species --> group
  DO na=1,ntyp

  tmp=TRIM(atm(na))
  h5len=LEN_TRIM(tmp)
  CALL esh5_write_species(na,tmp,h5len,atomic_number(tmp),zv(na))

  CALL iotk_write_attr (attr,"name",TRIM(atm(na)),first=.true.)
  CALL iotk_write_begin(iun, "group",ATTR=attr)
  ! charge
  CALL iotk_write_attr (attr,"name","charge",first=.true.)
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  write(iun,"(8X, F3.1)") zv(na)
  CALL iotk_write_end(iun, "parameter")
  ! atomic number
  CALL iotk_write_attr (attr,"name","atomicnumber",first=.true.)
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  write(iun,"(8X, I3)") atomic_number(TRIM(atm(na)))
  CALL iotk_write_end(iun, "parameter")
  !
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
  tau_r(1,na)=alat*tau(1,na)
  tau_r(2,na)=alat*tau(2,na)
  tau_r(3,na)=alat*tau(3,na)
  WRITE(iun,100) (tau_r(j,na),j=1,3)
  ENDDO
  !write(iun,100) tau
  CALL iotk_write_end(iun, "attrib")
  CALL iotk_write_end(iun, "particleset")

  !cartesian positions
  CALL esh5_write_positions(tau_r)
  CALL esh5_write_species_ids(ityp)

  CALL esh5_close_atoms()
  ! </particleset>

  ! <particleset name="e">
  CALL iotk_write_attr (attr,"name","e",first=.true.)
  CALL iotk_write_attr (attr,"random","yes")
  CALL iotk_write_attr (attr,"random_source","ion0")
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
  !!DO ik = 0, nk-1
  ik=0
   ! NOT create a xml input file for each k-point
   !  IF(nk .gt. 1) THEN
   !    tmp = TRIM( tmp_dir ) // TRIM( prefix ) //TRIM(iotk_index(ik))// '.wfs.xml'
   !  ELSE
   !    tmp = TRIM( tmp_dir ) // TRIM( prefix )// '.wfs.xml'
   !  ENDIF
   tmp = TRIM( tmp_dir ) // TRIM( prefix )// '.wfs.xml'
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
       CALL iotk_write_attr (attr,"sort","1")
       CALL iotk_write_attr (attr,"tilematrix","1 0 0 0 1 0 0 0 1")
       CALL iotk_write_attr (attr,"twistnum","0")
       CALL iotk_write_attr (attr,"source","ion0")
       CALL iotk_write_attr (attr,"version","0.10")
       CALL iotk_write_begin(iun, "determinantset",ATTR=attr)
          CALL iotk_write_attr (attr,"ecut",ecutwfc/2,first=.true.)
          ! basisset to overwrite cutoff to a smaller value
          !CALL iotk_write_begin(iun, "basisset",ATTR=attr)
          !   ! add grid to use spline on FFT grid
          !   CALL iotk_write_attr (attr,"dir","0",first=.true.)
          !   CALL iotk_write_attr (attr,"npts",nr1s)
          !   CALL iotk_write_attr (attr,"closed","no")
          !   CALL iotk_write_empty(iun, "grid",ATTR=attr)
          !   CALL iotk_write_attr (attr,"dir","1",first=.true.)
          !   CALL iotk_write_attr (attr,"npts",nr2s)
          !   CALL iotk_write_attr (attr,"closed","no")
          !   CALL iotk_write_empty(iun, "grid",ATTR=attr)
          !   CALL iotk_write_attr (attr,"dir","2",first=.true.)
          !   CALL iotk_write_attr (attr,"npts",nr3s)
          !   CALL iotk_write_attr (attr,"closed","no")
          !   CALL iotk_write_empty(iun, "grid",ATTR=attr)
          !CALL iotk_write_end(iun, "basisset")
          
          !CALL iotk_write_attr (attr,"href",TRIM(h5name),first=.true.)
          !CALL iotk_write_empty(iun, "coefficients",ATTR=attr)
  
          ! write the index of the twist angle
          !!!! remove twistIndex and twistAngle
          !using determinantset@twistnum
          !CALL iotk_write_attr (attr,"name","twistIndex",first=.true.)
          !CALL iotk_write_begin(iun, "h5tag",ATTR=attr)
          !write(iun,*) ik
          !CALL iotk_write_end(iun, "h5tag")

          !CALL iotk_write_attr (attr,"name","twistAngle",first=.true.)
          !CALL iotk_write_begin(iun, "h5tag",ATTR=attr)
          !g_red(1)=at(1,1)*xk(1,ik+1)+at(2,1)*xk(2,ik+1)+at(3,1)*xk(3,ik+1)
          !g_red(2)=at(1,2)*xk(1,ik+1)+at(2,2)*xk(2,ik+1)+at(3,2)*xk(3,ik+1)
          !g_red(3)=at(1,3)*xk(1,ik+1)+at(2,3)*xk(2,ik+1)+at(3,3)*xk(3,ik+1)
          !!write(iun,100) xk(1,ik+1),xk(2,ik+1),xk(3,ik+1)
          !write(iun,100) g_red(1),g_red(2),g_red(3)
          !CALL iotk_write_end(iun, "h5tag")
          !write(iun,'(a)') '<!-- Uncomment this out for bspline wavefunctions '
          !!CALL iotk_write_attr (attr,"name","eigenstates",first=.true.)
          !!CALL iotk_write_begin(iun, "h5tag",ATTR=attr)
          !!write(iun,'(a)') TRIM(eigname)
          !!CALL iotk_write_end(iun, "h5tag")
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

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! two-body jastro
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL iotk_write_attr (attr,"name","J2",first=.true.)
       CALL iotk_write_attr (attr,"type","Two-Body");
       CALL iotk_write_attr (attr,"function","Bspline");
       CALL iotk_write_attr (attr,"print","yes");
       CALL iotk_write_begin(iun, "jastrow",ATTR=attr)

         ! for uu
         CALL iotk_write_attr (attr,"speciesA","u",first=.true.)
         CALL iotk_write_attr (attr,"speciesB","u")
         !CALL iotk_write_attr (attr,"rcut","10")
         CALL iotk_write_attr (attr,"size","8")
         CALL iotk_write_begin(iun, "correlation",ATTR=attr)
           CALL iotk_write_attr (attr,"id","uu",first=.true.)
           CALL iotk_write_attr (attr,"type","Array")
           CALL iotk_write_begin(iun, "coefficients",ATTR=attr)
           write(iun,*) "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
           CALL iotk_write_end(iun, "coefficients")
         CALL iotk_write_end(iun, "correlation")

         ! for ud
         CALL iotk_write_attr (attr,"speciesA","u",first=.true.)
         CALL iotk_write_attr (attr,"speciesB","d")
         !CALL iotk_write_attr (attr,"rcut","10")
         CALL iotk_write_attr (attr,"size","8")
         CALL iotk_write_begin(iun, "correlation",ATTR=attr)
           CALL iotk_write_attr (attr,"id","ud",first=.true.)
           CALL iotk_write_attr (attr,"type","Array")
           CALL iotk_write_begin(iun, "coefficients",ATTR=attr)
           write(iun,*) "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
           CALL iotk_write_end(iun, "coefficients")
         CALL iotk_write_end(iun, "correlation")

       CALL iotk_write_end(iun, "jastrow")

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! one-body jastro
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL iotk_write_attr (attr,"name","J1",first=.true.)
       CALL iotk_write_attr (attr,"type","One-Body");
       CALL iotk_write_attr (attr,"function","Bspline");
       CALL iotk_write_attr (attr,"source","ion0");
       CALL iotk_write_attr (attr,"print","yes");
       CALL iotk_write_begin(iun, "jastrow",ATTR=attr)

       DO na=1,ntyp
         tmp=TRIM(atm(na))
         tmp_combo='e'//TRIM(atm(na))

         !h5len=LEN_TRIM(tmp)
         CALL iotk_write_attr (attr,"elementType",TRIM(tmp),first=.true.)
         !CALL iotk_write_attr (attr,"rcut","10")
         CALL iotk_write_attr (attr,"size","8")
         CALL iotk_write_begin(iun, "correlation",ATTR=attr)

         CALL iotk_write_attr (attr,"id",TRIM(tmp_combo),first=.true.)
         CALL iotk_write_attr (attr,"type","Array")
         CALL iotk_write_begin(iun, "coefficients",ATTR=attr)
         write(iun,*) "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
         CALL iotk_write_end(iun, "coefficients")
         CALL iotk_write_end(iun, "correlation")
       ENDDO
       CALL iotk_write_end(iun, "jastrow")
     
     CALL iotk_write_end(iun, "wavefunction")
  
     CALL iotk_close_write(iun)
  !ENDDO

  endif ! ionode

  g_global_to_local(:)=0
  DO ig=1,ngm
    g_global_to_local(ig_l2g(ig))=ig
  ENDDO

  gint_qmc(:,:) = 0
  DO ig=1, ngtot
    ig_c = g_global_to_local(igtog(ig))
    if (ig_c>0) then
        !g_cart(1,ig)=tpi/alat*g(1,ig_c)
        !g_cart(2,ig)=tpi/alat*g(2,ig_c)
        !g_cart(3,ig)=tpi/alat*g(3,ig_c)
        !g_qmc(1,ig)=at(1,1)*g(1,ig_c)+at(2,1)*g(2,ig_c)+at(3,1)*g(3,ig_c)
        !g_qmc(2,ig)=at(1,2)*g(1,ig_c)+at(2,2)*g(2,ig_c)+at(3,2)*g(3,ig_c)
        !g_qmc(3,ig)=at(1,3)*g(1,ig_c)+at(2,3)*g(2,ig_c)+at(3,3)*g(3,ig_c)
        gint_qmc(1,ig)=NINT(at(1,1)*g(1,ig_c)+at(2,1)*g(2,ig_c)+at(3,1)*g(3,ig_c))
        gint_qmc(2,ig)=NINT(at(1,2)*g(1,ig_c)+at(2,2)*g(2,ig_c)+at(3,2)*g(3,ig_c))
        gint_qmc(3,ig)=NINT(at(1,3)*g(1,ig_c)+at(2,3)*g(2,ig_c)+at(3,3)*g(3,ig_c))
        !WRITE(io,'(3(1x,f20.15))') g_cart(1,ig),g_cart(2,ig),g_cart(3,ig)
    endif
  ENDDO
  call mp_sum(gint_qmc, intra_pool_comm)

  gint_den(:,:) = 0
  DO ig=1,ngm
    gint_den(1,ig_l2g(ig))=NINT(at(1,1)*g(1,ig)+at(2,1)*g(2,ig)+at(3,1)*g(3,ig))
    gint_den(2,ig_l2g(ig))=NINT(at(1,2)*g(1,ig)+at(2,2)*g(2,ig)+at(3,2)*g(3,ig))
    gint_den(3,ig_l2g(ig))=NINT(at(1,3)*g(1,ig)+at(2,3)*g(2,ig)+at(3,3)*g(3,ig))
  ENDDO
  DO ispin=1,nspin
    den_g_global(:,ispin) = 0.d0
    DO ig=1,ngm
      den_g_global(ig_l2g(ig),ispin) = rho%of_g(ig,ispin)
    ENDDO
  ENDDO

  call mp_sum(gint_den, intra_pool_comm)
  call mp_sum(den_g_global, intra_pool_comm)

  n_rgrid(1)=nr1s
  n_rgrid(2)=nr2s
  n_rgrid(3)=nr3s

  save_complex=0
  if(ionode) then
    DO ik = 1, nktot
      !! evaluate the phase
      !phase(:) = (0.d0,0.d0)
      !if ( ig_(ik,ib)>0) phase( dffts%nl(ig_(ik,ib)) ) = (1.d0,0.d0)
      g_red(1)=at(1,1)*xk_full_list(1,ik)+at(2,1)*xk_full_list(2,ik)+at(3,1)*xk_full_list(3,ik)
      g_red(2)=at(1,2)*xk_full_list(1,ik)+at(2,2)*xk_full_list(2,ik)+at(3,2)*xk_full_list(3,ik)
      g_red(3)=at(1,3)*xk_full_list(1,ik)+at(2,3)*xk_full_list(2,ik)+at(3,3)*xk_full_list(3,ik)

      IF(g_red(1)*g_red(1)+g_red(2)*g_red(2)+g_red(3)*g_red(3)>1e-12) THEN
	  save_complex=1
      END IF
    END DO
  endif
  
  CALL mp_bcast(save_complex, ionode_id, world_comm )
  

  
!     WRITE(io,'(A10,3(1x,i6))') 'ngrid: ',n_rgrid(1:3) 

  !CALL esh5_open_electrons(nup, ndown,nspin,nk,nbnd,n_rgrid)!, save_complex)
  !CALL esh5_open_electrons(nup, ndown, nspin, nkfull, nbnd, n_rgrid)!, save_complex)

  if(pool_ionode) then
    if(ionode) then
      CALL esh5_open_electrons_base(nup, ndown, nspin, nkfull, nbnd, n_rgrid)!, save_complex)
    else
      CALL esh5_open_electrons(nup, ndown, nspin, nkfull, nbnd, n_rgrid)!, save_complex)
    endif
  endif
  
!   IF (write_psir) THEN
!     CALL esh5_write_psi_r_mesh(n_rgrid)
!   ENDIF

  !!NOT YET DECIDED
  !!CALL esh5_write_basis(g_qmc,g_cart,ngtot)
  !!CALL esh5_write_parameters(nelec_tot,nspin,nbnd,nkfull,ecutwfc/2,alat,at)
  !

  ALLOCATE (eigpacked(ngtot))
  ALLOCATE (eigval(nbnd))

!ionode writes all k-point and ev data
  if(ionode)then
    DO ik = 1, nkstot
      basekindex = ik + nbase
      ispin = 1
      if (basekindex > nktot) then
        ispin = 2
	basekindex = basekindex - nktot
      endif
      DO iks = 1,num_irrep(basekindex)  
	jks = xkfull_index(basekindex,iks)
	g_red(1)=at(1,1)*xk_full_list(1,jks)+at(2,1)*xk_full_list(2,jks)+at(3,1)*xk_full_list(3,jks)
	g_red(2)=at(1,2)*xk_full_list(1,jks)+at(2,2)*xk_full_list(2,jks)+at(3,2)*xk_full_list(3,jks)
	g_red(3)=at(1,3)*xk_full_list(1,jks)+at(2,3)*xk_full_list(2,jks)+at(3,3)*xk_full_list(3,jks)        
	
        CALL esh5_open_kpoint(jks)
        CALL esh5_write_kpoint_data(g_red,wk(basekindex)/num_irrep(basekindex),ngtot,iks,num_irrep(basekindex))

  !     only the 1 index kpoint will write this g vectors
	if(ik == 1) then
	    CALL esh5_write_gvectors_k(gint_qmc,ngtot)
	endif

! 	if (lsda) then
! 	  ispin = isk(ik)
! 	else
! 	  ispin=1
! 	endif
	
        CALL esh5_open_spin(ispin)
	DO ibnd = 1, nbnd
	  eigval(ibnd)=0.5*et(ibnd,ik)
	ENDDO
	CALL esh5_write_eigvalues(eigval)
	CALL esh5_close_spin()
	

        CALL esh5_close_kpoint()
      ENDDO
    ENDDO
  else
    DO ik = 1, nks
      basekindex = ik + nbase
      if (basekindex > nktot) then
	basekindex = basekindex - nktot
	ispin=2
      else
	ispin=1
      endif
      DO iks = 1,num_irrep(basekindex)  
	jks = xkfull_index(basekindex,iks)
	g_red(1)=at(1,1)*xk_full_list(1,jks)+at(2,1)*xk_full_list(2,jks)+at(3,1)*xk_full_list(3,jks)
	g_red(2)=at(1,2)*xk_full_list(1,jks)+at(2,2)*xk_full_list(2,jks)+at(3,2)*xk_full_list(3,jks)
	g_red(3)=at(1,3)*xk_full_list(1,jks)+at(2,3)*xk_full_list(2,jks)+at(3,3)*xk_full_list(3,jks)

	!! open kpoint 
	if(pool_ionode) CALL esh5_open_kpoint(jks)
! 	CALL esh5_write_kpoint_data(g_red,wk(ik)/num_irrep(basekindex),ngtot)
! 	if (lsda) then
! 	  ispin = isk(ik)
! 	else
! 	  ispin=1
! 	endif
	if(pool_ionode) CALL esh5_open_spin(ispin)
	if(pool_ionode) CALL esh5_close_spin()
      
	if(pool_ionode) CALL esh5_close_kpoint()

      ENDDO
    ENDDO  
  endif

100 FORMAT (3(1x,f20.15))

  if(save_complex /=1 .and. write_psir) ALLOCATE(psireal(nxxs))
  if(write_psir .or. expand_kp) then
    ALLOCATE(psitr(nxxs))
    IF(nproc_pool > 1) ALLOCATE(tmp_psic(nxxs))
  endif

!   if(ionode) print *,'PW2QMCPACK npw=',npw,'ngtot=',ngtot
  ! open real-space wavefunction on FFT grid
  !!CALL esh5_open_eigr(nr1s,nr2s,nr3s)
  !DO ik = 1, nk
  
  CALL start_clock ( 'big_loop' )
  if(nks .eq. 1) then ! treat 1 kpoint specially
    if(pool_ionode) write(*,"(A,I8,A)") '     k pool ', my_pool_id, ' has only 1 Kpoint. Bypass everything '
    ik=1
    basekindex = ik + nbase
    if (basekindex > nktot) then
      basekindex = basekindex - nktot
      ispin=2
    else
      ispin=1
    endif
    if(debug) write(6,*) "     starting davcio!"
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
    if(debug) write(6,*) "     davcio ends!"
    if(pool_ionode) CALL esh5_open_kpoint(basekindex)
    if(pool_ionode) CALL esh5_open_spin(ispin)
    DO ibnd = 1, nbnd
      eigpacked(:)=(0.d0,0.d0)
      eigpacked(igtomin(ig_l2g(igk_k(1:npw,ik))))=evc(1:npw,ibnd)
      IF (debug) CALL check_norm(npw, evc(1,ibnd), .true., jks, ibnd, "before collection")
      if (debug) write(6,"(A,1I5)") "     collecting band ", ibnd
      CALL mp_sum ( eigpacked , intra_pool_comm )
      if (debug) write(6,"(A,1I5)") "        writing band ", ibnd
      if (pool_ionode) THEN
        CALL check_norm(ngtot, eigpacked, .false., jks, ibnd, "after collection before writing")
        CALL esh5_write_psi_g(ibnd,eigpacked,ngtot)
      ENDIF
    enddo
    if(pool_ionode) CALL esh5_close_spin()
    if(pool_ionode) CALL esh5_close_kpoint() 
  else ! nk .neq. 1
    DO ik = 1, nks
    basekindex = ik + nbase
    if (basekindex > nktot) then
      basekindex = basekindex - nktot
      ispin=2
    else
      ispin=1
    endif
    DO iks = 1,num_irrep(basekindex)  
     jks = xkfull_index(basekindex,iks)
     isym = sym_list(jks)

     if(expand_kp) then
        call generate_symmetry_rotation(isym)
     endif

     if(pool_ionode) CALL esh5_open_kpoint(jks)

!      if(ionode) print *,'PW2QMCPACK ik,iks=',ik,iks

!      DO ispin = 1, nspin 
!         ikk = ik + nk*(ispin-1)
!       if (lsda) then
!         ispin = isk(ik)
!       else
!         ispin=1
!       endif

        !!! MAM: This could be outside the num_irrep group is ispin = 1,
        !!!      can I switch the order of esh5_open_spin and
        !!!      esh5_open_kpoint??? 
        CALL gk_sort (xk (1:3, ik), ngm, g, ecutwfc / tpiba2, npw, igk_k(:,ik), g2kin)
        CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
        CALL gk_sort (xk_full_list (1:3, jks), ngm, g, ecutwfc / tpiba2, npw_sym, igk_sym, g2kin_sym)
        if(npw .ne. npw_sym )  then
          write(*,*) 'Warning!!!: npw != npw_sym: ',npw,npw_sym
        endif 

        if(pool_ionode) CALL esh5_open_spin(ispin)

        DO ibnd = 1, nbnd !!transform G to R

! I should be able to do the rotation directly in G space, 
! but for now I'm doing it like this
           IF(expand_kp) then
             psic(:)=(0.d0,0.d0)
             psitr(:)=(0.d0,0.d0)
             tmp_evc(:) = (0.d0,0.d0) 
             IF(nproc_pool > 1) THEN
               ! 
               psic(dffts%nl(ig_l2g(igk_k(1:npw,ik))))=evc(1:npw,ibnd)
               
!                call errore ('pw2qmcpack','parallel version not fully implemented.',2)
               if(gamma_only) then
                      call errore ('pw2qmcpack','problems with gamma_only, not fully implemented.',2)
               endif
               !
               CALL invfft ('Wave', psic, dffts)
               !
!               call cgather_smooth(psic,psitr)
               call gather_grid(dffts,psic,psitr)
               tmp_psic(1:nxxs) = psitr(rir(1:nxxs))
!               call cscatter_smooth(tmp_psic,psic)
               call scatter_grid(dffts,tmp_psic,psic)
               !
               ! at this point, I have the rotated orbital in real space, might
               ! want to keep it stored somewhere for later use if write_psir 
               ! 
               CALL fwfft ('Wave', psic, dffts)
               !
               tmp_evc(1:npw_sym)=psic(dffts%nl(ig_l2g(igk_sym(1:npw_sym))))
               ! 
             ELSE ! nproc_pool <= 1
               ! 
               psic(dffts%nl(ig_l2g(igk_k(1:npw,ik))))=evc(1:npw,ibnd)
               if(gamma_only) then
                      call errore ('pw2qmcpack','problems with gamma_only, not fully implemented.',2)
               endif
               !
               CALL invfft ('Wave', psic, dffts)
               !
               psitr(1:nxxs) = psic(rir(1:nxxs))
               ! temporary hack to see if I have problems with inversion
               ! symmetry
               if(isym.lt.0 .AND. iks.gt.1 .AND. abs(isym).eq.abs(sym_list(xkfull_index(basekindex,iks-1)))   ) then  
                 psitr(1:nxxs) = CONJG(psitr(1:nxxs)) 
               endif
               !psitr(:) = psic(:)
               !
               CALL fwfft ('Wave', psitr, dffts)
               !
               tmp_evc(1:npw_sym)=psitr(dffts%nl(ig_l2g(igk_sym(1:npw_sym))))
               ! 
             ENDIF ! nprocpool 

             !mapping is different with expand_kp, revert to the slow method
             DO ig=1, ngtot
             ! now for all G vectors find the PW coefficient for this k-point
             found = .FALSE.
             !!! MMORALES: This is very inefficient, create a mapping in the beggining from g
             !!!           to the g grid used in qmcpack, and just set to -1 the elements
             !!!           outside the cutoff
               DO ig7 = 1, npw_sym
                 IF( ig_l2g(igk_sym(ig7)) == igtog(ig) )THEN
                   !!! FIX FIX FIX, In parallel, this is completely incorrect since each proc only
                   !has limited data, you have to do a sum reduce at the very end to the head node 
                   eigpacked(ig)=tmp_evc(ig7)
                   found = .TRUE.
                   GOTO 18
                 ENDIF
               ENDDO ! ig7
             ! if can't find the coefficient this is zero
             18            IF( .NOT. found ) eigpacked(ig)=(0.d0,0.d0)
             ENDDO ! ig
           ELSE ! expandkp = false
             !
             !tmp_evc(:) = evc(:,ibnd)
             eigpacked(:)=(0.d0,0.d0)
             eigpacked(igtomin(ig_l2g(igk_k(1:npw,ik))))=evc(1:npw,ibnd)
             IF (debug) CALL check_norm(npw, evc(1,ibnd), .true., jks, ibnd, "before collection")
             !
           ENDIF ! expandkp

           CALL mp_sum ( eigpacked , intra_pool_comm )
           if (pool_ionode) THEN
             CALL check_norm(ngtot, eigpacked, .false., jks, ibnd, "after collection before writing")
             CALL esh5_write_psi_g(ibnd,eigpacked,ngtot)
           ENDIF

           IF (write_psir) THEN
              psic(:)=(0.d0,0.d0)
              psic(dffts%nl(ig_l2g(igk_sym(1:npw_sym))))=tmp_evc(1:npw_sym)
              if(gamma_only) psic(dffts%nlm(ig_l2g(igk_sym(1:npw_sym)))) = CONJG(tmp_evc(1:npw_sym))
              !
              CALL invfft ('Wave', psic, dffts)
              !
              IF(nproc_pool > 1) THEN
                ! 
                tmp_psic=psic
!                call cgather_smooth(psic,tmp_psic)
                call gather_grid(dffts,psic,tmp_psic)
                psiCptr => tmp_psic
                ! 
              ELSE
                ! 
                psiCptr => psic
                ! 
              ENDIF
              !
              IF(save_complex .eq. 1) THEN
                 !
                   !psic(:)=psic(:)/omega
                   ii=1
                   DO igx=1,nr1s
                      DO igy=0,nr2s-1
                         DO igz=0,nr3s-1
                            psitr(ii)=psiCptr(igx+nr1s*(igy+igz*nr2s))/omega
                            ii=ii+1
                         ENDDO
                      ENDDO
                   ENDDO
                   if(pool_ionode) CALL esh5_write_psi_r(ibnd,psitr,save_complex)
                 !
              ELSE
                 !
                   ii=1
                   DO igx=1,nr1s
                      DO igy=0,nr2s-1
                         DO igz=0,nr3s-1
                            psireal(ii)=real(psiCptr(igx+nr1s*(igy+igz*nr2s)))/omega
                            ii=ii+1
                         ENDDO
                      ENDDO
                   ENDDO
                   if(pool_ionode) CALL esh5_write_psi_r(ibnd,psireal,save_complex)
                 !
              ENDIF
           ENDIF ! write_psir
           !! conversion and output complete for each band
        ENDDO ! ibnd
        if(pool_ionode) CALL esh5_close_spin()
!      ENDDO
     if(pool_ionode) CALL esh5_close_kpoint()
   ENDDO ! iks
  ENDDO ! ik

endif ! nk
CALL stop_clock( 'big_loop' )
#endif
  ! write charge density
  ! ignore spin for the time being
  !CALL esh5_write_rho(rho,rhog(1,1),ngm)

#if defined(__HDF5) || defined(__HDF5_C)
  CALL start_clock( 'write_h5' )
  if(ionode) then
    CALL esh5_open_density(gint_den,ngm_g,nr1s,nr2s,nr3s)
    DO ispin = 1, nspin
       CALL esh5_write_density_g(ispin,den_g_global(1,ispin))
    ENDDO
    CALL esh5_close_density()
  endif 
   
  if(pool_ionode) CALL esh5_close_electrons()
  if(pool_ionode) CALL esh5_close_file()

  CALL mp_barrier( intra_image_comm )
  CALL stop_clock( 'write_h5' )
!     glue h5 together
  CALL start_clock( 'glue_h5' )
  if(ionode) then
    if(npool>1) then
      h5name = TRIM( prefix ) // '.pwscf.h5'
      tmp = TRIM( tmp_dir )//TRIM( h5name )
      h5len = LEN_TRIM(tmp)
      call esh5_join_all(tmp,h5len,npool)
    endif
  endif
  CALL stop_clock( 'glue_h5' )
#endif

  IF( ALLOCATED(igtog) ) DEALLOCATE (igtog)
  IF( ALLOCATED(igtomin) ) DEALLOCATE (igtomin)
  IF( ALLOCATED(indx) ) DEALLOCATE (indx)
  IF( ALLOCATED(eigpacked) ) DEALLOCATE (eigpacked)
  !IF( ALLOCATED(g_qmc) ) DEALLOCATE (g_qmc)
  !IF( ALLOCATED(g_cart) ) DEALLOCATE (g_cart)
  IF( ALLOCATED(psireal) ) DEALLOCATE (psireal)
  IF( ALLOCATED(psitr) ) DEALLOCATE (psitr)
  IF( ALLOCATED(tmp_psic) ) DEALLOCATE (tmp_psic)
  IF( ALLOCATED(num_irrep) ) DEALLOCATE (num_irrep)
  IF( ALLOCATED(xkfull_index) ) DEALLOCATE (xkfull_index)
  IF( ALLOCATED(sym_list) ) DEALLOCATE (sym_list)
  IF( ALLOCATED(xk_full_list) ) DEALLOCATE (xk_full_list)
  IF( ALLOCATED(rir) ) DEALLOCATE (rir)
  IF( ALLOCATED(igk_sym) ) DEALLOCATE (igk_sym)
  IF( ALLOCATED(g2kin_sym) ) DEALLOCATE (g2kin_sym)
  !DEALLOCATE (phase)

  CONTAINS

  SUBROUTINE generate_symmetry_equivalent_list()
  ! 
  ! Code taken mostly from PW/exx.f90
  !
  !------------------------------------------------------------------------
  !
  USE kinds, ONLY: DP
  USE cell_base,  ONLY : at
  USE lsda_mod,   ONLY : nspin
  USE klist,      ONLY : xk
  USE io_global,  ONLY : stdout, ionode
  !
  USE klist,      ONLY : nkstot
  USE io_global,            ONLY : stdout
  USE wvfct,                ONLY : nbnd, npwx, npw, wg, et
  USE klist,                ONLY : wk, ngk, nks
  USE symm_base,            ONLY : nsym, s, ftau
  USE lsda_mod,             ONLY: lsda
  use fft_base,             ONLY : dffts
!  use fft_interfaces,       ONLY : invfft

  !
  IMPLICIT NONE
  !
  integer       :: is, ik, ikq, iq, ns ,  nktot
  logical       :: xk_not_found
  real (DP)     :: sxk(3), dxk(3), xk_cryst(3), xkk_cryst(3)
  logical :: exst
  REAL (DP)         :: eps =1.d-8

  !
  ! find all k-points equivalent by symmetry to the points in the k-list
  !
   
  if(lsda)then
    nktot=nkstot/2
  else
    nktot=nkstot
  endif
  
  nkfull = 0
  do ik =1, nktot
    !
    num_irrep(ik) = 0 
    !
    ! isym=1 is the identity
    do is=1,nsym
        xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)
        sxk(:) = s(:,1,is)*xk_cryst(1) + &
                 s(:,2,is)*xk_cryst(2) + &
                 s(:,3,is)*xk_cryst(3)
        ! add sxk to the auxiliary list if it is not already present
        xk_not_found = .true.
        do ikq=1, nkfull 
           if (xk_not_found ) then
              dxk(:) = sxk(:)-xk_full_list(:,ikq) - nint(sxk(:)-xk_full_list(:,ikq))
              if ( abs(dxk(1)).le.eps .and. &
                   abs(dxk(2)).le.eps .and. &
                   abs(dxk(3)).le.eps ) xk_not_found = .false.
           end if
        end do
        if (xk_not_found) then
           nkfull = nkfull + 1
           num_irrep(ik) = num_irrep(ik) + 1
           xkfull_index(ik,num_irrep(ik)) = nkfull
           xk_full_list(:,nkfull) = sxk(:)
           sym_list(nkfull) = is  
        end if

        sxk(:) = - sxk(:)
        xk_not_found = .true.
        do ikq=1, nkfull 
           if (xk_not_found ) then
              dxk(:) = sxk(:)-xk_full_list(:,ikq) - nint(sxk(:)-xk_full_list(:,ikq))
              if ( abs(dxk(1)).le.eps .and. &
                   abs(dxk(2)).le.eps .and. &
                   abs(dxk(3)).le.eps ) xk_not_found = .false.
           end if
        end do
        if (xk_not_found) then
           nkfull = nkfull + 1
           num_irrep(ik) = num_irrep(ik) + 1
           xkfull_index(ik,num_irrep(ik)) = nkfull
           xk_full_list(:,nkfull) = sxk(:)
           sym_list(nkfull) = -is
        end if

     end do
  end do
  !
  ! transform kp list to cartesian again
  do ik=1,nkfull
    dxk(:) = bg(:,1)*xk_full_list(1,ik) + &
             bg(:,2)*xk_full_list(2,ik) + &
             bg(:,3)*xk_full_list(3,ik)
    xk_full_list(:,ik) = dxk(:)
  enddo
  !
!   if(ionode) then
!     print *,'Symmetry Inequivalent list of k-points:'
!     print *,'Total number: ',nkstot
!     do ik =1, nkstot
!       WRITE(*,'(i6,3(1x,f20.15))') ik, xk(1:3,ik) 
!     enddo
!     print *,'Full list of k-points (crystal):'
!     print *,'Total number of k-points: ',nkfull
!     print *,'IRREP, N, SYM-ID, KP: '
!     do ik =1, nkstot
!       do ns=1,num_irrep(ik)
!         WRITE(*,'(i6,i6,i6,3(1x,f20.15))') ik,ns,sym_list(xkfull_index(ik,ns)) & 
!          ,xk_full_list(1:3,xkfull_index(ik,ns))
!       enddo
!     enddo
!   endif
  !
  ! check symm operations
  !
!   do ikq =1,nkfull
!     is = abs(sym_list(ikq))
!     if ( mod (s (2, 1, is) * dffts%nr1, dffts%nr2) .ne.0 .or. &
!        mod (s (3, 1, is) * dffts%nr1, dffts%nr3) .ne.0 .or. &
!        mod (s (1, 2, is) * dffts%nr2, dffts%nr1) .ne.0 .or. &
!        mod (s (3, 2, is) * dffts%nr2, dffts%nr3) .ne.0 .or. &
!        mod (s (1, 3, is) * dffts%nr3, dffts%nr1) .ne.0 .or. &
!        mod (s (2, 3, is) * dffts%nr3, dffts%nr2) .ne.0 ) then
!      call errore ('generate_symmetry_equivalent_list',' problems with grid',is)
!     end if
!   end do

  END SUBROUTINE generate_symmetry_equivalent_list 
  !
  SUBROUTINE generate_symmetry_rotation(is0) 
  USE kinds, ONLY: DP
  USE klist,      ONLY : xk
  USE io_global,  ONLY : stdout, ionode
  !
  USE io_global,            ONLY : stdout
  USE symm_base,            ONLY : nsym, s, ftau
  use fft_base,             ONLY : dffts

  !
  IMPLICIT NONE
  !
  integer, intent(in)  :: is0
  !
  integer       :: i,j,k, ir, ri, rj, rk, is
  logical :: exst
  REAL (DP)         :: eps =1.d-6

  !
  do ir=1, nxxs
    rir(ir) = ir
  end do
  is = abs(is0)
  do k = 1, dffts%nr3
   do j = 1, dffts%nr2
    do i = 1, dffts%nr1
      call ruotaijk (s(1,1,is), ftau(1,is), i, j, k, &
              dffts%nr1,dffts%nr2,dffts%nr3, ri, rj , rk )
      ir =   i + ( j-1)*dffts%nr1x + ( k-1)*dffts%nr1x*dffts%nr2x
      rir(ir) = ri + (rj-1)*dffts%nr1x + (rk-1)*dffts%nr1x*dffts%nr2x
    end do
   end do
  end do
  !
  END SUBROUTINE generate_symmetry_rotation 
  !
END SUBROUTINE compute_qmcpack
