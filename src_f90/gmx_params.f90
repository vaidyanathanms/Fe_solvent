! Parameters for GMX trajectories

!####################################################################
! Version 2.0: November-05-2025
! Author: Vaidyanathan Sethuraman
! Email: vm5@ornl.gov
!####################################################################

!####################################################################
! Analyzes GROMACS outputs in .gro format
! Can do analysis based on COM or on individual atoms
! Auxiliary files: analyze_fe.f90; anainp.txt
!####################################################################

MODULE PARAMS_GMX

  USE OMP_LIB
  IMPLICIT NONE

  ! Atom, processor and time details
  INTEGER :: ntotatoms, ntotatomtypes, ntotmols
  INTEGER :: nframes, freqfr, nfrcntr, skipfr
  REAL    :: start_time
  INTEGER :: nproc
  INTEGER :: ioncnt,  iontype
  INTEGER :: c_ioncnt, c_iontype
  INTEGER :: nmol_totcom, ncom_types
  INTEGER :: nclust_types
  INTEGER, PARAMETER :: com_type = 199 ! Set here
  
  ! Structural quantities
  INTEGER :: rdffreq,rmaxbin,npairs,rdfpaircnt
  INTEGER :: maxneighsize, neighfreq
  INTEGER :: ntotion_centers,totmult_centers
  INTEGER :: maxsize_species
  INTEGER :: neighfreq_COM,maxneighsize_COM
  REAL    :: rneigh_cut,rcatan_cut,rclust_cut,rcatCOM_cut
  REAL    :: rvolavg,rdomcut,rbinval,rvolval
  REAL    :: rCOMneigh_cut
  
  ! Dynamic quantities
  REAL    :: delta_t,q_targ,q_tol,q_targ_max,q_targ_min
  REAL    :: q_bin,  qval, q_hi, q_lo, q_hi2, q_lo2
  INTEGER :: q_nmax, q_possN
  
  ! All flags
  INTEGER :: box_from_file_flag, box_type_flag
  INTEGER :: rdfcalc, rgcalc, rgall, rgavg, rgfreq
  INTEGER :: rdfcalc_flag
  INTEGER :: catan_neighcalc_flag
  INTEGER :: clust_calc_flag, clust_time_flag
  INTEGER :: comflag,catCOM_neighcalc_flag
  INTEGER :: ion_dynflag, cion_dynflag, com_dynflag
  INTEGER :: ion_diff,cion_diff,com_diff
  INTEGER :: catan_autocfflag, catCOM_autocfflag
  INTEGER :: name_to_type_map_flag
  INTEGER :: multclust_calc_flag
  INTEGER :: dynfsktflag
  
  ! File names and unit numbers
  CHARACTER(LEN = 256) :: ana_fname, dum_fname
  CHARACTER(LEN = 256) :: data_fname,traj_fname,log_fname,box_fname
  CHARACTER(LEN = 256) :: rdf_fname
  INTEGER, PARAMETER :: anaread = 2, logout = 3
  INTEGER, PARAMETER :: trajread = 15, boxread = 20, inpread = 100
  INTEGER, PARAMETER :: dumwrite = 50
  INTEGER, PARAMETER :: max_char = 5
  INTEGER, PARAMETER :: clustwrite = 150
  
  !Math constants
  REAL*8, PARAMETER :: pival  = 3.14159265359
  REAL*8, PARAMETER :: pi2val = 2.0*pival

  !Trajectory file read details
  REAL :: box_xl,box_yl,box_zl, boxval
  REAL :: box_xf,box_yf,box_zf
  INTEGER :: nbox_steps, timestep
  REAL*8 :: act_time

  !Required global arrays
  CHARACTER(max_char),ALLOCATABLE :: name_arr(:) 
  REAL,ALLOCATABLE,DIMENSION(:) :: boxx_arr,boxy_arr,boxz_arr
  REAL,ALLOCATABLE,DIMENSION(:,:) :: masses
  REAL*8,ALLOCATABLE,DIMENSION(:) :: com_masses
  REAL*8,ALLOCATABLE,DIMENSION(:,:) :: rxyz_gmx,comxyz_gmx
  INTEGER,ALLOCATABLE,DIMENSION(:) :: type_arr,comtyp_arr
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: aidvals
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ion_IDTYP_arr&
       &,countion_IDTYP_arr, com_IDTYP_arr
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: allionids,multionids

  !Required arrays - Statics Quantities
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: pairs_rdf
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: mclust_type_arr
  REAL*8,ALLOCATABLE,DIMENSION(:) :: clust_avg, spec_avg
  REAL*8,ALLOCATABLE,DIMENSION(:,:):: rdfarray
  REAL*8,ALLOCATABLE,DIMENSION(:,:):: mclust_rcut_arr
  REAL*8,ALLOCATABLE,DIMENSION(:) :: cat_an_neighavg,an_cat_neighavg

  !Required Arrays - Dynamic Quantities
  INTEGER*8, ALLOCATABLE, DIMENSION(:,:) :: q_nlist
  REAL*8, ALLOCATABLE, DIMENSION(:)   :: tarr_gmx
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: trx_gmx,try_gmx,trz_gmx
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: itrx_gmx,itry_gmx,itrz_gmx
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ctrx_gmx,ctry_gmx,ctrz_gmx
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: comtx_gmx,comty_gmx,comtz_gmx

END MODULE PARAMS_GMX
