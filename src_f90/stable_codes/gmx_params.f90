! Parameters for GMX trajectories

MODULE PARAMS_GMX

  USE OMP_LIB
  IMPLICIT NONE

  ! Atom, processor and time details
  INTEGER :: ntotatoms, ntotatomtypes
  INTEGER :: nframes, freqfr, nfrcntr, skipfr
  REAL    :: start_time
  INTEGER :: nproc
  INTEGER :: ioncnt,  iontype
  INTEGER :: c_ioncnt, c_iontype
  INTEGER :: p_ioncnt, p_iontype, npoly_types
  INTEGER :: nclust_types
  
  ! Structural quantities
  INTEGER :: rdffreq,rmaxbin,npairs,rdfpaircnt
  INTEGER :: maxneighsize, neighfreq
  INTEGER :: ntotion_centers,totmult_centers
  INTEGER :: maxsize_species
  REAL    :: rneigh_cut,rcatan_cut,rclust_cut,rcatpol_cut
  REAL    :: rvolavg,rdomcut,rbinval,rvolval

  ! Dynamic quantities
  REAL    :: delta_t
  
  ! All flags
  INTEGER :: box_from_file_flag, box_type_flag
  INTEGER :: rdfcalc, rgcalc, rgall, rgavg, rgfreq
  INTEGER :: rdfcalc_flag
  INTEGER :: catan_neighcalc_flag
  INTEGER :: clust_calc_flag, clust_time_flag
  INTEGER :: polyflag
  INTEGER :: ion_dynflag, cion_dynflag, pion_dynflag
  INTEGER :: ion_diff,cion_diff,pion_diff
  INTEGER :: catan_autocfflag, catpol_autocfflag
  INTEGER :: name_to_type_map_flag
  INTEGER :: multclust_calc_flag
  
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
  REAL*8,ALLOCATABLE,DIMENSION(:,:) :: rxyz_lmp
  INTEGER,ALLOCATABLE,DIMENSION(:) :: type_arr,polytyp_arr
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: aidvals
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ionarray,counterarray
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: polyionarray
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: allionids,multionids

  !Required arrays - Statics Quantities
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: pairs_rdf
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: mclust_type_arr
  REAL*8,ALLOCATABLE,DIMENSION(:) :: clust_avg, spec_avg
  REAL*8,ALLOCATABLE,DIMENSION(:,:):: rdfarray
  REAL*8,ALLOCATABLE,DIMENSION(:,:):: mclust_rcut_arr
  REAL*8,ALLOCATABLE,DIMENSION(:) :: cat_an_neighavg,an_cat_neighavg

  !Required Arrays - Dynamic Quantities
  REAL*8, ALLOCATABLE, DIMENSION(:)   :: tarr_lmp
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: trx_lmp,try_lmp,trz_lmp
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: itrx_lmp,itry_lmp,itrz_lmp
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ctrx_lmp,ctry_lmp,ctrz_lmp
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ptrx_lmp,ptry_lmp,ptrz_lmp

END MODULE PARAMS_GMX
