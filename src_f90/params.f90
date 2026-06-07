!---------------To analyze properties of polymer-ion systems---------
!---------------Version 2: July-15-2024------------------------------
!---------------Main File: analyze.f90-------------------------------
!********************************************************************

MODULE ANALYZE_ALL

  USE OMP_LIB
  IMPLICIT NONE

  ! Required Input Variables
  INTEGER :: initdist
  INTEGER :: nframes, skipfr, freqfr, nfrcntr
  INTEGER :: nchains, atperchain
  INTEGER :: nproc

  !Structural adn dynamics analysis averages/counters/inputs
  INTEGER :: rdffreq,rmaxbin,npairs
  INTEGER :: rgfreq, rdfpaircnt
  REAL    :: rvolavg,rdomcut,rbinval,rvolval
  REAL    :: re2ave, re4ave, rg2ave, rg4ave, b2ave
  REAL    :: delta_t

  !General required vars required for computing properties
  INTEGER :: npoly_types
  INTEGER :: ioncnt,  iontype
  INTEGER :: c_ioncnt, c_iontype
  INTEGER :: p_ioncnt, p_iontype
  INTEGER :: maxneighsize, neighfreq
  INTEGER :: ntotion_centers
  REAL    :: rneigh_cut,rcatan_cut,rcatpol_cut1,rcatpol_cut2
  
  ! All analysis flags
  INTEGER :: polyflag
  INTEGER :: rdfcalc, rgcalc, rgall, rgavg
  INTEGER :: catan_neighcalc
  INTEGER :: bfrdf_calc
  INTEGER :: clust_calc
  INTEGER :: ion_dynflag, cion_dynflag, pion_dynflag
  INTEGER :: ion_diff, cion_diff, pion_diff
  INTEGER :: catan_autocfflag, catpol_autocfflag

  ! File names and unit numbers
  CHARACTER(LEN = 256) :: ana_fname,data_fname,traj_fname,log_fname
  CHARACTER(LEN = 256) :: rdf_fname, dum_fname
  INTEGER, PARAMETER :: anaread = 2,   logout = 3
  INTEGER, PARAMETER :: dumwrite = 50
  INTEGER, PARAMETER :: inpread = 100, rgwrite = 200,rgavgwrite = 300

  !Math Constants
  REAL*8, PARAMETER :: pival  = 3.14159265359
  REAL*8, PARAMETER :: pi2val = 2.0*pival

  !Global analysis variables and arrays
  INTEGER :: atomflag, velflag, bondflag, anglflag, dihdflag,imprflag
  INTEGER :: ntotatoms, ntotbonds, ntotangls,ntotdihds,ntotimprs
  INTEGER :: ntotatomtypes,ntotbondtypes,ntotangltypes,ntotdihdtypes&
       &,ntotimprtypes

  !Lammps trajectory file read details
  REAL :: box_xl,box_yl,box_zl, boxval
  INTEGER*8 :: timestep

  !Required Arrays - LAMMPS
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rxyz_lmp, vel_xyz, charge_lmp&
       &,masses
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bond_lmp, angl_lmp,&
       & dihd_lmp, impr_lmp,aidvals
  CHARACTER,ALLOCATABLE,DIMENSION(:,:) :: keywords
  REAL,ALLOCATABLE,DIMENSION(:):: boxx_arr, boxy_arr,boxz_arr

  !Required Arrays - Dynamic Quantities
  INTEGER*8, ALLOCATABLE, DIMENSION(:) :: tarr_lmp
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: itrx_lmp,itry_lmp,itrz_lmp

END MODULE ANALYZE_ALL

!--------------------------------------------------------------------
