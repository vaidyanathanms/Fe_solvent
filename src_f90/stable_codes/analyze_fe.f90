! Main program for analyzing GROMACS trajectories

!####################################################################
! Version 2.0: November-05-2025
! Author: Vaidyanathan Sethuraman
! Email: vm5@ornl.gov
!####################################################################

!####################################################################
! Analyzes GROMACS outputs in .gro format
! Can do analysis based on COM or on individual atoms
! Auxiliary files: gmx_params.f90; anainp.txt
!####################################################################
PROGRAM ANALYZE_GMXTRAJ

  USE PARAMS_GMX
  IMPLICIT NONE
  
  ! Print headers
  PRINT *, "Analysis of Fe-salt systems .."
  PRINT *, "Starting OMP Threads .."
  
!$OMP PARALLEL
  nproc = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
  PRINT *, "Number of threads: ", nproc
  
  CALL DEFAULTVALUES()
  CALL READ_ANA_INP_FILE()
  CALL READ_DATAFILE()
  CALL SORT_ION_CION_ARR()
  CALL ALLOCATE_COM_ARRAYS()
  CALL ALLOCATE_ANALYSIS_ARRAYS()
  CALL ANALYZE_TRAJECTORYFILE() 
  CALL ALLOUTPUTS()
  CALL DEALLOCATE_ARRAYS()

  PRINT *, "Completed all analysis .."

END PROGRAM ANALYZE_GMXTRAJ

!---------------------------------------------------------------------

SUBROUTINE READ_ANA_INP_FILE()

  USE PARAMS_GMX
  IMPLICIT NONE
  
  INTEGER :: nargs,ierr,logflag,AllocateStatus,i,j
  INTEGER :: ncut_offs,type_a,type_b, a_ind, b_ind
  REAL    :: rcab_cut_val
  CHARACTER(100) :: fname_pref
  CHARACTER(256) :: dumchar
  CHARACTER(max_char) :: aname
  
  CALL DEFAULTVALUES()

  nargs = IARGC()
  IF(nargs .NE. 1) STOP "Input incorrect"

  logflag = 0

  CALL GETARG(nargs,ana_fname)

  OPEN(unit = anaread,file=trim(ana_fname),action="read",status="old"&
       &,iostat=ierr)
  
  IF(ierr /= 0) THEN

     PRINT *, trim(ana_fname), "not found"
     STOP

  END IF

  DO

     READ(anaread,*,iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT

     ! Read file and trajectory details
     IF(dumchar == 'datafile') THEN
        
        READ(anaread,*,iostat=ierr) data_fname

     ELSEIF(dumchar == 'trajectory_file') THEN

        READ(anaread,*,iostat=ierr) traj_fname

     ELSEIF(dumchar == 'nframes') THEN

        READ(anaread,*,iostat=ierr) nframes

     ELSEIF(dumchar == 'start_time') THEN

        READ(anaread,*,iostat=ierr) start_time

     ELSEIF(dumchar == 'freqfr') THEN

        READ(anaread,*,iostat=ierr) freqfr

     ! Ion/counterion/COM type definitions
     ELSEIF(dumchar == 'ion_type') THEN
        
        READ(anaread,*,iostat=ierr) iontype
        
     ELSEIF(dumchar == 'cion_type') THEN

        READ(anaread,*,iostat=ierr) c_iontype

        
     ELSEIF(dumchar == 'name_to_type_map') THEN

        READ(anaread,*,iostat=ierr) ntotatomtypes
        
        ALLOCATE(name_arr(ntotatomtypes),stat=AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate name_arr"
        ALLOCATE(type_arr(ntotatomtypes),stat=AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate type_arr"

        DO i = 1, ntotatomtypes

           READ(anaread,*,iostat=ierr) name_arr(i), type_arr(i)

        END DO
        
        name_to_type_map_flag = 1

     ELSEIF(dumchar == 'com_types') THEN 
        
        READ(anaread,*,iostat=ierr) ncom_types
        
        ALLOCATE(comtyp_arr(ncom_types),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate comtyp_arr"
        
        READ(anaread,*,iostat=ierr) (comtyp_arr(i),i=1,ncom_types)
        comflag = 1

     !Here onwards static properties
     ELSEIF(dumchar == 'compute_rdf') THEN

        rdfcalc_flag = 1
        READ(anaread,*,iostat=ierr) rdffreq, rmaxbin, rdomcut,npairs
        
        ALLOCATE(pairs_rdf(npairs,3),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate pairs_rdf"
      
        DO i = 1,npairs

           READ(anaread,*,iostat=ierr) pairs_rdf(i,1), pairs_rdf(i,2)

           IF(pairs_rdf(i,1) == com_type .OR. pairs_rdf(i,2) ==&
                & com_type) THEN

              IF(comflag /= 1) THEN

                 PRINT *, "Pair type is COM_type without defining com_&
                      &types..."
                 PRINT *, "ERROR: Define com_types..."
                 STOP

              END IF
              
           END IF
           
        END DO

     ELSEIF(dumchar == 'compute_clust') THEN

        clust_calc_flag = 1
        READ(anaread,*,iostat=ierr) rclust_cut, clust_time_flag
        
     ELSEIF(dumchar == 'compute_catanneigh') THEN

        catan_neighcalc_flag = 1
        READ(anaread,*,iostat=ierr) neighfreq,maxneighsize,rneigh_cut
                   
     ELSEIF(dumchar == 'compute_multtype_clust') THEN

        multclust_calc_flag = 1
        READ(anaread,*,iostat=ierr) nclust_types
        ncut_offs = INT(nclust_types*(nclust_types+1)/2)
        
        ALLOCATE(mclust_type_arr(nclust_types,2),stat=AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate mclust_type_arr"
        mclust_type_arr = 0 ! DO NOT INITIALIZE TO ANY OTHER NUMBER

        DO i = 1, nclust_types
           READ(anaread,*,iostat=ierr) type_a
           mclust_type_arr(i,1) = type_a
        END DO ! 2nd element is the number of each type
           
        ALLOCATE(mclust_rcut_arr(nclust_types,nclust_types),stat&
             &=AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate mclust_rcut_arr"

        DO i = 1, ncut_offs

           READ(anaread,*,iostat=ierr) type_a, type_b, rcab_cut_val
           CALL MAP_TYPE_TO_INDEX(type_a,a_ind)
           CALL MAP_TYPE_TO_INDEX(type_b,b_ind)
           mclust_rcut_arr(a_ind,b_ind) = rcab_cut_val
           mclust_rcut_arr(b_ind,a_ind) = rcab_cut_val

        END DO

     ! Here onwards dynamic properties
     ELSEIF(dumchar == 'compute_iondiff') THEN

        READ(anaread,*,iostat=ierr) ion_diff, delta_t
        ion_dynflag = 1
                
     ELSEIF(dumchar == 'compute_ciondiff') THEN

        READ(anaread,*,iostat=ierr) cion_diff, delta_t
        cion_dynflag = 1

     ELSEIF(dumchar == 'compute_catanrestime') THEN
        
        READ(anaread,*,iostat=ierr) rcatan_cut
        catan_autocfflag = 1
        ion_dynflag = 1; cion_dynflag = 1

     ELSEIF(dumchar == 'compute_dynfskt') THEN

        READ(anaread,*,iostat=ierr) q_targ_min,q_targ_max,q_bin,q_tol&
             &,q_nmax
        ion_dynflag = 1; cion_dynflag = 1; dynfsktflag = 1
        
     ! Here onwards with respect to COM of the polymer/solvent
     ! For COM of polymer/solvent

     ELSEIF(dumchar == 'compute_comdiff') THEN 
        IF (comflag .NE. 1) STOP "Should define com_types first .."
        READ(anaread,*,iostat=ierr) com_diff, delta_t
        com_dynflag = 1

     ELSEIF(dumchar == 'compute_catcomrestime') THEN
        IF (comflag .NE. 1) STOP "Should define com_types first .."
        READ(anaread,*,iostat=ierr) rcatCOM_cut
        catCOM_autocfflag = 1; com_dynflag =1; ion_dynflag = 1

     ! Read log filename
     ELSEIF(dumchar == 'log_file') THEN

        READ(anaread,*,iostat=ierr) log_fname
        logflag  = 1

     ELSE
        
        PRINT *, "unknown keyword: ", trim(dumchar)
        STOP

     END IF

  END DO

  IF(logflag == 0) THEN
     
     WRITE(fname_pref,'(A11,I0,I0,A1)') "log_",iontype,c_iontype,'_'
     log_fname  = trim(adjustl(fname_pref))//trim(adjustl(traj_fname))

  END IF
  
  OPEN(unit = logout,file=trim(log_fname),action="write",status="repla&
       &ce",iostat=ierr)

  PRINT *, "Analysis input file read finished .."

  CALL SANITY_CHECK_IONTYPES()
  
END SUBROUTINE READ_ANA_INP_FILE

!--------------------------------------------------------------------

SUBROUTINE DEFAULTVALUES()

  USE PARAMS_GMX
  IMPLICIT NONE

  ! Initialize flags
  rgall = 0; rgcalc = 0; rdfcalc = 0
  ion_dynflag = 0; cion_dynflag = 0; com_dynflag = 0
  ion_diff = 0; cion_diff = 0
  com_diff = 0; comflag = 0
  catan_autocfflag = 0; catCOM_autocfflag = 0
  dynfsktflag = 0
  
  ! Initialize iontypes
  c_iontype = -1; iontype = -1

  !Initialize system quantities
  ncom_types = 0; ioncnt = 0; c_ioncnt = 0; nmol_totcom= 0

  ! Initialize distributions and frequencies
  rdffreq = 0; rgfreq = 1

  ! Initialzie structural quantities
  rdomcut = 10.0;  rmaxbin = 100; rbinval = REAL(rdomcut)&
       &/REAL(rmaxbin)
  rcatan_cut = 0.0; rneigh_cut = 0.0
  
  ! Initialize structural averages
  rvolavg = 0; rgavg = 0

  ! Initialize dynamical quantities
  rcatCOM_cut = 0.0
  q_targ_min = 0.0; q_targ_max = 0.0; q_bin = 0.0; q_tol = 0.0
  q_nmax = 5; q_targ = 0
  
END SUBROUTINE DEFAULTVALUES

!--------------------------------------------------------------------

SUBROUTINE READ_DATAFILE()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,j,ierr,k,u,AllocateStatus,imax,aid
  INTEGER :: flag, cntr, nwords
  INTEGER :: molid,atype,ix,iy,iz,mtype
  REAL    :: charge,rx,ry,rz,massval
  REAL    :: xlo,xhi,ylo,yhi,zlo,zhi
  CHARACTER(LEN=256):: headline
  LOGICAL :: ok_t, ok_step
  CHARACTER(LEN=max_char):: molname,aname
  CHARACTER(len=256) :: rline,dumchar

  OPEN(unit=inpread,file = trim(data_fname),action =&
       & 'read', status='old',iostat=ierr) 
  
  IF(ierr .NE. 0) STOP "Data file not found"

  WRITE(logout,*) "Datafile used: ", trim(adjustl(data_fname))

  ntotatoms = 0
  READ(inpread,'(A)') headline
  READ(inpread,*) ntotatoms

  CALL GET_VALUE_REAL(trim(adjustl(headline)),'t',act_time,ok_t)
  CALL GET_VALUE_INT(trim(adjustl(headline)),'step',timestep&
       &,ok_step)

  IF(ok_t == .FALSE. .OR. ok_step == .FALSE.) THEN
     
     PRINT *, "Error in reading datafile .."
     PRINT *, trim(adjustl(headline))
     STOP
     
  END IF

  PRINT *, "STATISTICS"
  PRINT *, "Number of atoms/atomtypes: " , ntotatoms,ntotatomtypes
  flag = 0; cntr = 0

  CALL ALLOCATE_TOPO_ARRAYS()
         
  masses = -1

  DO j = 1,ntotatoms
        
     READ(inpread,'(i5,2a5,i5,3f8.3,3f8.4)') molid,molname,aname&
          &,aid,rxyz_gmx(aid,1),rxyz_gmx(aid,2),rxyz_gmx(aid,3)
 
     CALL MAP_ANAME_TO_ATYPE(aname,atype,k,0)
     
     aidvals(j,1)     = j
     aidvals(j,2)     = molid
     aidvals(j,3)     = atype
     rxyz_gmx(j,1)    = rx
     rxyz_gmx(j,2)    = ry
     rxyz_gmx(j,3)    = rz
     
     IF(masses(k,1) == -1) THEN
        
        CALL ASSIGN_MASSES(aname,atype,j,massval)
        masses(k,1) = INT(atype)
        masses(k,2) = massval
        
     END IF
        
  END DO

  ntotmols = MAXVAL(aidvals(:,2))
  PRINT *, "Total number of distinct molecules ..", ntotmols
  WRITE(logout,*) "Total number of distinct molecules ..", ntotmols
  
  PRINT *, "Writing mass and type data..."
  DO i = 1,ntotatomtypes
     WRITE(logout,*) name_arr(i),type_arr(i),masses(i,1),masses(i,2)
  END DO
     
  PRINT *, "Datafile read completed..."

  IF(comflag) THEN
     PRINT *, "Generating COM arrays .."
     CALL SORT_COM_ARR()
  END IF
  
  CLOSE(inpread)

END SUBROUTINE READ_DATAFILE

!--------------------------------------------------------------------

SUBROUTINE MAP_ANAME_TO_ATYPE(aname,atype,kcnt,tval)

  USE PARAMS_GMX
  IMPLICIT NONE
  
  CHARACTER(len=*), INTENT(IN) :: aname
  REAL, INTENT(IN) :: tval
  INTEGER, INTENT(OUT) :: atype,kcnt
  INTEGER :: atom_find_flag
  
  atom_find_flag = -1

  DO kcnt = 1,ntotatomtypes

     IF(trim(adjustl(aname)) == trim(adjustl(name_arr(kcnt)))) THEN
        
        atype = type_arr(kcnt)
        atom_find_flag = 1
        EXIT
        
     END IF
     
  END DO

  IF(atom_find_flag == -1) THEN

     
     PRINT *, "ERROR: Unknown atom name ", aname, " at t= ", tval
     STOP

  END IF
    

END SUBROUTINE MAP_ANAME_TO_ATYPE

!--------------------------------------------------------------------

SUBROUTINE MAP_TYPE_TO_INDEX(atype,a_index)

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: atype
  INTEGER, INTENT(OUT) :: a_index
  INTEGER :: i
  
  a_index = -1

  DO i = 1,nclust_types

     IF(atype == mclust_type_arr(i,1)) THEN

        a_index = i
        EXIT

     END IF

  END DO

  IF(a_index == -1) THEN

     PRINT *, "Unknown atom type in cut-off data", atype
     STOP

  END IF
     
END SUBROUTINE MAP_TYPE_TO_INDEX

!--------------------------------------------------------------------

SUBROUTINE ASSIGN_MASSES(aname,atype,aid,massval)

  USE PARAMS_GMX
  IMPLICIT NONE

  CHARACTER(len=*), INTENT(IN) :: aname
  INTEGER, INTENT(IN) :: atype, aid
  REAL, INTENT(OUT) :: massval

  massval = 1.0 ! Assign unity as default
  IF(trim(adjustl(aname)) == 'Al') THEN
     
     massval = 26.981539

  ELSEIF(trim(adjustl(aname)) == 'Cl') THEN
     
     massval = 35.453
     
  ELSEIF(trim(adjustl(aname)) == 'H' .OR. trim(adjustl(aname)) ==&
       & 'HW' .OR. trim(adjustl(aname)) == 'HW1' .OR.&
       & trim(adjustl(aname)) == 'HW2') THEN
     
     massval = 1.00784
     
  ELSEIF(trim(adjustl(aname)) == 'C' .OR. trim(adjustl(aname)) == '1C&
       &' .OR. trim(adjustl(aname)) == '2C') THEN
     
     massval = 12.011
     
  ELSEIF(trim(adjustl(aname)) == 'O' .OR. trim(adjustl(aname)) ==&
       & '1O' .OR. trim(adjustl(aname)) == '2O' .OR.&
       & trim(adjustl(aname)) == '3O' .OR. trim(adjustl(aname)) == '4&
       &O' .OR. trim(adjustl(aname)) == 'OW') THEN
     
     massval = 15.999
     
  ELSEIF(trim(adjustl(aname)) == 'N' .OR. trim(adjustl(aname)) ==&
        & '1N') THEN
     
     massval = 14.0067

  ELSEIF(trim(adjustl(aname)) == 'S' .OR.  trim(adjustl(aname)) ==&
        & '1S' .OR. trim(adjustl(aname)) == '2S') THEN
     
     massval = 32.0650
     
  ELSEIF(trim(adjustl(aname)) == 'F' .OR. trim(adjustl(aname)) ==&
       & '1F' .OR. trim(adjustl(aname)) == '2F' .OR.&
       & trim(adjustl(aname)) == '3F' .OR. trim(adjustl(aname)) == '4&
       &F' .OR. trim(adjustl(aname)) == '5F' .OR.&
       & trim(adjustl(aname)) == '6F') THEN
     
     massval = 18.998
     
  ELSEIF(trim(adjustl(aname)) == 'FE' .OR. trim(adjustl(aname)) ==  '&
       &FE2P' .OR. trim(adjustl(aname)) == 'FE3P' .OR.&
       & trim(adjustl(aname)) == 'FE2') THEN
     
     massval = 55.845
     
  ELSE
     
     PRINT *, trim(adjustl(aname)), " and ID ", aid, " not found!"
     PRINT *, "WARNING: Assigning unit mass to ", atype
     
  END IF
  
END SUBROUTINE ASSIGN_MASSES

!--------------------------------------------------------------------

SUBROUTINE FIND_MASS_FROM_TYPE(atype,dmassval)

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: atype
  REAL, INTENT(OUT) :: dmassval
  INTEGER :: mcnt, mass_find_flag 

  mass_find_flag = -1

  DO mcnt = 1, ntotatomtypes

     IF(atype == masses(mcnt,1)) THEN

        dmassval = masses(mcnt,2)
        mass_find_flag = 1
        EXIT

     END IF

  END DO

  IF(mass_find_flag == -1) THEN

     PRINT *, "ERROR: Unidentified type for finding mass .."
     PRINT *, atype, type_arr
     STOP

  END IF
        

END SUBROUTINE FIND_MASS_FROM_TYPE

!--------------------------------------------------------------------


SUBROUTINE SORT_ION_CION_ARR()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,j,a1type,cnt,AllocateStatus,ntotion_cnt,aid,molid
  CHARACTER(100) :: fname_pref
  INTEGER, DIMENSION(1:ntotatoms,2) :: dumionarr,dumcionarr

  dumionarr = -1; dumcionarr = -1
  cnt = 0
  ntotion_cnt = 0

  DO i = 1,ntotatoms

     a1type = aidvals(i,3)

     IF(a1type == iontype) THEN
        ioncnt = ioncnt + 1
        dumionarr(ioncnt,1) = i
        dumionarr(ioncnt,2) = a1type

     ELSEIF(a1type == c_iontype) THEN
        c_ioncnt = c_ioncnt + 1
        dumcionarr(c_ioncnt,1) = i
        dumcionarr(c_ioncnt,2) = a1type
        
     END IF    

  END DO

  ! Always identify ion and counter-ion types
  PRINT *, "Number of atoms of ion type: ", ioncnt
  PRINT *, "Number of atoms of cntion type: ", c_ioncnt

  ALLOCATE(ion_IDTYP_arr(ioncnt,2),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate ion_IDTYP_arr"
  ALLOCATE(countion_IDTYP_arr(c_ioncnt,2),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate countion_IDTYP_arr"

  ! Load ion array
  i = 0

  DO WHILE(dumionarr(i+1,1) .NE. -1)

     i = i + 1
     ion_IDTYP_arr(i,1) = dumionarr(i,1)
     ion_IDTYP_arr(i,2) = dumionarr(i,2)

  END DO

    IF(i .NE. ioncnt) THEN
     PRINT *, i, ioncnt
     STOP "Wrong total count in ion_IDTYP_arr"
  END IF

  DO i = 1,ioncnt

     IF(ion_IDTYP_arr(i,1) == -1 .OR. ion_IDTYP_arr(i,2) == -1) THEN

        PRINT *, i,ion_IDTYP_arr(i,1), ion_IDTYP_arr(i,2)
        PRINT *, "Something wrong in assigning ion_IDTYP_arr"
        STOP

     END IF

     IF(ion_IDTYP_arr(i,2) .NE. iontype) THEN

        PRINT *, i,ion_IDTYP_arr(i,1), ion_IDTYP_arr(i,2)
        PRINT *, "Something wrong in ion_IDTYP_arr type"
        STOP

     END IF

  END DO

  WRITE(fname_pref,'(A11,I0,A4)') "iontype_",iontype,'.txt'
  dum_fname  = trim(adjustl(fname_pref))
  OPEN(unit = 93,file=dum_fname,action="write",status="replace")

  WRITE(93,*) "Reference type/count: ", iontype, ioncnt

  DO i = 1,ioncnt
     WRITE(93,'(3(I0,1X))') i, ion_IDTYP_arr(i,1), ion_IDTYP_arr(i,2)
  END DO

  CLOSE(93)

  ! Load counterion array

  i = 0

  DO WHILE(dumcionarr(i+1,1) .NE. -1)

     i = i + 1
     countion_IDTYP_arr(i,1) = dumcionarr(i,1)
     countion_IDTYP_arr(i,2) = dumcionarr(i,2)

  END DO

  IF(i .NE. c_ioncnt) THEN
     PRINT *, i, c_ioncnt
     STOP "Wrong total count in countion_IDTYP_arr"
  END IF

  DO i = 1,c_ioncnt

     IF(countion_IDTYP_arr(i,1) == -1 .OR. countion_IDTYP_arr(i,2) ==&
          & -1) THEN

        PRINT *, i,countion_IDTYP_arr(i,1), countion_IDTYP_arr(i,2)
        PRINT *, "Something wrong in assigning countion_IDTYP_arr"
        STOP

     END IF

     IF(countion_IDTYP_arr(i,2) .NE. c_iontype) THEN

        PRINT *, i,countion_IDTYP_arr(i,1), countion_IDTYP_arr(i,2)
        PRINT *, "Something wrong in counterion_IDTYP_arr type"
        STOP

     END IF

  END DO
  
  WRITE(fname_pref,'(A11,I0,A4)') "ciontype_",c_iontype,'.txt'
  dum_fname  = trim(adjustl(fname_pref))
  OPEN(unit = 93,file=dum_fname,action="write",status="replace")

  WRITE(93,*) "Reference type/count: ", c_iontype, c_ioncnt

  DO i = 1,c_ioncnt

     WRITE(93,'(3(I0,1X))') i, countion_IDTYP_arr(i,1),&
          & countion_IDTYP_arr(i,2)

  END DO

  CLOSE(93)

! COM of polymer/solvent ion array required only for diffusion systems
   ! Cluster calc requires to add iontype and c_iontype in same array
  IF (clust_calc_flag) THEN

     ntotion_centers = ioncnt + c_ioncnt
     PRINT *, "Total number of ion centers: ", ntotion_centers
     cnt = 1

     ALLOCATE(allionids(ntotion_centers,2),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate allionids"
     ALLOCATE(clust_avg(ntotion_centers),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate clust_avg"

     allionids = 0 ! Initial allocation

     DO i = 1,ntotatoms

        a1type = aidvals(i,3)

        IF(a1type == iontype .OR. a1type == c_iontype) THEN

           allionids(cnt,1) = i
           allionids(cnt,2) = a1type
           cnt = cnt + 1

        END IF

     END DO

  ELSE

     ALLOCATE(allionids(1,2),stat = AllocateStatus)
     DEALLOCATE(allionids)
     ALLOCATE(clust_avg(1),stat = AllocateStatus)
     DEALLOCATE(clust_avg)

  END IF
     
END SUBROUTINE SORT_ION_CION_ARR
  
!--------------------------------------------------------------------

SUBROUTINE SORT_COM_ARR()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER i, j
  INTEGER, DIMENSION(1:ntotmols,2) :: dumcomarr
  REAL, DIMENSION(1:ntotmols) :: dumcom_massarr
  INTEGER :: AllocateStatus
  REAL :: massval
  
  dumcomarr = -1; dumcom_massarr = 0; massval = 0

  print *, comtyp_arr
  
  
  DO i = 1, ntotatoms

     IF(ANY(aidvals(i,3) == comtyp_arr)) THEN
        
        IF (.NOT. ANY(aidvals(i,2) ==  dumcomarr(:,1))) THEN
           ! replace type by com_type which is set to 199
           nmol_totcom = nmol_totcom + 1
           dumcomarr(nmol_totcom,1) = aidvals(i,2)
           dumcomarr(nmol_totcom,2) = com_type
           CALL FIND_MASS_FROM_TYPE(aidvals(i,3),massval)
           dumcom_massarr(nmol_totcom) = dumcom_massarr(nmol_totcom) +&
                & massval
        ELSE
           CALL FIND_MASS_FROM_TYPE(aidvals(i,3),massval)
           dumcom_massarr(nmol_totcom) = dumcom_massarr(nmol_totcom) +&
                & massval
        
        END IF

     END IF
     
  END DO

  PRINT *, "Number of COM points (molecule centers): ",nmol_totcom
  WRITE(logout,*) "Number of COM points (molecule centers): "&
       &,nmol_totcom

  ALLOCATE(com_IDTYP_arr(nmol_totcom,2),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate com_IDTYP_arr"
  ALLOCATE(com_masses(nmol_totcom),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate com_IDTYP_arr"
  
  i = 0

  DO WHILE(dumcomarr(i+1,1) .NE. -1) 
     
     i = i + 1
     com_IDTYP_arr(i,1) = dumcomarr(i,1)
     com_IDTYP_arr(i,2) = dumcomarr(i,2)
     com_masses(i) = dumcom_massarr(i)
     
  END DO
  
  IF(i .NE. nmol_totcom) THEN
     
     PRINT *, i, nmol_totcom
     STOP "Wrong total count in com_IDTYP_arr"
     
  END IF
  
  DO i = 1,nmol_totcom
     
     IF(com_IDTYP_arr(i,1) == -1 .OR. com_IDTYP_arr(i,2) == -1)&
          & THEN
        
        PRINT *, i,com_IDTYP_arr(i,1), com_IDTYP_arr(i,2)
        PRINT *, "Something wrong in assigning com_IDTYP_arr"
        STOP
        
     END IF
     
     IF(com_IDTYP_arr(i,2) .NE. com_type) THEN
        
        PRINT *, i,com_type,com_IDTYP_arr(i,1), com_IDTYP_arr(i,2)
        PRINT *, "Something wrong in com_IDTYP_arr"
        STOP
        
     END IF
     
  END DO
  
  
  OPEN(unit = 93,file="COMlist.txt",action="write",status="repl&
       &ace")
  
  WRITE(93,*) "Reference new-type/count: ", com_type, nmol_totcom
  WRITE(93,*) "#  ","ID  ","New-type", "Total Mass"
  DO i = 1,nmol_totcom
 
     WRITE(93,'(3(I0,1X),F14.8)') i,com_IDTYP_arr(i,1)&
          &,com_IDTYP_arr(i,2),com_masses(i)
     
  END DO
  
  CLOSE(93)
  
END SUBROUTINE SORT_COM_ARR

!--------------------------------------------------------------------

SUBROUTINE SANITY_CHECK_IONTYPES()

  USE PARAMS_GMX
  IMPLICIT NONE

  IF(ion_dynflag .OR. catan_neighcalc_flag) THEN

     IF(iontype == -1) THEN

        PRINT *, "ion type undefined for neigh or diff calculation"
        STOP

     END IF

  END IF

  IF(cion_dynflag .OR. catan_neighcalc_flag) THEN

     IF(c_iontype == -1) THEN

        PRINT *, "counter-ion type undefined for neigh or diff calculation"
        STOP

     END IF

  END IF

END SUBROUTINE SANITY_CHECK_IONTYPES

!--------------------------------------------------------------------

SUBROUTINE ANALYZE_TRAJECTORYFILE()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: aid,ierr,atchk,atype,jumpfr,jout
  INTEGER :: at_cnt,molid,kdummy
  REAL :: xlo,xhi,ylo,yhi,zlo,zhi
  CHARACTER(LEN=256):: headline
  LOGICAL :: ok_t, ok_step
  CHARACTER(LEN=max_char):: molname,aname

  INQUIRE(file=traj_fname, exist=ierr)
  IF (ierr == 0) THEN
     WRITE(*,*) 'ERROR: file not found: ', trim(traj_fname)
     STOP
  END IF

  OPEN(unit = trajread,file =traj_fname,action="read",status="old"&
       &,iostat=ierr)

  IF(ierr /= 0) STOP "trajectory file not found"

  PRINT *, "Trajectory file used: ",trim(adjustl(traj_fname))
  WRITE(logout,*) "Trajectory file used: "&
       &,trim(adjustl(traj_fname))
  
  PRINT *, "Analyzing trajectory file..."
  PRINT *, "Beginning to analyze ", nframes, " frames.."
  WRITE(logout,*) "Beginning to analyze ", nframes, " frames.."
  CALL STRUCT_INIT()

  nfrcntr = 0
  
  DO 

     READ(trajread,'(A)') headline
     READ(trajread,*) atchk

     CALL GET_VALUE_REAL(trim(adjustl(headline)),'t',act_time,ok_t)
     CALL GET_VALUE_INT(trim(adjustl(headline)),'step',timestep&
          &,ok_step)

     CALL PRINT_TRAJ_STATS_ERR(headline,ok_t,ok_step,atchk)

     IF(act_time .LT. start_time) THEN

        DO at_cnt = 1,atchk+1 ! +1 -> For the box-line
           READ(trajread,*)
        END DO
        CYCLE
        
     END IF
     
     IF(nfrcntr == 0) THEN
        PRINT *, "Starting time: ", act_time
        WRITE(logout,*) "Starting time: ", act_time
     END IF

     nfrcntr = nfrcntr + 1
     IF(nfrcntr .GT. nframes) EXIT
     tarr_gmx(nfrcntr) = act_time

     IF(mod(nfrcntr,100) == 0) PRINT *, "Analyzed ", nfrcntr, " frames&
          &; Current time (ps): ", act_time
     
     DO at_cnt = 1,atchk
        
        READ(trajread,'(i5,2a5,i5,3f8.3,3f8.4)') molid,molname,aname&
             &,aid,rxyz_gmx(aid,1),rxyz_gmx(aid,2),rxyz_gmx(aid,3)

        CALL MAP_ANAME_TO_ATYPE(aname,atype,kdummy,act_time)
        
        IF(atype .NE. aidvals(aid,3)) THEN
           
           PRINT *, "Incorrect atom ids"
           PRINT *, timestep, act_time
           PRINT *, at_cnt,trim(adjustl(aname)),atype,molid,aidvals(aid&
                &,3)
           STOP
           
        END IF

        ! Store to time-dependent array for dynamics
        IF(ion_dynflag == 1 .AND. atype == iontype) THEN
                 
           CALL MAP_ION_CION_REFTYPE(aid,atype,jout)
           itrx_gmx(jout,nfrcntr) = rxyz_gmx(aid,1)
           itry_gmx(jout,nfrcntr) = rxyz_gmx(aid,2)
           itrz_gmx(jout,nfrcntr) = rxyz_gmx(aid,3)

        ELSEIF(cion_dynflag == 1 .AND. atype == c_iontype) THEN

           CALL MAP_ION_CION_REFTYPE(aid,atype,jout)
           ctrx_gmx(jout,nfrcntr) = rxyz_gmx(aid,1)
           ctry_gmx(jout,nfrcntr) = rxyz_gmx(aid,2)
           ctrz_gmx(jout,nfrcntr) = rxyz_gmx(aid,3)

        END IF
        
     END DO

     READ(trajread,*) box_xl, box_yl, box_zl
     boxx_arr(nfrcntr) = box_xl
     boxy_arr(nfrcntr) = box_yl
     boxz_arr(nfrcntr) = box_zl

     IF(comflag) CALL COMPUTE_COM_ARR(nfrcntr)
     IF(nfrcntr == 1) PRINT *, "Beginning statics analysis..."

     CALL STRUCT_MAIN(nfrcntr)

     DO jumpfr = 1,freqfr
        
        READ(trajread,*)
        READ(trajread,*) atchk

        DO at_cnt = 1,atchk+1

           READ(trajread,*) 

        END DO

     END DO

  END DO

  CLOSE(trajread)

  PRINT *, "Trajectory read completed .."
  PRINT *, "Last frame analyzed ..", tarr_gmx(nfrcntr-1)
  PRINT *, "Total frames analyzed ..", nfrcntr-1
  WRITE(logout,*) "Total frames analyzed ..", nfrcntr-1
  WRITE(logout,*) "Last frame analyzed ..", tarr_gmx(nfrcntr-1)
    
  PRINT *, "Beginning dynamical analysis..."
  CALL DYNAMICS_MAIN()

END SUBROUTINE ANALYZE_TRAJECTORYFILE

!--------------------------------------------------------------------

SUBROUTINE PRINT_TRAJ_STATS_ERR(headline,ok_t,ok_step,atchk)

  USE PARAMS_GMX
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN):: headline
  LOGICAL, INTENT(IN) :: ok_t, ok_step
  INTEGER, INTENT(IN) :: atchk

  IF(ok_t == .FALSE. .OR. ok_step == .FALSE.) THEN
     
     PRINT *, "Error in reading trajectory .."
     PRINT *, trim(adjustl(headline))
     STOP
     
  END IF

  IF(atchk .GT. ntotatoms) THEN
     
     PRINT *, "More atoms found in trajectory than datafile.."
     PRINT *, trim(adjustl(headline))
     PRINT *, atchk, ntotatoms
     STOP
     
  END IF

END SUBROUTINE PRINT_TRAJ_STATS_ERR

!--------------------------------------------------------------------

SUBROUTINE MAP_ION_CION_REFTYPE(jin,atype,jout)
! Maps atomid into the corresponding place in array
  USE PARAMS_GMX

  IMPLICIT NONE

  INTEGER :: i
  INTEGER, INTENT(IN):: jin,atype
  INTEGER, INTENT(OUT) :: jout

  jout = -1

  IF(atype == iontype) THEN

     DO i = 1,ioncnt
        
        IF(jin == ion_IDTYP_arr(i,1)) THEN

           jout = i

           EXIT

        END IF

     END DO

  ELSEIF(atype == c_iontype) THEN

     DO i = 1,c_ioncnt
        
        IF(jin == countion_IDTYP_arr(i,1)) THEN

           jout = i

           EXIT

        END IF

     END DO

  END IF
  
  IF(jout == -1) THEN
     
     PRINT *, jin, atype
     STOP "Could not find a match"

  END IF


END SUBROUTINE MAP_ION_CION_REFTYPE

!--------------------------------------------------------------------

SUBROUTINE MAP_COM_REFTYPE(ref_aid, ref_molid, ref_atype, jout)

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i
  INTEGER, INTENT(IN):: ref_aid, ref_molid, ref_atype
  INTEGER, INTENT(OUT) :: jout

  DO i = 1,nmol_totcom
     
     IF(ref_molid == com_IDTYP_arr(i,1)) THEN
        
        jout = i

        ! Sanity check
        IF(.NOT. ANY(ref_atype == comtyp_arr)) THEN
           
           PRINT *, "ERROR: Unknown atom id in COM molecule "
           PRINT *, ref_aid, ref_atype, ref_molid, comtyp_arr
           STOP
           
        END IF

        EXIT
        
     END IF
     
  END DO

  IF(jout == -1) THEN
       
     PRINT *, ref_aid, ref_atype
     STOP "Could not find a match"

  END IF


END SUBROUTINE MAP_COM_REFTYPE

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_COM_ARR(tcntr)

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: ii, jout
  INTEGER, INTENT(IN) :: tcntr
  REAL :: massval
  
!$OMP PARALLEL DO
  DO ii = 1, nmol_totcom
     comxyz_gmx(ii,1) = 0.0
     comxyz_gmx(ii,2) = 0.0
     comxyz_gmx(ii,3) = 0.0
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(ii,jout,massval)
  DO ii = 1, ntotatoms
     
     IF(ANY(aidvals(ii,3) == comtyp_arr)) THEN 
        
        CALL MAP_COM_REFTYPE(aidvals(ii,1),aidvals(ii,2),aidvals(ii&
             &,3),jout)
        CALL FIND_MASS_FROM_TYPE(aidvals(ii,3),massval)

        comxyz_gmx(jout,1) = comxyz_gmx(jout,1) + massval*rxyz_gmx(ii,1)
        comxyz_gmx(jout,2) = comxyz_gmx(jout,2) + massval*rxyz_gmx(ii,2)
        comxyz_gmx(jout,3) = comxyz_gmx(jout,3) + massval*rxyz_gmx(ii,3)

     END IF

  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
  DO ii = 1,nmol_totcom

     comxyz_gmx(ii,1) = comxyz_gmx(ii,1)/REAL(com_masses(ii))
     comxyz_gmx(ii,2) = comxyz_gmx(ii,2)/REAL(com_masses(ii))
     comxyz_gmx(ii,3) = comxyz_gmx(ii,3)/REAL(com_masses(ii))

  END DO
!$OMP END PARALLEL DO

  IF(com_dynflag) CALL GENERATE_COM_TIME_ARR(nfrcntr)
  
END SUBROUTINE COMPUTE_COM_ARR

!--------------------------------------------------------------------

SUBROUTINE GENERATE_COM_TIME_ARR(tcntr)

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: tcntr
  INTEGER :: imol

!$OMP PARALLEL DO PRIVATE(imol)
  DO imol = 1, nmol_totcom
     
     comtx_gmx(imol,tcntr) = comxyz_gmx(imol,1)
     comty_gmx(imol,tcntr) = comxyz_gmx(imol,2)
     comtz_gmx(imol,tcntr) = comxyz_gmx(imol,3)
     
  END DO
!$OMP END PARALLEL DO

END SUBROUTINE GENERATE_COM_TIME_ARR

!--------------------------------------------------------------------

SUBROUTINE STRUCT_INIT()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,j,t1,t2,norm,acnt,fcnt,a1id,molid,flagch,flagpr,jmax
  INTEGER :: AllocateStatus

  IF(rdfcalc_flag) THEN

     rdfarray = 0.0
     rbinval = rdomcut/REAL(rmaxbin)

     DO i = 1, npairs

        t1 = 0; t2 = 0
        
        IF(pairs_rdf(i,1) /= com_type .AND. pairs_rdf(i,2) /=&
             & com_type) THEN

           DO j = 1,ntotatoms

              IF(aidvals(j,3) == pairs_rdf(i,1)) t1 = t1+1
              IF(aidvals(j,3) == pairs_rdf(i,2)) t2 = t2+1
              
           END DO

           IF(pairs_rdf(i,1) == pairs_rdf(i,2)) THEN
              pairs_rdf(i,3) = t1*(t1-1) !g_AA(r)
           ELSE
              pairs_rdf(i,3) = t1*t2 !g_AB(r)
           END IF

        ELSEIF(pairs_rdf(i,1) == pairs_rdf(i,2) .AND. pairs_rdf(i,1) &
             &== com_type) THEN

           pairs_rdf(i,3) = nmol_totcom*(nmol_totcom-1) !g_COM-COM(r)
           
        ELSEIF(pairs_rdf(i,1) == com_type .AND. pairs_rdf(i,2) &
             & /= com_type) THEN

           DO j = 1,ntotatoms

              IF(aidvals(j,3) == pairs_rdf(i,2)) t1 = t1+1
              
           END DO
           
           pairs_rdf(i,3) = t1*nmol_totcom !g_COM-A(r)
           
        ELSEIF(pairs_rdf(i,2) == com_type .AND. pairs_rdf(i,1) &
             & /= com_type) THEN

           DO j = 1,ntotatoms

              IF(aidvals(j,3) == pairs_rdf(i,1)) t1 = t1+1
              
           END DO

           pairs_rdf(i,3) = t1*nmol_totcom !g_A-COM(r)

        ELSE

           PRINT *, "ERROR: Unknown types for RDF calcalations ..."
           PRINT *, pairs_rdf(:,1), pairs_rdf(:,2)
           STOP

        END IF

     END DO

  END IF
     

END SUBROUTINE STRUCT_INIT

!--------------------------------------------------------------------


SUBROUTINE STRUCT_MAIN(tval)

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER, INTENT(IN):: tval
  INTEGER :: t1, t2
  INTEGER :: clock_rate, clock_max
  CHARACTER(100) :: fname_pref
  INTEGER :: AllocateStatus
  
  IF(rdfcalc_flag) THEN

     IF(tval == 1) THEN

        PRINT *, "Checking RDF calculations ..."
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL COMPUTE_RDF(tval)
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for RDF analysis: ',REAL(t2&
             &-t1)/REAL(clock_rate), ' seconds'

     ELSEIF (mod(tval-1,rdffreq)==0) THEN

        CALL COMPUTE_RDF(tval)

     END IF

  END IF

  IF(catan_neighcalc_flag) THEN

     IF(tval == 1) THEN
        
        PRINT *, "Checking Cation-Anion Neighbor calculations ..."
        cat_an_neighavg = 0.0; an_cat_neighavg=0.0
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL CAT_AN_NEIGHS()
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for neighbor analysis: ',REAL(t2&
             &-t1)/REAL(clock_rate), ' seconds'

     ELSEIF(mod(tval,neighfreq) == 0) THEN

        CALL CAT_AN_NEIGHS()

     END IF

  END IF

  IF(clust_calc_flag) THEN

     IF(tval == 1) THEN

        PRINT *, "Checking Binary cluster calculations ..."
        IF(clust_time_flag) THEN
           
           WRITE(fname_pref,'(A11,I0,A4)') "clusttime_",c_iontype,'.tx&
                &t'       
           dum_fname  = trim(adjustl(fname_pref))
           OPEN(unit = clustwrite,file=dum_fname,action="write"&
                &,status="replace")
        END IF

        clust_avg = 0
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL BINARY_CLUSTER_ANALYSIS(tval)
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for cluster analysis= ',REAL(t2&
             &-t1)/REAL(clock_rate), ' seconds'
     ELSE

        CALL BINARY_CLUSTER_ANALYSIS(tval)

     END IF

  END IF

  IF(multclust_calc_flag) THEN

     IF(tval == 1) THEN
       
        PRINT *, "Checking multi-cluster calculations ..."
        CALL SETUP_MULTITYPE_CLUSTER_ARR()
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL MULTITYPE_CLUSTER_ANALYSIS(tval)
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for cluster analysis= ',REAL(t2&
             &-t1)/REAL(clock_rate), ' seconds'
     ELSE

        CALL MULTITYPE_CLUSTER_ANALYSIS(tval)

     END IF

  ELSE

     ALLOCATE(multionids(1,1),stat = AllocateStatus)
     DEALLOCATE(multionids)
     ALLOCATE(spec_avg(1),stat = AllocateStatus)
     DEALLOCATE(spec_avg)

  END IF

END SUBROUTINE STRUCT_MAIN

!--------------------------------------------------------------------

SUBROUTINE SETUP_MULTITYPE_CLUSTER_ARR()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,j, AllocateStatus, center_cnt
  
  CALL COUNT_TYPES_FOR_MULTITYPE_CLUSTER()
  maxsize_species = product(mclust_type_arr(:,2))
  
  ALLOCATE(multionids(totmult_centers,2),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate multionids"
  ALLOCATE(spec_avg(maxsize_species),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate spec_avg"
  
  multionids = 0
  center_cnt = 0 ! counter for multionids
  spec_avg   = 0.0
  
  DO i = 1, ntotatoms
     
     DO j = 1,nclust_types
        
        IF(aidvals(i,3) == mclust_type_arr(j,1)) THEN
           
           center_cnt = center_cnt + 1
           multionids(center_cnt,1) = aidvals(i,1)
           multionids(center_cnt,2) = aidvals(i,3)
           
           EXIT
           
        END IF
        
     END DO
     
  END DO

  IF(center_cnt .NE. totmult_centers) THEN
     
     PRINT *, "Unequal number in counting centers for mult-clust"
     PRINT *, center_cnt, totmult_centers
     STOP
     
  END IF

END SUBROUTINE SETUP_MULTITYPE_CLUSTER_ARR

!--------------------------------------------------------------------
  
SUBROUTINE COUNT_TYPES_FOR_MULTITYPE_CLUSTER()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,j

  totmult_centers = 0
  
  DO i = 1,ntotatoms

     DO j = 1,nclust_types

        IF(aidvals(i,3) == mclust_type_arr(j,1)) THEN

           mclust_type_arr(j,2) = mclust_type_arr(j,2) + 1
           totmult_centers = totmult_centers + 1
           EXIT

        END IF

     END DO

  END DO

END SUBROUTINE COUNT_TYPES_FOR_MULTITYPE_CLUSTER

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RDF(iframe)

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,a1type,a2type,ibin,a1id,a2id,paircnt,AllocateStatus
  REAL :: rxval,ryval,rzval,rval
  INTEGER :: a1ref,a2ref
  INTEGER,ALLOCATABLE, DIMENSION(:,:) :: dumrdfarray

  rvolval = box_xl*box_yl*box_zl
  rvolavg = rvolavg + rvolval

  ALLOCATE(dumrdfarray(0:rmaxbin-1,npairs),stat=AllocateStatus)
  IF(AllocateStatus/=0) STOP "dumrdfarray not allocated"
  dumrdfarray = 0

  
!$OMP PARALLEL 
!$OMP DO PRIVATE(i,j,a1type,a2type,a1id,a2id,rval,rxval,ryval,rzval&
!$OMP& ,ibin,paircnt,a1ref,a2ref) REDUCTION(+:dumrdfarray)
  DO paircnt = 1,npairs

     a1ref = pairs_rdf(paircnt,1); a2ref = pairs_rdf(paircnt,2)

     ! Both do not correspond to COM - type
     IF(a1ref /= com_type .AND. a2ref /= com_type) THEN

        DO i = 1,ntotatoms
           
           a1id   = aidvals(i,1)
           a1type = aidvals(i,3)
        
           DO j = 1,ntotatoms
              
              a2id   = aidvals(j,1)
              a2type = aidvals(j,3)
           
              ! Remove identical IDs when computing g_AA(r)
              IF(a1id == a2id .AND. a1ref == a2ref) CYCLE
           
              
              IF(a1type == a1ref .AND. a2type == a2ref) THEN
              
                 rxval = rxyz_gmx(a1id,1) - rxyz_gmx(a2id,1)
                 ryval = rxyz_gmx(a1id,2) - rxyz_gmx(a2id,2)
                 rzval = rxyz_gmx(a1id,3) - rxyz_gmx(a2id,3)
                 
                 rxval = rxval - box_xl*ANINT(rxval/box_xl)
                 ryval = ryval - box_yl*ANINT(ryval/box_yl)
                 rzval = rzval - box_zl*ANINT(rzval/box_zl)
              
                 rval = sqrt(rxval**2 + ryval**2 + rzval**2)
                 ibin = FLOOR(rval/rbinval)
              
              
                 IF(ibin .LT. rmaxbin) THEN
                 
                    dumrdfarray(ibin,paircnt) = dumrdfarray(ibin&
                         &,paircnt) + 1
                    
                 END IF
                 
              END IF
           
           END DO

        END DO

        
     ELSEIF(a1ref == com_type .AND. a2ref /= com_type) THEN
        !a1type == com_type

        DO i = 1,nmol_totcom
           
           a1id   = i
           a1type = com_type
           
           DO j = 1,ntotatoms
              
              a2id   = aidvals(j,1)
              a2type = aidvals(j,3)
              
              IF(a2type == a2ref) THEN
              
                 rxval = comxyz_gmx(a1id,1) - rxyz_gmx(a2id,1)
                 ryval = comxyz_gmx(a1id,2) - rxyz_gmx(a2id,2)
                 rzval = comxyz_gmx(a1id,3) - rxyz_gmx(a2id,3)
                 
                 rxval = rxval - box_xl*ANINT(rxval/box_xl)
                 ryval = ryval - box_yl*ANINT(ryval/box_yl)
                 rzval = rzval - box_zl*ANINT(rzval/box_zl)
              
                 rval = sqrt(rxval**2 + ryval**2 + rzval**2)
                 ibin = FLOOR(rval/rbinval)
              
              
                 IF(ibin .LT. rmaxbin) THEN
                 
                    dumrdfarray(ibin,paircnt) = dumrdfarray(ibin&
                         &,paircnt) + 1
                    
                 END IF
                 
              END IF
           
           END DO

        END DO

     ELSEIF(a1ref /= com_type .AND. a2ref == com_type) THEN
        !a2ref = com_type

        DO i = 1,ntotatoms
           
           a1id     = aidvals(i,1)
           a1type   = aidvals(i,3)
           
           DO j = 1,nmol_totcom
              
              a2id  = j
              a2type = com_type
              
              IF(a1type == a1ref) THEN
              
                 rxval = rxyz_gmx(a1id,1) - comxyz_gmx(a2id,1)
                 ryval = rxyz_gmx(a1id,2) - comxyz_gmx(a2id,2)
                 rzval = rxyz_gmx(a1id,3) - comxyz_gmx(a2id,3)
                 
                 rxval = rxval - box_xl*ANINT(rxval/box_xl)
                 ryval = ryval - box_yl*ANINT(ryval/box_yl)
                 rzval = rzval - box_zl*ANINT(rzval/box_zl)
              
                 rval = sqrt(rxval**2 + ryval**2 + rzval**2)
                 ibin = FLOOR(rval/rbinval)
              
              
                 IF(ibin .LT. rmaxbin) THEN
                 
                    dumrdfarray(ibin,paircnt) = dumrdfarray(ibin&
                         &,paircnt) + 1
                    
                 END IF
                 
              END IF
              
           END DO
           
        END DO

     ELSEIF(a1ref == com_type .AND. a1ref == com_type) THEN
        ! Both a1ref and a2ref are COM_types
        
        DO i = 1,nmol_totcom
           
           a1id   = i
           a1type = com_type
           DO j = 1,nmol_totcom
              
              a2id   = j
              a2type = com_type
              ! Remove identical IDs when computing g_COM-COM(r)
              IF(a1id == a2id) CYCLE
                         
              rxval = comxyz_gmx(a1id,1) - comxyz_gmx(a2id,1)
              ryval = comxyz_gmx(a1id,2) - comxyz_gmx(a2id,2)
              rzval = comxyz_gmx(a1id,3) - comxyz_gmx(a2id,3)
              
              rxval = rxval - box_xl*ANINT(rxval/box_xl)
              ryval = ryval - box_yl*ANINT(ryval/box_yl)
              rzval = rzval - box_zl*ANINT(rzval/box_zl)
              
              rval = sqrt(rxval**2 + ryval**2 + rzval**2)
              ibin = FLOOR(rval/rbinval)
              
              
              IF(ibin .LT. rmaxbin) THEN
                 
                 dumrdfarray(ibin,paircnt) = dumrdfarray(ibin&
                      &,paircnt) + 1
                 
              END IF
              
           END DO
           
        END DO

     ELSE
        
        PRINT *, "ERROR: Unknown types in pairs_rdf "
        PRINT *, pairs_rdf(paircnt,1), pairs_rdf(paircnt,2)
        STOP
        
     END IF
             
  END DO
!$OMP END DO

!$OMP DO PRIVATE(i,j)
  DO j = 1,npairs

     DO i = 0,rmaxbin-1

        rdfarray(i,j) = rdfarray(i,j) + REAL(dumrdfarray(i,j))&
             &*rvolval/(REAL(pairs_rdf(j,3)))

     END DO

  END DO
!$OMP END DO

!$OMP END PARALLEL

  DEALLOCATE(dumrdfarray)

END SUBROUTINE COMPUTE_RDF

!--------------------------------------------------------------------

SUBROUTINE CAT_AN_NEIGHS()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,j,a1id,a2id,neigh_cnt,tid
  INTEGER,DIMENSION(1:maxneighsize,0:nproc-1) :: cat_an_neigh_inst&
       &,an_cat_neigh_inst
  REAL :: rxval, ryval, rzval, rval

  cat_an_neigh_inst = 0; an_cat_neigh_inst = 0

!$OMP PARALLEL PRIVATE(i,j,a1id,a2id,rxval,ryval,rzval,rval,neigh_cnt,tid)
!$OMP DO
  DO i = 1,ioncnt

     neigh_cnt = 0
     a1id = ion_IDTYP_arr(i,1)
     tid = OMP_GET_THREAD_NUM()

     DO j = 1,c_ioncnt

        a2id = countion_IDTYP_arr(j,1)

        rxval = rxyz_gmx(a1id,1) - rxyz_gmx(a2id,1)
        ryval = rxyz_gmx(a1id,2) - rxyz_gmx(a2id,2)
        rzval = rxyz_gmx(a1id,3) - rxyz_gmx(a2id,3)

        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)

        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        
        IF(rval .LT. rneigh_cut) THEN

           neigh_cnt = neigh_cnt + 1

        END IF

     END DO

     IF(neigh_cnt + 1 .GT. maxneighsize) THEN

        PRINT *, "Neighbor count exceeded max size"
        PRINT *, neigh_cnt, maxneighsize
        STOP

     END IF

     cat_an_neigh_inst(neigh_cnt+1,tid) = cat_an_neigh_inst(neigh_cnt&
          &+1,tid) + 1

  END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(i,j,a1id,a2id,rxval,ryval,rzval,rval,neigh_cnt,tid)
!$OMP DO
  DO i = 1,c_ioncnt

     neigh_cnt = 0
     a1id = countion_IDTYP_arr(i,1)
     tid = OMP_GET_THREAD_NUM()

     DO j = 1,ioncnt

        a2id = ion_IDTYP_arr(j,1)

        rxval = rxyz_gmx(a1id,1) - rxyz_gmx(a2id,1)
        ryval = rxyz_gmx(a1id,2) - rxyz_gmx(a2id,2)
        rzval = rxyz_gmx(a1id,3) - rxyz_gmx(a2id,3)

        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)

        rval = sqrt(rxval**2 + ryval**2 + rzval**2)

        IF(rval .LT. rneigh_cut) THEN

           neigh_cnt = neigh_cnt + 1

        END IF

     END DO

     IF(neigh_cnt + 1 .GT. maxneighsize) THEN

        PRINT *, "Neighbor count exceeded max size"
        PRINT *, neigh_cnt, maxneighsize
        STOP

     END IF

     an_cat_neigh_inst(neigh_cnt+1,tid) = an_cat_neigh_inst(neigh_cnt&
          &+1,tid) + 1

  END DO
!$OMP END DO

!$OMP DO 
  DO  i = 1,maxneighsize
     DO j = 0,nproc-1
        cat_an_neighavg(i) = cat_an_neighavg(i) + cat_an_neigh_inst(i&
             &,j)
        an_cat_neighavg(i) = an_cat_neighavg(i) + an_cat_neigh_inst(i&
             &,j)
     END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE CAT_AN_NEIGHS

!--------------------------------------------------------------------

SUBROUTINE SOLVENT_SEPARATED_ION_PAIRS()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,j,a1id,a2id,neigh_cnt,tid
  INTEGER,DIMENSION(1:maxneighsize,0:nproc-1) :: cat_an_neigh_inst&
       &,an_cat_neigh_inst
  REAL :: rxval, ryval, rzval, rval

  cat_an_neigh_inst = 0; an_cat_neigh_inst = 0

!$OMP PARALLEL PRIVATE(i,j,a1id,a2id,rxval,ryval,rzval,rval,neigh_cnt,tid)
!$OMP DO
  DO i = 1,ioncnt

     neigh_cnt = 0
     a1id = ion_IDTYP_arr(i,1)
     tid = OMP_GET_THREAD_NUM()

     DO j = 1,c_ioncnt

        a2id = countion_IDTYP_arr(j,1)

        rxval = rxyz_gmx(a1id,1) - rxyz_gmx(a2id,1)
        ryval = rxyz_gmx(a1id,2) - rxyz_gmx(a2id,2)
        rzval = rxyz_gmx(a1id,3) - rxyz_gmx(a2id,3)

        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)

        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        
        IF(rval .LT. rneigh_cut) THEN

           neigh_cnt = neigh_cnt + 1

        END IF

     END DO

     IF(neigh_cnt + 1 .GT. maxneighsize) THEN

        PRINT *, "Neighbor count exceeded max size"
        PRINT *, neigh_cnt, maxneighsize
        STOP

     END IF

     cat_an_neigh_inst(neigh_cnt+1,tid) = cat_an_neigh_inst(neigh_cnt&
          &+1,tid) + 1

  END DO
!$OMP END DO
!$OMP END PARALLEL
  

END SUBROUTINE SOLVENT_SEPARATED_ION_PAIRS

!--------------------------------------------------------------------

SUBROUTINE BINARY_CLUSTER_ANALYSIS(frnum)

  USE PARAMS_GMX
  IMPLICIT NONE

!Ref Sevick et.al ., J Chem Phys 88 (2)
!!$  INTEGER, DIMENSION(ntotion_centers,ntotion_centers) :: all_direct,&
!!$       & all_neigh
!!$  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: catan_direct
  
  INTEGER :: i,j,k,a2ptr,a1id,a2id,itype,jtype,jptr,idum,jflag,jcnt&
       &,iflag,jtot,jind,jprev
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: all_direct,all_neigh
  INTEGER, DIMENSION(1:ntotion_centers) :: union_all,scnt,all_linked
  REAL :: rxval, ryval, rzval, rval
  INTEGER, INTENT(IN) :: frnum
  INTEGER :: AllocateStatus

  ALLOCATE(all_direct(ntotion_centers,ntotion_centers),stat =&
       & AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate all_direct"

  ALLOCATE(all_neigh(ntotion_centers,ntotion_centers),stat =&
       & AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate all_neigh"

!$OMP PARALLEL 

!$OMP DO PRIVATE(i,j)
  DO i = 1,ntotion_centers

     scnt(i) = 0; all_linked(i)  = 0
     union_all(i) = -1

     DO j = 1,ntotion_centers

        IF(i .NE. j) THEN
           all_direct(i,j) = 0
!!$           catan_direct(i,j) = 0
        END IF

        IF(i == j) THEN
           all_direct(i,j) = 1
!!$           catan_direct(i,j) = 0
        END IF

        all_neigh(i,j) = 0
        
     END DO

  END DO
!$OMP END DO

!Create Direct connectivity matrix
!all_direct - does not distinguish between Li and P neigh
!catan_direct - neighbors with sequence cat-an-cat-an.. or an-cat-an-cat...

!$OMP DO PRIVATE(i,j,a1id,a2ptr,a2id,rxval,ryval,rzval,rval,itype&
!$OMP& ,jptr,jtype)  
  DO i = 1,ntotion_centers

     a1id = allionids(i,1)
     a2ptr = 1
     itype = aidvals(a1id,3)
     jptr  = 1
     all_neigh(i,i) = a1id

     DO j = 1,ntotion_centers

        a2id = allionids(j,1)

        rxval = rxyz_gmx(a1id,1) - rxyz_gmx(a2id,1)
        ryval = rxyz_gmx(a1id,2) - rxyz_gmx(a2id,2)
        rzval = rxyz_gmx(a1id,3) - rxyz_gmx(a2id,3)

        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)

        rval = sqrt(rxval**2 + ryval**2 + rzval**2)

        IF(rval .LT. rclust_cut .AND. a1id .NE. a2id) THEN

           all_direct(i,j) = 1
           all_neigh(i,j)  = a2id
!           all_neigh(i,a2ptr) = a2id
!           a2ptr = a2ptr + 1

           jtype = aidvals(a2id,3)

!           IF(itype .NE. jtype) THEN

!              catan_direct(i,j) = 1
!              catan_neigh(i,j)  = a2id
!              catan_neigh(i,jptr+1) = a2id
!              itype = jtype
!              jptr  = jptr + 1

!           END IF

        END IF

     END DO

  END DO

!$OMP END DO  

  
!Check for symmetry
  IF(frnum == 1) THEN
!$OMP DO
     DO i = 1,ntotion_centers

        DO j = 1,ntotion_centers

           IF(all_direct(i,j) .NE. all_direct(j,i)) STOP "Unsymmetric&
                & all_direct"

          IF(all_neigh(i,j) .NE. 0) THEN

              IF(all_neigh(i,j) .NE. all_neigh(j,j) .OR. all_neigh(j&
                   &,i) .NE. all_neigh(i,i)) THEN

                 PRINT *, i,j,all_direct(i,j),all_direct(j,i)&
                      &,all_neigh(j,i),all_neigh(i,i)
                 STOP "Unsymmetric neighbor list"

              END IF

           END IF

        END DO

     END DO
!$OMP END DO
  END IF

!$OMP END PARALLEL        

  !Intersection

  DO i = 1,ntotion_centers-1 !Ref row counter

     iflag = 0
     idum  = i

     DO WHILE(iflag == 0 .AND. union_all(i) == -1)

        jflag = 0
        k    = 1 !Column counter
        j    = idum+1 !Other row counter

        DO WHILE(jflag == 0 .AND. k .LE. ntotion_centers)

           IF((all_direct(i,k) == all_direct(j,k)).AND. all_direct(i&
                &,k)== 1) THEN

              jflag = 1
!!$              jprev = 0

              DO jcnt = 1,ntotion_centers


!!$                 IF(all_direct(j,jcnt) == 1) jprev = 1

                 !Replace highest row by union of two rows

                 all_direct(j,jcnt) = all_direct(i,jcnt) .OR.&
                      & all_direct(j,jcnt)

!!$                 IF((all_direct(j,jcnt) == 1 .AND. all_direct(i,jcnt)&
!!$                      &==1) .AND. jprev == 0) THEN
!!$                    
!!$                    all_neigh(j,jcnt) = all_neigh(i,jcnt)
!!$                    jprev = 0 !Other condition is already
!!$                    ! incorporated before
!!$                 END IF
!!$                 
              END DO

              union_all(i) = 1 !One match implies the low ranked row
              ! is present in high ranked row

           ELSE

              k = k + 1

           END IF

        END DO

        IF(union_all(i) == 1) THEN

           iflag = 1

        ELSE

           idum  = idum + 1

        END IF

        IF(idum == ntotion_centers) iflag = 1

     END DO

  END DO

!Count
  jtot = 0
!$OMP PARALLEL PRIVATE(i,j,jind) 
!$OMP DO
  DO i = 1,ntotion_centers

     IF(union_all(i) == -1) THEN

        jind = 0

        DO j = 1,ntotion_centers

           IF(all_direct(i,j) == 1) jind = jind + 1

        END DO

        scnt(jind) = scnt(jind) + 1
        all_linked(i) = jind

     END IF

  END DO
!$OMP END DO

!$OMP DO

  DO i = 1,ntotion_centers

     clust_avg(i) = clust_avg(i) + scnt(i)

  END DO
!$OMP END DO

!$OMP END PARALLEL

  IF(frnum == 1) THEN
     OPEN(unit =90,file ="scnt.txt",action="write",status="replace")
     IF(clust_time_flag) WRITE(clustwrite,'(3(I0,1X),F14.8)')&
          & ntotion_centers, iontype, c_iontype, rclust_cut
  END IF

  jtot = 0

  DO i = 1,ntotion_centers

     IF(frnum == 1) WRITE(90,*) i,scnt(i)
     jtot = jtot + all_linked(i)
     
  END DO

  IF(clust_time_flag) WRITE(clustwrite,*) frnum, scnt(:)

  IF(jtot .NE. ntotion_centers) THEN

     PRINT *, "Sum of ions not equal to total ion centers"
     PRINT *, jtot, ntotion_centers
     STOP

  END IF

  IF(frnum == 1) CLOSE(90)

  IF(frnum == 1) THEN

     OPEN(unit =90,file ="all_neigh.txt",action="write",status="replace")

     DO i = 1,ntotion_centers

        IF(union_all(i) == -1) THEN

           WRITE(90,*) i,all_linked(i)

           DO j = 1,ntotion_centers

              IF(all_direct(i,j) == 1) WRITE(90,*) j,allionids(j,1),&
                   & allionids(j,2),all_direct(i,j)

           END DO

        END IF

     END DO

     CLOSE(90)

  END IF

  DEALLOCATE(all_direct)
  DEALLOCATE(all_neigh)
  
END SUBROUTINE BINARY_CLUSTER_ANALYSIS

!--------------------------------------------------------------------

SUBROUTINE MULTITYPE_CLUSTER_ANALYSIS(frnum)

  USE PARAMS_GMX
  IMPLICIT NONE

!Extending Ref Sevick et.al ., J Chem Phys 88 (2)

  INTEGER :: i,j,k,a2ptr,a1id,a2id,itype,jtype,jptr,idum,jflag,jcnt&
       &,iflag,jtot,jind,jprev,spec_ind,stride,i_index,j_index
  INTEGER, DIMENSION(1:totmult_centers,1:totmult_centers) ::&
       & all_direct,all_neigh
  INTEGER, DIMENSION(1:totmult_centers) :: union_all,scnt,all_linked
  INTEGER, DIMENSION(1:maxsize_species) :: sum_species
  INTEGER, DIMENSION(1:nclust_types) :: sum_atoms
  REAL :: rxval, ryval, rzval, rval, rcut_ij
  INTEGER, INTENT(IN) :: frnum

!$OMP PARALLEL
!$OMP DO PRIVATE(i,j)
  DO i = 1,totmult_centers

     scnt(i) = 0; all_linked(i)  = 0
     union_all(i) = -1

     DO j = 1,totmult_centers
        
        IF(i == j) THEN
           all_direct(i,j) = 1
        ELSE
           all_direct(i,j) = 0
        END IF

        all_neigh(i,j) = 0
        
     END DO

  END DO
!$OMP END DO

!Create Direct connectivity matrix
!all_direct - does not distinguish between different molecules


!$OMP DO PRIVATE(i,j,a1id,a2ptr,a2id,rxval,ryval,rzval,rval,itype&
!$OMP& ,jptr,jtype,rcut_ij,i_index,j_index)  
  DO i = 1,totmult_centers

     a1id  = multionids(i,1)
     a2ptr = 1
     itype = aidvals(a1id,3)
     jptr  = 1
     all_neigh(i,i) = a1id

     DO j = 1,totmult_centers

        a2id = multionids(j,1)
        jtype = aidvals(a2id,3)
        
        rxval = rxyz_gmx(a1id,1) - rxyz_gmx(a2id,1)
        ryval = rxyz_gmx(a1id,2) - rxyz_gmx(a2id,2)
        rzval = rxyz_gmx(a1id,3) - rxyz_gmx(a2id,3)

        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)

        rval = sqrt(rxval**2 + ryval**2 + rzval**2)

        CALL MAP_TYPE_TO_INDEX(itype,i_index)
        CALL MAP_TYPE_TO_INDEX(jtype,j_index)
        rcut_ij = mclust_rcut_arr(i_index,j_index)
        
        IF(rval .LT. rcut_ij .AND. a1id .NE. a2id) THEN

           all_direct(i,j) = 1
           all_neigh(i,j)  = a2id

        END IF

     END DO

  END DO

!$OMP END DO  

  
  
!Check for symmetry
  IF(frnum == 1) THEN
!$OMP DO
     DO i = 1,totmult_centers

        DO j = 1,totmult_centers

           IF(all_direct(i,j) .NE. all_direct(j,i)) THEN

              PRINT *, i, j, all_direct(i,j), all_direct(j,i)
              STOP "Unsymmetric all_direct"

           END IF

           IF(all_neigh(i,j) .NE. 0) THEN

              IF(all_neigh(i,j) .NE. all_neigh(j,j) .OR. all_neigh(j&
                   &,i) .NE. all_neigh(i,i)) THEN

                 PRINT *, i,j,all_direct(i,j),all_direct(j,i)&
                      &,all_neigh(j,i),all_neigh(i,i)
                 STOP "Unsymmetric neighbor list"

              END IF

           END IF

        END DO

     END DO
!$OMP END DO
  END IF

!$OMP END PARALLEL        

  !Intersection

  DO i = 1,totmult_centers-1 !Ref row counter

     iflag = 0
     idum  = i

     DO WHILE(iflag == 0 .AND. union_all(i) == -1)

        jflag = 0
        k    = 1 !Column counter
        j    = idum+1 !Other row counter

        DO WHILE(jflag == 0 .AND. k .LE. totmult_centers)

           IF((all_direct(i,k) == all_direct(j,k)).AND. all_direct(i&
                &,k)== 1) THEN

              jflag = 1

              DO jcnt = 1,totmult_centers

                 !Replace highest row by union of two rows
                 all_direct(j,jcnt) = all_direct(i,jcnt) .OR.&
                      & all_direct(j,jcnt)

              END DO

              union_all(i) = 1 !One match implies the low ranked row
              ! is present in high ranked row

           ELSE

              k = k + 1

           END IF

        END DO

        IF(union_all(i) == 1) THEN

           iflag = 1

        ELSE

           idum  = idum + 1

        END IF

        IF(idum == totmult_centers) iflag = 1

     END DO

  END DO

!Count
  jtot = 0
  sum_species = 0

  

  !** sum_atoms(i)**
  !sum_atoms(i) corresponds to the number of occurences of type "i" in
  !each row of all_direct(i,j) when union_all(i) = -1

  !** sum_species(i) **
  !sum_species(i) is a 1D array obtained by converting the nD
  !sum_atoms(i)

!!$  print *, multionids(:,2)

!!$  print *, "all independent rows"
!!$  do i = 1,totmult_centers
!!$     if (union_all(i) == -1) print *, i
!!$  end do
  
!$OMP PARALLEL 
!$OMP DO PRIVATE(i,j,k,jind,sum_atoms,spec_ind,stride,j_index) &
!$OMP& REDUCTION(+:sum_species) 

  DO i = 1,totmult_centers

     IF(union_all(i) == -1) THEN

        jind = 0
        sum_atoms = 0
        
        DO j = 1,totmult_centers

           IF(all_direct(i,j) == 1) THEN

              jind = jind + 1 ! No-identity preserved

              ! With identity preserved
              ! multionids(j,2): type of j atom in multionids
              CALL MAP_TYPE_TO_INDEX(multionids(j,2),j_index)

              !Sanity check
              IF(j_index > nclust_types) THEN

                 PRINT *, "Unphysical j_index", j_index, i,&
                      & nclust_types
                 STOP

              END IF

              sum_atoms(j_index) = sum_atoms(j_index) + 1
              
           END IF

        END DO

!!$        print *, "all_direct row", all_direct(i,:)
!!$        print *, "sum_atoms", i,jind,multionids(i,2),sum_atoms

        ! Note, spec_ind cannot be 0 since at least the element will
        ! be bonded to itself
        
        ! Need to convert the nD array to 1D array
        ! mclust_type_arr(k,2): amount of type k
        spec_ind = 0; stride = 1
        DO k = 1,nclust_types
           
           spec_ind = spec_ind + sum_atoms(k) * stride
           stride = stride * mclust_type_arr(k,2)
           
        END DO

        IF(spec_ind == 0) THEN

           PRINT *, "Unphysical all_direct matrix", i, sum_atoms,&
                & all_direct(i,:)
           STOP

        END IF
        
        sum_species(spec_ind) = sum_species(spec_ind) + 1

           
        scnt(jind) = scnt(jind) + 1
        all_linked(i) = jind

     END IF

  END DO
!$OMP END DO

!$OMP DO PRIVATE(i)

  DO i = 1,maxsize_species

     spec_avg(i) = spec_avg(i) + sum_species(i)

  END DO
!$OMP END DO

  
!$OMP END PARALLEL

  IF(frnum == 1) THEN
     OPEN(unit =90,file ="scnt.txt",action="write",status="replace")


     OPEN(unit =92,file ="species.txt",action="write",status&
          &="replace")
     
     DO i = 1, maxsize_species

        WRITE(92,*) i, sum_species(i)

     END DO

     CLOSE(92)
     
  END IF

  jtot = 0

  DO i = 1, totmult_centers

     IF(frnum == 1) WRITE(90,*) i,scnt(i)
     jtot = jtot + all_linked(i)
     
  END DO
  
  IF(jtot .NE. totmult_centers) THEN

     PRINT *, "Sum of centers not equal to total molecules"
     PRINT *, jtot, totmult_centers
     STOP

  END IF

  IF(frnum == 1) CLOSE(90)

  IF(frnum == 1) THEN

     OPEN(unit =90,file ="all_neigh.txt",action="write",status="replace")

     DO i = 1,totmult_centers

        IF(union_all(i) == -1) THEN

           WRITE(90,*) i,all_linked(i)

           DO j = 1,totmult_centers

              IF(all_direct(i,j) == 1) WRITE(90,*) i,j,multionids(i&
                   &,2),multionids(j,2)

           END DO

        END IF

     END DO

     CLOSE(90)

  END IF

END SUBROUTINE MULTITYPE_CLUSTER_ANALYSIS

!--------------------------------------------------------------------

SUBROUTINE ALLOUTPUTS()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,ierr

  IF (nframes == 0 .AND. nfrcntr .GT. 0) nframes = nfrcntr
  PRINT *, "Number of frames from start to end: ", nframes/(freqfr+1)
  PRINT *, "Frequency of Frames: ", freqfr + 1
  PRINT *, "Total number of Frames analyzed: ", nfrcntr

  WRITE(logout,*) "Number of frames from start to end: ", nframes&
       &/(freqfr+1)
  WRITE(logout,*) "Frequency of Frames: ", freqfr+1
  WRITE(logout,*) "Total number of Frames analyzed: ", nfrcntr

  IF(rdfcalc_flag) THEN
     PRINT *, "Writing RDFs .."
     CALL OUTPUT_ALLRDF()
  END IF

  IF(catan_neighcalc_flag) THEN
     PRINT *, "Writing neighbors .."
     CALL OUTPUT_ALLNEIGHBORS()
  END IF

  IF(clust_calc_flag) THEN

     PRINT *, "Writing binary-cluster outputs .."
     CALL OUTPUT_BINARY_CLUSTERS()

  END IF

  IF(multclust_calc_flag) THEN

     PRINT *, "Writing COM-cluster outputs .."
     CALL OUTPUT_MULTI_CLUSTERS()

  END IF

  
END SUBROUTINE ALLOUTPUTS

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_BINARY_CLUSTERS()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,ierr
  CHARACTER(100) :: fname_pref
  
  WRITE(fname_pref,'(A11,I0,I0,A1)') "clust_",iontype,c_iontype,'_'
  dum_fname = trim(adjustl(fname_pref))//trim(adjustl(traj_fname))

  OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
       &,status="replace",iostat=ierr)
  
  IF(ierr /= 0) PRINT *, "Unknown clust_filename"
  
  DO i = 1,ntotion_centers
     
     WRITE(dumwrite,'(I0,1X,F14.8,1X)') i, REAL(clust_avg(i))&
          &/REAL(nframes)
     
  END DO
  CLOSE(dumwrite)

  IF(clust_time_flag) CLOSE(clustwrite)
  
END SUBROUTINE OUTPUT_BINARY_CLUSTERS

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_MULTI_CLUSTERS()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,ierr,k,stride
  INTEGER :: atom_index,index_remain
  
  dum_fname = "multiclust_"//trim(adjustl(traj_fname))//".dat"

  OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
       &,status="replace",iostat=ierr)
  
  IF(ierr /= 0) PRINT *, "Unknown multiclust_filename"

  ! Unflatten and write only the non-zero species
  
  DO i = 1,maxsize_species

     IF(spec_avg(i) .NE. 0) THEN

        index_remain = i; stride = 1

        WRITE(dumwrite,'(I0,1X)',advance="no") i
        
        DO k = 1,nclust_types

           IF(k>1) stride = stride*mclust_type_arr(k-1,2)
           atom_index = MOD(index_remain/stride,mclust_type_arr(k,2))

           WRITE(dumwrite,'(I0,1X)',advance="no") atom_index

        END DO
             
        WRITE(dumwrite,'(F14.8,1X)') REAL(spec_avg(i))/REAL(nframes)

     END IF
        
  END DO
  CLOSE(dumwrite)
  
END SUBROUTINE OUTPUT_MULTI_CLUSTERS

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLNEIGHBORS()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,frnorm,ierr
  CHARACTER(100) :: fname_pref
  REAL :: totcat_an_neigh,totan_cat_neigh

  IF(neighfreq == 1) frnorm = nframes
  IF(neighfreq .NE. 1) frnorm = nframes/neighfreq + 1

  totcat_an_neigh = 0.0; totan_cat_neigh = 0.0

  IF(catan_neighcalc_flag) THEN
!$OMP PARALLEL DO REDUCTION(+:totcat_an_neigh,totan_cat_neigh) PRIVATE(i)

     DO i = 1,maxneighsize

        totcat_an_neigh = totcat_an_neigh + REAL(cat_an_neighavg(i))
        totan_cat_neigh = totan_cat_neigh + REAL(an_cat_neighavg(i))

     END DO

!$OMP END PARALLEL DO

     WRITE(fname_pref,'(A11,I0,I0,A1)') "catanneigh_",iontype&
          &,c_iontype,'_'
     dum_fname = trim(adjustl(fname_pref))//trim(adjustl(traj_fname))
     OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
          &,status="replace")

     DO i = 1,maxneighsize

        WRITE(dumwrite,'(I0,1X,4(F14.8,1X))') i-1,&
             & REAL(cat_an_neighavg(i))/REAL(frnorm),100.0&
             &*REAL(cat_an_neighavg(i))/totcat_an_neigh&
             &,REAL(an_cat_neighavg(i))/REAL(frnorm),100.0&
             &*REAL(an_cat_neighavg(i))/totan_cat_neigh

     END DO

     CLOSE(dumwrite)
     
  END IF

END SUBROUTINE OUTPUT_ALLNEIGHBORS

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLRDF()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,j,ierr
  REAL, PARAMETER :: vconst = 4.0*pival/3.0
  REAL :: rlower,rupper,nideal,rdffrnorm,acrnorm

  IF(rdfcalc_flag) THEN

     rdffrnorm = INT(nfrcntr/rdffreq)
     rvolavg = rvolavg/REAL(rdffrnorm)
     PRINT *, "Average volume of box", rvolavg

     IF(rdfcalc_flag) THEN
        dum_fname = "rdf_"//trim(adjustl(traj_fname))
        OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
             &,status="replace",iostat=ierr)

        IF(ierr /= 0) THEN
           PRINT *, "Could not open", trim(dum_fname)
        END IF

        WRITE(dumwrite,'(A,8X)',advance="no") "r"

        DO j = 1,npairs

           WRITE(dumwrite,'(I0,A1,I0,8X)',advance="no") pairs_rdf(j&
                &,1),'-',pairs_rdf(j,2)

        END DO

        WRITE(dumwrite,*)

        DO i = 0,rmaxbin-1

           rlower = real(i)*rbinval
           rupper = rlower + rbinval
           nideal = vconst*(rupper**3 - rlower**3)

           WRITE(dumwrite,'(F16.5,2X)',advance="no") 0.5*rbinval&
                &*(REAL(2*i+1))

           DO j = 1,npairs

              WRITE(dumwrite,'(F16.9,1X)',advance="no")rdfarray(i,j)&
               &/(rdffrnorm*nideal)

           END DO

           WRITE(dumwrite,*)

        END DO

        CLOSE(dumwrite)

     END IF

  END IF

END SUBROUTINE OUTPUT_ALLRDF

!--------------------------------------------------------------------

SUBROUTINE DYNAMICS_MAIN()

  USE PARAMS_GMX
  IMPLICIT NONE

  IF(ion_diff) THEN
     PRINT *, "Beginning ion diffusion calculation..."
     CALL DIFF_IONS()
     PRINT *, "Finished ion diffusion calculation..."
  END IF

  IF(cion_diff) THEN
     PRINT *, "Beginning counter-ion diffusion calculation..."
     CALL DIFF_COUNTERIONS()
     PRINT *, "Finished counter-ion diffusion calculation..."
  END IF

  IF(com_diff) THEN
     PRINT *, "Beginning COM diffusion calculation..."
     CALL DIFF_COMGROUP()
     PRINT *, "Finished COM diffusion calculation..."
  END IF

  
  IF(catan_autocfflag) THEN
     PRINT *, "Beginning cat-an residence time calculation..."
     CALL RESIDENCE_TIME_ANCAT()
     PRINT *, "Finished cat-an residence time calculation..."
  END IF

  IF(catCOM_autocfflag) THEN
     PRINT *, "Beginning cat-COM residence time calculation..."
     CALL RESIDENCE_TIME_CATCOM()
     PRINT *, "Finished cat-COM residence time calculation..."
  END IF

  IF(dynfsktflag) THEN
     
     IF(q_targ_min < 0) STOP "q_targ_max should be greater than 0"
     IF(q_targ_max < q_targ_min) STOP "q_targ_max > q_targ_min"
     IF(q_tol > q_bin)  STOP "q_tol should be lesser than q_bin"

     q_targ = q_targ_min

     DO WHILE (q_targ <= q_targ_max)

        PRINT *, "Qtarget: ", q_targ
        
        PRINT *, "Making q-target bounds..."
        CALL MAKE_QTARG_BOUNDS()
     
        PRINT *, "Making q-vectors..."
        CALL MAKE_QVECS()

        IF(q_possN == 0) THEN
           
           PRINT *, "No integer q-vectors found near target |q|"
           PRINT *, q_targ, q_lo, q_hi, q_tol
        ELSE
           
           PRINT *, "Beginning ion self-intermediate scattering..."
           CALL DYNAMIC_SELF_ION_FSKT()
           PRINT *, "Finished ion self-intermediate scattering..."
           
           PRINT *, "Beginning countion self-intermediate scattering..."
           CALL DYNAMIC_SELF_COUNTION_FSKT()
           PRINT *, "Finished countion self-intermediate scattering..."
           
        END IF

        DEALLOCATE(q_nlist) ! deallocate for next q_nlist data
        q_targ = q_targ + q_bin

     END DO
        
  END IF
     

END SUBROUTINE DYNAMICS_MAIN

!--------------------------------------------------------------------

SUBROUTINE DIFF_IONS()

  USE PARAMS_GMX

  IMPLICIT NONE

  INTEGER :: i,j,tinc,ifin,tim,aid,atype,ierr,jout
  REAL    :: rxcm, rycm, rzcm
  REAL, DIMENSION(0:nframes-1) :: gxarr,gyarr,gzarr

  dum_fname = "iondiff_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file = dum_fname, status="replace",action&
       &="write",iostat = ierr)

  IF(ierr /= 0) STOP "Ion diffusion file not found"

  WRITE(dumwrite,'(2(I0,1X),F14.8)') ioncnt, iontype, delta_t

! Ion Diffusion Analysis

  DO i = 0,nframes-1

     gxarr(i) = 0.0
     gyarr(i) = 0.0
     gzarr(i) = 0.0

  END DO

!$OMP PARALLEL DO PRIVATE(tinc,ifin,i,tim,j,rxcm,rycm,rzcm,aid)&
!$OMP&  REDUCTION(+:gxarr,gyarr,gzarr)
  DO tinc = 0, nframes-1

     ifin = nframes - tinc

     DO i = 1,ifin

        tim = i + tinc

        DO j = 1,ioncnt

           rxcm = itrx_gmx(j,tim) - itrx_gmx(j,i)
           rycm = itry_gmx(j,tim) - itry_gmx(j,i)
           rzcm = itrz_gmx(j,tim) - itrz_gmx(j,i)

           gxarr(tinc) = gxarr(tinc) + rxcm**2
           gyarr(tinc) = gyarr(tinc) + rycm**2
           gzarr(tinc) = gzarr(tinc) + rzcm**2


        END DO

     END DO

     gxarr(tinc) = gxarr(tinc)/(ifin*ioncnt)
     gyarr(tinc) = gyarr(tinc)/(ifin*ioncnt)
     gzarr(tinc) = gzarr(tinc)/(ifin*ioncnt)


  END DO
!$OMP END PARALLEL DO

  DO i = 0, nframes-1

     WRITE(dumwrite,"(4(F14.5,1X))") tarr_gmx(i+1), gxarr(i)&
          &,gyarr(i), gzarr(i)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE DIFF_IONS

!--------------------------------------------------------------------

SUBROUTINE DIFF_COUNTERIONS()

  USE PARAMS_GMX

  IMPLICIT NONE

  INTEGER :: i,j,tinc,ifin,tim,aid,atype,ierr,jout
  REAL    :: rxcm, rycm, rzcm
  REAL, DIMENSION(0:nframes-1) :: gxarr,gyarr,gzarr

  dum_fname = "countiondiff_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file = dum_fname, status="replace",action&
       &="write",iostat = ierr)

  IF(ierr /= 0) STOP "Counter-ion diffusion file not found"

  WRITE(dumwrite,'(2(I0,1X),F14.8)') c_ioncnt, c_iontype, delta_t

! Ion Diffusion Analysis

  DO i = 0,nframes-1

     gxarr(i) = 0.0
     gyarr(i) = 0.0
     gzarr(i) = 0.0

  END DO

! To do shifted time average for segmental diffusion

!$OMP PARALLEL DO PRIVATE(tinc,ifin,i,tim,j,rxcm,rycm,rzcm,aid)&
!$OMP&  REDUCTION(+:gxarr,gyarr,gzarr)
  DO tinc = 0, nframes-1

     ifin = nframes - tinc
     
     DO i = 1,ifin
        
        tim = i + tinc
      
        DO j = 1,c_ioncnt

           rxcm = ctrx_gmx(j,tim) - ctrx_gmx(j,i)
           rycm = ctry_gmx(j,tim) - ctry_gmx(j,i)
           rzcm = ctrz_gmx(j,tim) - ctrz_gmx(j,i)

           gxarr(tinc) = gxarr(tinc) + rxcm**2
           gyarr(tinc) = gyarr(tinc) + rycm**2
           gzarr(tinc) = gzarr(tinc) + rzcm**2
                      
        END DO

     END DO
     
     gxarr(tinc) = gxarr(tinc)/(REAL(ifin*c_ioncnt))
     gyarr(tinc) = gyarr(tinc)/(REAL(ifin*c_ioncnt))
     gzarr(tinc) = gzarr(tinc)/(REAL(ifin*c_ioncnt))

     
  END DO
!$OMP END PARALLEL DO

  DO i = 0, nframes-1

     WRITE(dumwrite,"(4(F14.5,1X))") tarr_gmx(i+1), gxarr(i)&
          &,gyarr(i),gzarr(i)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE DIFF_COUNTERIONS

!--------------------------------------------------------------------

SUBROUTINE DIFF_COMGROUP()

  USE PARAMS_GMX

  IMPLICIT NONE

  INTEGER :: i,j,tinc,ifin,tim,aid,atype,ierr,jout
  REAL    :: rxcm, rycm, rzcm
  REAL, DIMENSION(0:nframes-1) :: gxarr,gyarr,gzarr

  dum_fname = "COMdiff_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file = dum_fname, status="replace",action&
       &="write",iostat = ierr)

  IF(ierr /= 0) STOP "COM diffusion file not found"

  WRITE(dumwrite,'(2(I0,1X),F14.8)') ioncnt, iontype, delta_t

! Ion Diffusion Analysis

  DO i = 0,nframes-1

     gxarr(i) = 0.0
     gyarr(i) = 0.0
     gzarr(i) = 0.0

  END DO

!$OMP PARALLEL DO PRIVATE(tinc,ifin,i,tim,j,rxcm,rycm,rzcm,aid)&
!$OMP&  REDUCTION(+:gxarr,gyarr,gzarr)
  DO tinc = 0, nframes-1

     ifin = nframes - tinc

     DO i = 1,ifin

        tim = i + tinc

        DO j = 1,nmol_totcom

           rxcm = comtx_gmx(j,tim) - comtx_gmx(j,i)
           rycm = comty_gmx(j,tim) - comty_gmx(j,i)
           rzcm = comtz_gmx(j,tim) - comtz_gmx(j,i)

           gxarr(tinc) = gxarr(tinc) + rxcm**2
           gyarr(tinc) = gyarr(tinc) + rycm**2
           gzarr(tinc) = gzarr(tinc) + rzcm**2


        END DO

     END DO

     gxarr(tinc) = gxarr(tinc)/(ifin*nmol_totcom)
     gyarr(tinc) = gyarr(tinc)/(ifin*nmol_totcom)
     gzarr(tinc) = gzarr(tinc)/(ifin*nmol_totcom)


  END DO
!$OMP END PARALLEL DO

  DO i = 0, nframes-1

     WRITE(dumwrite,"(4(F14.5,1X))") tarr_gmx(i+1), gxarr(i)&
          &,gyarr(i), gzarr(i)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE DIFF_COMGROUP

!--------------------------------------------------------------------

!Ref: Borodin and Smith
!Macromolecules Vol: 39, No: 4, 1620-1629, 2006
! FROM THE PERSPECTIVE OF CATION
SUBROUTINE RESIDENCE_TIME_CATCOM()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,j,tval,a1id,a2id,ierr,ifin,tinc,tim
  REAL :: rxval,ryval,rzval,rval
  INTEGER,DIMENSION(1:nmol_totcom,nframes) :: autocf
  REAL,DIMENSION(0:nframes-1) :: tplot_cf

  dum_fname = "autocorrion_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file = dum_fname, status="replace",action&
       &="write",iostat = ierr)

  IF(ierr /= 0) STOP "ionpair residence time file not found"

!$OMP PARALLEL
!$OMP DO PRIVATE(i,j)
  DO j = 1,nframes

     DO i = 1,nmol_totcom

        autocf(i,j) = 0

     END DO
     
     tplot_cf(j-1) = 0

  END DO
!$OMP END DO

!$OMP DO PRIVATE(tval,i,j,a1id,a2id,rxval,ryval,rzval,rval)
  DO tval = 1,nframes 

     DO i = 1,ioncnt !populate autocorrelation fn array

        j = 1; a1id = ion_IDTYP_arr(i,1)

        IF(aidvals(a1id,3) .NE. iontype) THEN
           
           PRINT *, "Wrong atom type"
           PRINT *, tval, a1id, aidvals(a1id,3), iontype,&
                & ion_IDTYP_arr(i,1)
           STOP
           
        END IF

        ! See at least one anion COM is near each cation
        DO WHILE(j .LE. nmol_totcom)

           rxval = comtx_gmx(j,tval) - itrx_gmx(i,tval) 
           ryval = comty_gmx(j,tval) - itry_gmx(i,tval) 
           rzval = comtz_gmx(j,tval) - itrz_gmx(i,tval) 
           
           rxval = rxval - box_xl*ANINT(rxval/box_xl)
           ryval = ryval - box_yl*ANINT(ryval/box_yl)
           rzval = rzval - box_zl*ANINT(rzval/box_zl)
           
           rval = sqrt(rxval**2 + ryval**2 + rzval**2)
           
           IF(rval .LT. rcatCOM_cut) THEN
              
              autocf(i,tval) = 1
              j = nmol_totcom+1
              
           ELSE

              j = j +1

           END IF

        END DO

     END DO

  END DO
!$OMP END DO

!$OMP DO PRIVATE(tinc,ifin,tim,i,j)

  DO tinc = 0, nframes-1 !compute spectral product

     ifin = nframes - tinc
     
     DO i = 1,ifin
        
        tim = i + tinc
      
        DO j = 1,nmol_totcom

           tplot_cf(tinc) = tplot_cf(tinc) + REAL(autocf(j,tim)&
                &*autocf(j,i))
           
        END DO

     END DO

     tplot_cf(tinc) = tplot_cf(tinc)/REAL(ifin*ioncnt)
     
  END DO
!$OMP END DO
!$OMP END PARALLEL

  DO tinc = 0, nframes-1

     WRITE(dumwrite,"(I0,1X,F14.6)") tinc, tplot_cf(tinc)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE RESIDENCE_TIME_CATCOM

!--------------------------------------------------------------------

SUBROUTINE MAKE_QTARG_BOUNDS()

  USE PARAMS_GMX
  IMPLICIT NONE

  q_hi  = q_targ + q_tol; q_lo = q_targ - q_tol
  q_hi2 = q_hi*q_hi; q_lo2 = q_lo*q_lo

END SUBROUTINE MAKE_QTARG_BOUNDS

!--------------------------------------------------------------------

SUBROUTINE MAKE_QVECS()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER, ALLOCATABLE :: q_temp(:,:)
  INTEGER :: q_nx, q_ny, q_nz, q_cap
  INTEGER :: AllocateStatus, idx
  REAL :: fx, fy, fz, qx, qy, qz, qmag

  q_cap = (2*q_nmax+1)**3; q_possN = 0
  ALLOCATE(q_nlist(3,q_cap), stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate aidvals"

  q_nlist = 0.0; q_possN = 0

  ! Use the first frame to fix nx,ny,nz
  fx = pi2val/boxx_arr(1)
  fy = pi2val/boxy_arr(1)
  fz = pi2val/boxz_arr(1)

!$OMP PARALLEL DEFAULT(NONE) SHARED(q_nlist,fx,fy,fz,q_lo2,q_hi2,q_nmax)&
!$OMP& SHARED(q_possN) PRIVATE(q_nx,q_ny,q_nz,qx,qy,qz,qmag,idx)

!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  DO q_nx = -q_nmax, q_nmax
     DO q_ny = -q_nmax, q_nmax
        DO q_nz = -q_nmax, q_nmax

           IF(q_nx == 0 .AND. q_ny == 0 .AND. q_nz == 0) CYCLE

           qx = fx*q_nx; qy = fy*q_ny; qz = fz*q_nz
           qmag = qx*qx + qy*qy + qz*qz

           IF(qmag >= q_lo2 .AND. qmag <= q_hi2) THEN
              !$OMP ATOMIC CAPTURE
              idx = q_possN; q_possN = q_possN + 1
              !$OMP END ATOMIC
              idx = idx + 1
              q_nlist(1,idx) = q_nx
              q_nlist(2,idx) = q_ny
              q_nlist(3,idx) = q_nz

           END IF

        END DO
     END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  ! Allocate and deallocate q_nlist if big cap is present
  IF(q_possN /= q_cap) THEN

     ALLOCATE(q_temp(3,q_possN))
     q_temp(:,1:q_possN) = q_nlist(:,1:q_possN)
     DEALLOCATE(q_nlist)
     
     ALLOCATE(q_nlist(3,q_possN))
     q_nlist = q_temp
     DEALLOCATE(q_temp)
     
  END IF
  
END SUBROUTINE MAKE_QVECS
   
!--------------------------------------------------------------------

SUBROUTINE DYNAMIC_SELF_ION_FSKT()

  USE PARAMS_GMX
  IMPLICIT NONE

  REAL, DIMENSION(0:nframes-1) :: Fskt
  INTEGER :: qinc, ierr
  INTEGER :: i,j,tinc,tim,ifin
  REAL :: dx,dy,dz,qx,qy,qz,qdotr,sum_cos,accum
  REAL :: fx0, fy0, fz0
  CHARACTER(LEN=20) :: q_str

  WRITE(q_str,'(F5.2)') q_targ
  dum_fname = "Fsktion_"//trim(adjustl(q_str))//"_"&
       &//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file = dum_fname, status="replace",action&
       &="write",iostat = ierr)

  IF(ierr /= 0) STOP "Ion Fskt file not found"

  WRITE(dumwrite,'(2(I0,1X),5(F14.8,1X))') ioncnt, iontype, delta_t,&
       & q_hi, q_lo, q_tol, q_bin
  
  Fskt = 0.0

!$OMP PARALLEL DO PRIVATE(tinc,ifin,fx0,fy0,fz0,i,accum,tim,qinc,&
!$OMP& qx,qy,qz,j,dx,dy,dz,qdotr,sum_cos) REDUCTION(+:Fskt)
  DO tinc = 0, nframes-1
     
     ifin = nframes - tinc  
     fx0 = pi2val/boxx_arr(tinc+1)
     fy0 = pi2val/boxy_arr(tinc+1)
     fz0 = pi2val/boxz_arr(tinc+1)

     DO i = 1,ifin

        accum = 0.0
        tim = i + tinc
        
        DO qinc = 1, q_possN
           
           qx = fx0 * q_nlist(1,qinc)
           qy = fy0 * q_nlist(2,qinc)
           qz = fz0 * q_nlist(3,qinc)

           sum_cos = 0.0
           
           DO j = 1, ioncnt

              dx = itrx_gmx(j,tim) - itrx_gmx(j,i)
              dy = itry_gmx(j,tim) - itry_gmx(j,i)
              dz = itrz_gmx(j,tim) - itrz_gmx(j,i)

              qdotr = qx*dx + qy*dy + qz*dz
              sum_cos = sum_cos + cos(qdotr)
              
           END DO

           accum = accum + (sum_cos / real(ioncnt))
           
        END DO
        
     END DO

     Fskt(tinc) = Fskt(tinc) + (accum / real(q_possN*ifin))
     
  END DO
!$OMP END PARALLEL DO

  DO i = 0,nframes-1

     WRITE(dumwrite,"(2(F14.5,1X))") tarr_gmx(i+1), Fskt(i)

  END DO

  CLOSE(dumwrite)
  
END SUBROUTINE DYNAMIC_SELF_ION_FSKT

!--------------------------------------------------------------------

SUBROUTINE DYNAMIC_SELF_COUNTION_FSKT()

  USE PARAMS_GMX
  IMPLICIT NONE

  REAL, DIMENSION(0:nframes-1) :: Fskt
  INTEGER :: qinc, ierr
  INTEGER :: i,j,tinc,tim,ifin
  REAL :: dx,dy,dz,qx,qy,qz,qdotr,sum_cos,accum
  REAL :: fx0, fy0, fz0
  CHARACTER(LEN=20) :: q_str

  WRITE(q_str,'(F5.2)') q_targ
  dum_fname = "Fsktcountion_"//trim(adjustl(q_str))//"_"&
       &//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file = dum_fname, status="replace",action&
       &="write",iostat = ierr)

  IF(ierr /= 0) STOP "counter-ion Fskt file not found"

  WRITE(dumwrite,'(2(I0,1X),5(F14.8,1X))') c_ioncnt, c_iontype,&
       & delta_t, q_hi, q_lo, q_tol, q_bin
  
  Fskt = 0.0

!$OMP PARALLEL DO PRIVATE(tinc,ifin,fx0,fy0,fz0,i,accum,tim,qinc,&
!$OMP& qx,qy,qz,j,dx,dy,dz,qdotr,sum_cos) REDUCTION(+:Fskt)
  DO tinc = 0, nframes-1
     
     ifin = nframes - tinc  
     fx0 = pi2val/boxx_arr(tinc+1)
     fy0 = pi2val/boxy_arr(tinc+1)
     fz0 = pi2val/boxz_arr(tinc+1)

     DO i = 1,ifin

        accum = 0.0
        tim = i + tinc
        
        DO qinc = 1, q_possN
           
           qx = fx0 * q_nlist(1,qinc)
           qy = fy0 * q_nlist(2,qinc)
           qz = fz0 * q_nlist(3,qinc)

           sum_cos = 0.0
           
           DO j = 1, c_ioncnt

              dx = ctrx_gmx(j,tim) - ctrx_gmx(j,i)
              dy = ctry_gmx(j,tim) - ctry_gmx(j,i)
              dz = ctrz_gmx(j,tim) - ctrz_gmx(j,i)

              qdotr = qx*dx + qy*dy + qz*dz
              sum_cos = sum_cos + cos(qdotr)
              
           END DO

           accum = accum + (sum_cos / real(c_ioncnt))
           
        END DO
        
     END DO

     Fskt(tinc) = Fskt(tinc) + (accum / real(q_possN*ifin))
     
  END DO
!$OMP END PARALLEL DO

  DO i = 0,nframes-1

     WRITE(dumwrite,"(2(F14.5,1X))") tarr_gmx(i+1), Fskt(i)

  END DO

  CLOSE(dumwrite)

  
END SUBROUTINE DYNAMIC_SELF_COUNTION_FSKT

!--------------------------------------------------------------------

!!$SUBROUTINE DYNAMIC_CROSS_ION_WATER_FSKT()
!!$
!!$  USE PARAMS_GMX
!!$  IMPLICIT NONE
!!$
!!$  
!!$  
!!$
!!$END SUBROUTINE DYNAMIC_CROSS_FSKT

!--------------------------------------------------------------------


!Ref: Borodin and Smith
!Macromolecules Vol: 39, No: 4, 1620-1629, 2006
!Here we look at the residence time on the ANION
!FROM THE PERSPECTIVE OF ANION
SUBROUTINE RESIDENCE_TIME_ANCAT()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: i,j,tval,a1id,a2id,ierr,ifin,tinc,tim
  REAL :: rxval,ryval,rzval,rval
  INTEGER,DIMENSION(1:c_ioncnt,nframes) :: autocf
  REAL,DIMENSION(0:nframes-1) :: tplot_cf

  dum_fname = "autocorrcion_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file = dum_fname, status="replace",action&
       &="write",iostat = ierr)

  IF(ierr /= 0) STOP "c-ion residence time file not found"

!$OMP PARALLEL
!$OMP DO PRIVATE(i,j)
  DO j = 1,nframes

     DO i = 1,c_ioncnt

        autocf(i,j) = 0

     END DO
     
     tplot_cf(j-1) = 0

  END DO
!$OMP END DO

!$OMP DO PRIVATE(tval,i,j,a1id,a2id,rxval,ryval,rzval,rval)
  DO tval = 1,nframes 

     DO i = 1,c_ioncnt !populate autocorrelation fn array

        j = 1; a1id = countion_IDTYP_arr(i,1)
        
        IF(aidvals(a1id,3) .NE. c_iontype) THEN
           
           PRINT *, "Wrong counter-ion type"
           PRINT *, tval, a1id, aidvals(a1id,3), c_iontype,&
                & countion_IDTYP_arr(i,1)

        END IF

        DO WHILE(j .LE. ioncnt)

           a2id = ion_IDTYP_arr(j,1)

           IF(aidvals(a2id,3) .NE. iontype) THEN
           
              PRINT *, "Wrong ion type"
              PRINT *, tval, a2id, aidvals(a2id,3), iontype&
                   &,ion_IDTYP_arr(j,1)
              
           END IF
           
           rxval = ctrx_gmx(i,tval) - itrx_gmx(j,tval) 
           ryval = ctry_gmx(i,tval) - itry_gmx(j,tval) 
           rzval = ctrz_gmx(i,tval) - itrz_gmx(j,tval) 
           
           rxval = rxval - box_xl*ANINT(rxval/box_xl)
           ryval = ryval - box_yl*ANINT(ryval/box_yl)
           rzval = rzval - box_zl*ANINT(rzval/box_zl)
           
           rval = sqrt(rxval**2 + ryval**2 + rzval**2)
           
           IF(rval .LT. rcatan_cut) THEN
              
              autocf(i,tval) = 1
              j = ioncnt+1
              
           ELSE

              j = j +1

           END IF

        END DO

     END DO

  END DO
!$OMP END DO

!$OMP DO PRIVATE(tinc,ifin,tim,i,j)

  DO tinc = 0, nframes-1 !compute spectral product

     ifin = nframes - tinc
     
     DO i = 1,ifin
        
        tim = i + tinc
      
        DO j = 1,c_ioncnt

           tplot_cf(tinc) = tplot_cf(tinc) + REAL(autocf(j,tim)&
                &*autocf(j,i))
           
        END DO

     END DO

     tplot_cf(tinc) = tplot_cf(tinc)/REAL(ifin*c_ioncnt)
     
  END DO
!$OMP END DO
!$OMP END PARALLEL

  DO tinc = 0, nframes-1

     WRITE(dumwrite,"(I0,1X,F14.6)") tinc, tplot_cf(tinc)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE RESIDENCE_TIME_ANCAT

!--------------------------------------------------------------------

SUBROUTINE GET_VALUE_REAL(line, key, val, success)
  
  USE PARAMS_GMX
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN)  :: line, key
  REAL(8),      INTENT(OUT) :: val
  LOGICAL,      INTENT(OUT) :: success
  CHARACTER(len=:), ALLOCATABLE :: ex_s
  INTEGER :: ios, start_val, stop_val, buf

  val = 0.0

  CALL EXTRACT_VALSTRING(line,key,start_val,stop_val,buf,success)
  IF (.NOT. success) THEN
     PRINT *, "Could not extract time value .."
     PRINT *, line
     STOP
  END IF
  
  ALLOCATE(CHARACTER(LEN=buf) :: ex_s)
  ex_s = line(start_val:stop_val)
  

  READ(ex_s, *, iostat=ios) val
  success = (ios == 0)
  
END SUBROUTINE GET_VALUE_REAL
!--------------------------------------------------------------------

SUBROUTINE GET_VALUE_INT(line, key, val, success)

  USE PARAMS_GMX
  IMPLICIT NONE

  CHARACTER(*), INTENT(IN)  :: line, key
  INTEGER,      INTENT(OUT) :: val
  LOGICAL,      INTENT(OUT) :: success
  CHARACTER(LEN=:), ALLOCATABLE :: ex_s
  INTEGER :: ios,start_val, stop_val, buf

  val = 0
  CALL EXTRACT_VALSTRING(line,key,start_val,stop_val,buf,success)
  IF (.NOT. success) THEN
     PRINT *, "Could not extract step value .."
     PRINT *, line
     STOP
  END IF

  ALLOCATE(CHARACTER(LEN=buf) :: ex_s)
  ex_s = line(start_val:stop_val)
  
  READ(ex_s, *, iostat=ios) val
  success = (ios == 0)
  
END SUBROUTINE GET_VALUE_INT

!--------------------------------------------------------------------
! Generic code to extract string from ChatGPT
SUBROUTINE EXTRACT_VALSTRING(line,key,start_val,stop_val,buf,success)

  USE PARAMS_GMX
  IMPLICIT NONE

  CHARACTER(*), INTENT(IN)  :: line, key
  LOGICAL, INTENT(OUT) :: success
  INTEGER :: ex_p, ex_i, ex_j, ex_n
  INTEGER, INTENT(OUT) ::  start_val, stop_val, buf
  CHARACTER(len=*), PARAMETER :: delims = ' ,;'//achar(9)  ! space,
  ! comma, semicolon, tab
  success = .false.
  
  ex_n = len_trim(line)
  IF(ex_n == 0) RETURN
  
  ! Find "key="
  ex_p = index(line(1:ex_n), trim(key)//'=')
  IF (ex_p == 0) RETURN
 
  ! Start after "key="
  ex_i = ex_p + len_trim(key) + 1
  
  ! Skip leading spaces before the value
  DO WHILE (ex_i <= ex_n .AND. line(ex_i:ex_i) == ' ')
     ex_i = ex_i + 1
  END DO
  IF (ex_i > ex_n) RETURN
  
  ! Find next delimiter (or end)  
  ex_j = SCAN(line(ex_i:ex_n), delims)
  
  IF (ex_j == 0) THEN
     start_val = ex_i; stop_val = ex_n
  ELSE
     start_val = ex_i; stop_val = ex_i+ex_j-2
  END IF

  DO WHILE (start_val <= stop_val .AND. line(start_val:start_val) == '&
       & ')
     start_val = start_val + 1
  END DO
  IF(start_val > stop_val) RETURN

  buf = stop_val - start_val +1

  success = .true.
    
END SUBROUTINE EXTRACT_VALSTRING

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_TOPO_ARRAYS()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: AllocateStatus

  ! Allocate LAMMPS structure

  ALLOCATE(aidvals(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate aidvals"
  ALLOCATE(rxyz_gmx(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate rxyz_gmx"
  ALLOCATE(masses(ntotatomtypes,2),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate masses"

  ! Allocate box details

  ALLOCATE(boxx_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxx_arr"
  ALLOCATE(boxy_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxy_arr"
  ALLOCATE(boxz_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxz_arr"
  
  PRINT *, "Successfully allocated memory for topology"

END SUBROUTINE ALLOCATE_TOPO_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_COM_ARRAYS()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: AllocateStatus

  ! Allocate for COM

  IF(comflag) THEN
     ALLOCATE(comxyz_gmx(nmol_totcom,3),stat=AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate comxtz_gmx"
  ELSE
     ALLOCATE(comxyz_gmx(1,1),stat=AllocateStatus)
     DEALLOCATE(comxyz_gmx)
     ALLOCATE(com_IDTYP_arr(1,2),stat=AllocateStatus)
     DEALLOCATE(com_IDTYP_arr)
     ALLOCATE(comtyp_arr(1),stat = AllocateStatus)
     DEALLOCATE(comtyp_arr)
     ALLOCATE(com_masses(1),stat = AllocateStatus)
     DEALLOCATE(com_masses)
  END IF

END SUBROUTINE ALLOCATE_COM_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_ANALYSIS_ARRAYS()

  USE PARAMS_GMX
  IMPLICIT NONE

  INTEGER :: AllocateStatus

  ! Allocate for statics

  IF(rdfcalc_flag) THEN
     ALLOCATE(rdfarray(0:rmaxbin-1,npairs),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate rdfarray"
  ELSE
     ALLOCATE(rdfarray(1,1),stat = AllocateStatus)
     DEALLOCATE(rdfarray)
  END IF

  
  IF(catan_neighcalc_flag) THEN
     ALLOCATE(cat_an_neighavg(1:maxneighsize),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate cat_an_neighavg"
     ALLOCATE(an_cat_neighavg(1:maxneighsize),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate an_cat_neighavg"
  ELSE
     ALLOCATE(cat_an_neighavg(1),stat = AllocateStatus)
     ALLOCATE(an_cat_neighavg(1),stat = AllocateStatus)
     DEALLOCATE(cat_an_neighavg)
     DEALLOCATE(an_cat_neighavg)
  END IF


  ! Allocate for dynamics 

  ALLOCATE(tarr_gmx(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate tarr_gmx"

  IF(ion_dynflag) THEN
     ALLOCATE(itrx_gmx(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate itrx_gmx"
     ALLOCATE(itry_gmx(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate itry_gmx"
     ALLOCATE(itrz_gmx(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate itrz_gmx"
  ELSE
     ALLOCATE(itrx_gmx(1,1))
     ALLOCATE(itry_gmx(1,1))
     ALLOCATE(itrz_gmx(1,1))
     DEALLOCATE(itrx_gmx)
     DEALLOCATE(itry_gmx)
     DEALLOCATE(itrz_gmx)
  END IF

  IF(cion_dynflag) THEN
     ALLOCATE(ctrx_gmx(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ctrx_gmx"
     ALLOCATE(ctry_gmx(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ctry_gmx"
     ALLOCATE(ctrz_gmx(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ctrz_gmx"
  ELSE
     ALLOCATE(ctrx_gmx(1,1))
     ALLOCATE(ctry_gmx(1,1))
     ALLOCATE(ctrz_gmx(1,1))
     DEALLOCATE(ctrx_gmx)
     DEALLOCATE(ctry_gmx)
     DEALLOCATE(ctrz_gmx)
  END IF
  
  IF(com_dynflag) THEN
     ALLOCATE(comtx_gmx(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate comtx_gmx"
     ALLOCATE(comty_gmx(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate comty_gmx"
     ALLOCATE(comtz_gmx(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate comtz_gmx"
  ELSE
     ALLOCATE(comtx_gmx(1,1))
     ALLOCATE(comty_gmx(1,1))
     ALLOCATE(comtz_gmx(1,1))
     DEALLOCATE(comtx_gmx)
     DEALLOCATE(comty_gmx)
     DEALLOCATE(comtz_gmx)
  END IF

  IF(dynfsktflag == 0) THEN
     ALLOCATE(q_nlist(1,1),stat=AllocateStatus)
     DEALLOCATE(q_nlist)
  END IF
  
  PRINT *, "Successfully allocated memory for analyis"

END SUBROUTINE ALLOCATE_ANALYSIS_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE DEALLOCATE_ARRAYS()

  USE PARAMS_GMX

  IMPLICIT NONE

  !Global arrays
  DEALLOCATE(aidvals)
  DEALLOCATE(rxyz_gmx)
  DEALLOCATE(boxx_arr)
  DEALLOCATE(boxy_arr)
  DEALLOCATE(boxz_arr)
  DEALLOCATE(masses)
  
  !Statics calculations arrays
  IF(rdfcalc_flag) DEALLOCATE(rdfarray)

  !Dynamic calculations arrays
  IF(ion_dynflag) THEN
     DEALLOCATE(itrx_gmx)
     DEALLOCATE(itry_gmx)
     DEALLOCATE(itrz_gmx)
  END IF

  IF(cion_dynflag) THEN
     DEALLOCATE(ctrx_gmx)
     DEALLOCATE(ctry_gmx)
     DEALLOCATE(ctrz_gmx)
  END IF

  ! COM arrays
  IF(comflag) THEN
     DEALLOCATE(comxyz_gmx)
  END IF
  
  IF(com_dynflag) THEN
     DEALLOCATE(comtx_gmx)
     DEALLOCATE(comty_gmx)
     DEALLOCATE(comtz_gmx)
  END IF

  
END SUBROUTINE DEALLOCATE_ARRAYS

!----------------------------------------------------------------
