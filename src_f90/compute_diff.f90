PROGRAM MAIN

  USE ANALYZE_ALL
  IMPLICIT NONE

  !$OMP PARALLEL
  nproc = OMP_GET_NUM_THREADS()
  PRINT *, "Number of threads are: ", nproc
  !$OMP END PARALLEL

  skipfr = 0
  ntotatoms = 90
  ntotbonds = 0
  ntotangls = 0
  ntotdihds = 0
  ntotatomtypes = 1
  ntotbondtypes = 0
  ntotdihdtypes = 0
  ntotimprtypes = 0

  ! Frame, molecules and processor details
  nframes = 3300; skipfr = 0; freqfr = 0; nfrcntr = 0
  nchains = 0; atperchain = 0

  
  CALL DEFAULTVALUES()
  CALL ALLOCATE_TOPO_ARRAYS()
  CALL ALLOCATE_ANALYSIS_ARRAYS()
  CALL ANALYZE_TRAJECTORYFILE()

  
END PROGRAM MAIN

!--------------------------------------------------------------------

SUBROUTINE DEFAULTVALUES()

  USE ANALYZE_ALL
  IMPLICIT NONE


  ! Initialize flags
  polyflag = 0
  rgall = 0; rgcalc = 0; rdfcalc = 0
  ion_dynflag = 0; cion_dynflag = 0; pion_dynflag = 0
  ion_diff = 1; cion_diff = 0; pion_diff = 0
  bfrdf_calc = 0
  catan_autocfflag = 0; catpol_autocfflag = 0

  ! Initialize iontypes
  c_iontype = -1; p_iontype = -1; iontype = -1

  !Initialize system quantities
  npoly_types = 0; ioncnt = 0; c_ioncnt = 0; p_ioncnt= 0

  ! Initialize distributions and frequencies
  rdffreq = 0; rgfreq = 1

  ! Initialzie structural quantities
  rdomcut = 10.0;  rmaxbin = 100; rbinval = REAL(rdomcut)&
       &/REAL(rmaxbin)
  rcatan_cut = 0.0; rneigh_cut = 0.0

  ! Initialize structural averages
  rvolavg = 0; rgavg = 0

  ! Initialize dynamical quantities
  rcatpol_cut1 = 0.0; rcatpol_cut2 = 0.0
END SUBROUTINE DEFAULTVALUES

!--------------------------------------------------------------------


SUBROUTINE ANALYZE_TRAJECTORYFILE()

  USE ANALYZE_ALL

  IMPLICIT NONE

  INTEGER :: i,j,aid,ierr,atchk,atype,jumpfr,jout
  REAL :: xlo,xhi,ylo,yhi,zlo,zhi

  traj_fname = "com_0.1.lammpstrj"
  OPEN(unit = 15,file =traj_fname,action="read",status="old"&
       &,iostat=ierr)

  IF(ierr /= 0) STOP "trajectory file not found"

  PRINT *, "Trajectory file used: ",trim(adjustl(traj_fname))
  WRITE(logout,*) "Trajectory file used: "&
       &,trim(adjustl(traj_fname))

  
  PRINT *, "Analyzing trajectory file..."

  print *, skipfr, nframes
  DO i = 1,skipfr

     DO j = 1,ntotatoms+9

        READ(15,*) 

     END DO

     IF(mod(i,100) == 0) PRINT *, "Skipped ", i, "frames"

  END DO

  DO i = 1,nframes

     nfrcntr = nfrcntr + 1
     IF(mod(i,100) == 0) PRINT *, "Processing ", i+1,"th frame"

     READ(15,*)
     READ(15,*) timestep

     tarr_lmp(i) = timestep
     
     READ(15,*) 
     READ(15,*) atchk

     READ(15,*) 
     READ(15,*) xlo, xhi
     READ(15,*) ylo, yhi
     READ(15,*) zlo, zhi

     READ(15,*)

     box_xl = xhi - xlo
     box_yl = yhi - ylo
     box_zl = zhi - zlo
     
     boxx_arr(i)  = box_xl
     boxy_arr(i)  = box_yl
     boxz_arr(i)  = box_zl

     DO j = 1,atchk

        READ(15,*) aid,atype,rxyz_lmp(aid,1),rxyz_lmp(aid,2)&
             &,rxyz_lmp(aid,3)

        itrx_lmp(aid,i) = rxyz_lmp(aid,1)
        itry_lmp(aid,i) = rxyz_lmp(aid,2)
        itrz_lmp(aid,i) = rxyz_lmp(aid,3)

     END DO
         
     DO jumpfr = 1,freqfr

        READ(15,*)
        READ(15,*)        
        READ(15,*)
 
        READ(15,*) atchk

        DO j = 1,atchk+5

           READ(15,*) 

        END DO

     END DO

  END DO

  CLOSE(15)

  PRINT *, "Beginning dynamical analysis..."
  CALL DIFF_IONS()

END SUBROUTINE ANALYZE_TRAJECTORYFILE

!--------------------------------------------------------------------


SUBROUTINE DIFF_IONS()

  USE ANALYZE_ALL

  IMPLICIT NONE

  INTEGER :: i,j,tinc,ifin,tim,aid,atype,ierr,jout
  REAL    :: rxcm, rycm, rzcm
  REAL, DIMENSION(0:nframes-1) :: gxarr,gyarr,gzarr

  dum_fname = "iondiff_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file = dum_fname, status="replace",action&
       &="write",iostat = ierr)

  IF(ierr /= 0) STOP "Ion diffusion file not found"

  WRITE(dumwrite,'(2(I0,1X),F14.8)') ntotatoms, 1, delta_t
  print *, "starting diff calc"
  
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

        DO j = 1,ntotatoms

           rxcm = itrx_lmp(j,tim) - itrx_lmp(j,i)
           rycm = itry_lmp(j,tim) - itry_lmp(j,i)
           rzcm = itrz_lmp(j,tim) - itrz_lmp(j,i)

           gxarr(tinc) = gxarr(tinc) + rxcm**2
           gyarr(tinc) = gyarr(tinc) + rycm**2
           gzarr(tinc) = gzarr(tinc) + rzcm**2


        END DO

     END DO

     gxarr(tinc) = gxarr(tinc)/(ifin*ntotatoms)
     gyarr(tinc) = gyarr(tinc)/(ifin*ntotatoms)
     gzarr(tinc) = gzarr(tinc)/(ifin*ntotatoms)


  END DO
!$OMP END PARALLEL DO

  DO i = 0, nframes-1

     WRITE(dumwrite,"(I10,1X,3(F14.5,1X))") tarr_lmp(i+1), gxarr(i)&
          &,gyarr(i), gzarr(i)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE DIFF_IONS

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_TOPO_ARRAYS()

  USE ANALYZE_ALL
  IMPLICIT NONE

  INTEGER :: AllocateStatus

! Allocate LAMMPS structure

  ALLOCATE(aidvals(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate aidvals"
  ALLOCATE(rxyz_lmp(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate rxyz_lmp"
  ALLOCATE(charge_lmp(ntotatoms,1),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate charge_lmp"
  ALLOCATE(vel_xyz(ntotatoms,4),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate vel_xyz"

  IF(ntotbonds /= 0) THEN
     ALLOCATE(bond_lmp(ntotbonds,4),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate bond_lmp"
  ELSE
     PRINT *, "Warning: No bonds - Not correct for bonded systems"
     ALLOCATE(bond_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(bond_lmp)
  END IF
  
  IF(ntotangls /= 0) THEN
     ALLOCATE(angl_lmp(ntotangls,5),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate angl_lmp"
  ELSE
     ALLOCATE(angl_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(angl_lmp)
  END IF
     
  IF(ntotdihds /= 0) THEN
     ALLOCATE(dihd_lmp(ntotdihds,6),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate dihd_lmp"
  ELSE
     ALLOCATE(dihd_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(dihd_lmp)
  END IF
  
  IF(ntotimprs /= 0) THEN
     ALLOCATE(impr_lmp(ntotimprs,6),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate zlmp"
  ELSE
     ALLOCATE(impr_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(impr_lmp)
  END IF

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

SUBROUTINE ALLOCATE_ANALYSIS_ARRAYS()

  USE ANALYZE_ALL
  IMPLICIT NONE

  INTEGER :: AllocateStatus

! Allocate for dynamics 

  ALLOCATE(tarr_lmp(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate tarr_lmp"
  

  ALLOCATE(itrx_lmp(ntotatoms,nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate itrx_lmp"
  ALLOCATE(itry_lmp(ntotatoms,nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate itry_lmp"
  ALLOCATE(itrz_lmp(ntotatoms,nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate itrz_lmp"

  PRINT *, "Successfully allocated memory for analyis"

END SUBROUTINE ALLOCATE_ANALYSIS_ARRAYS

!--------------------------------------------------------------------
