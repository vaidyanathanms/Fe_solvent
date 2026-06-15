! Program to read GROMACS trr files
! See https://github.com/wesbarnett/libgmxfort
! https://github.com/chen3262/xdrfile_fortran/blob/master/trr-interface.f90

PROGRAM READ_TRR_PROG
  
  USE, INTRINSIC :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer
  USE TRR
  
  IMPLICIT NONE
  CHARACTER (len=1024) :: fntemp, filename
  REAL, ALLOCATABLE :: pos(:,:), vel(:,:), force(:,:)
  INTEGER :: NATOMS, STEP, STAT, i
  REAL :: box(3,3), lambda, time, box_trans(3,3)
  TYPE(C_PTR) :: xd_c
  TYPE(xdrfile), pointer :: xd
  LOGICAL :: ex

  ! Set the file name for C.
  CALL GETARG(1,fntemp)
  filename = trim(fntemp)//C_NULL_CHAR
  
  inquire(file=trim(filename),exist=ex)
  
  if (ex .eqv. .false.) then
     write(0,*)
     write(0,'(a)') " Error: "//trim(filename)//" does not exist."
     write(0,*)
     stop
  end if

  STAT = READ_TRR_NATOMS(filename,NATOMS)
  allocate(pos(3,NATOMS))
  allocate(vel(3,NATOMS))
  allocate(force(3,NATOMS))
  
  ! Open the file for reading. Convert C pointer to Fortran pointer.
  xd_c = xdrfile_open(filename,"r")
  call c_f_pointer(xd_c,xd)

  STAT = read_trr(xd,NATOMS,STEP,time,lambda,box_trans,pos,vel,force)
  
  do while ( STAT == 0 )
     
     ! C is row-major, whereas Fortran is column major. Hence the following.
     box = transpose(box_trans)
     
     ! Just an example to show what was read in (Modify this part with your own fortran code)
     write(*,'(a,f12.6,a,i0)') " Time (ps): ", time, "  Step: ", STEP
     write(*,'(a,f12.6,a,i0)') " Lambda: ", lambda, "  No. Atoms: ", NATOMS
     
     do i=1,NATOMS
        write(*,'(3f8.3,3f8.4,3f8.4)') pos(:,i), vel(:,i), force(:,i)
     end do
     ! This is the same order as found in the GRO format fyi
     write(*,'(3f9.5)') box(1,1), box(2,2), box(3,3)
     
     STAT = read_trr(xd,NATOMS,STEP,time,lambda,box_trans,pos,vel,force)
     
  end do

  STAT = xdrfile_close(xd)
  deallocate(pos)
  deallocate(vel)
  deallocate(force)
  
END PROGRAM READ_TRR_PROG
