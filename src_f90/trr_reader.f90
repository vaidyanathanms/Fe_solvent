SUBROUTINE TRR_READER(fname,write_xyz,xyz_out)

  USE iso_fortran_env, only: int8, int32, int64, real32, real64
  IMPLICIT NONE

  character(len=:), allocatable :: fname, xyz_out
  logical, save :: little = is_little_endian()
  logical :: write_xyz = .false.
  integer :: i, ios
  integer :: unit
  integer(int64) :: fsize, pos
  integer(int32) :: ir_size, e_size, box_size, vir_size, pres_size, top_size, sym_size
  integer(int32) :: x_size, v_size, f_size
  integer(int32) :: natoms, step, nre
  real(real64)   :: time_d, lambda_d
  real(real32)   :: time_f, lambda_f
  integer :: rbytes
  real(real64), allocatable :: boxD(:,:), xD(:,:), vD(:,:), fD(:,:)
  real(real32), allocatable :: boxS(:,:), xS(:,:), vS(:,:), fS(:,:)
  logical :: have_box, have_x, have_v, have_f
  integer :: nframes
  integer :: argn
  
  call parse_args(fname, write_xyz, xyz_out)
  call open_stream(fname, unit, fsize)

  nframes = 0
  pos = 0_int64

  do
     if (.not. can_read_more(unit)) exit

     ! --- Read TRR header fields (XDR big-endian int32) ---
     if (.not. read_i32(unit, ir_size)) exit
     if (.not. read_i32(unit, e_size))  exit
     if (.not. read_i32(unit, box_size)) exit
     if (.not. read_i32(unit, vir_size)) exit
     if (.not. read_i32(unit, pres_size)) exit
     if (.not. read_i32(unit, top_size)) exit
     if (.not. read_i32(unit, sym_size)) exit
     if (.not. read_i32(unit, x_size)) exit
     if (.not. read_i32(unit, v_size)) exit
     if (.not. read_i32(unit, f_size)) exit
     
     if (.not. read_i32(unit, natoms)) exit
     if (.not. read_i32(unit, step))   exit
     if (.not. read_i32(unit, nre))    exit
     
     ! Determine real size (4 or 8) from available block sizes
     rbytes = determine_real_bytes(natoms, x_size, v_size, f_size, box_size)
     
     ! Read time, lambda in the detected precision
     select case (rbytes)
     case (4)
        if (.not. read_r32(unit, time_f))    exit
        if (.not. read_r32(unit, lambda_f))  exit
        time_d = real(time_f, real64)
        lambda_d = real(lambda_f, real64)
     case (8)
        if (.not. read_r64(unit, time_d))    exit
        if (.not. read_r64(unit, lambda_d))  exit
     case default
        write(*,*) 'ERROR: Unable to infer real size (not 4 or 8).'
        exit
     end select
     
     have_box = (box_size > 0)
     have_x   = (x_size   > 0)
     have_v   = (v_size   > 0)
     have_f   = (f_size   > 0)
     
     ! --- Read box (3x3) if present ---
     if (have_box) then
        if (rbytes == 4) then
           allocate(boxS(3,3))
           if (.not. read_r32_array(unit, boxS, 9)) call fail('Failed reading box (S)')
        else
           allocate(boxD(3,3))
           if (.not. read_r64_array(unit, boxD, 9)) call fail('Failed reading box (D)')
        end if
     end if

     ! --- Skip other matrix blocks if present (virial, pressure, etc.) ---
     call skip_bytes(unit, vir_size)
     call skip_bytes(unit, pres_size)
     call skip_bytes(unit, ir_size)
     call skip_bytes(unit, e_size)
     call skip_bytes(unit, top_size)
     call skip_bytes(unit, sym_size)
     
     ! --- Read coordinates / velocities / forces ---
     if (have_x) then
        if (rbytes == 4) then
           allocate(xS(3, natoms))
           if (.not. read_r32_array(unit, xS, 3*natoms)) call fail('Failed reading X (S)')
        else
           allocate(xD(3, natoms))
           if (.not. read_r64_array(unit, xD, 3*natoms)) call fail('Failed reading X (D)')
        end if
     end if
     if (have_v) then
        if (rbytes == 4) then
           allocate(vS(3, natoms))
           if (.not. read_r32_array(unit, vS, 3*natoms)) call fail('Failed reading V (S)')
        else
           allocate(vD(3, natoms))
           if (.not. read_r64_array(unit, vD, 3*natoms)) call fail('Failed reading V (D)')
        end if
     end if
     if (have_f) then
        if (rbytes == 4) then
           allocate(fS(3, natoms))
           if (.not. read_r32_array(unit, fS, 3*natoms)) call fail('Failed reading F (S)')
        else
           allocate(fD(3, natoms))
           if (.not. read_r64_array(unit, fD, 3*natoms)) call fail('Failed reading F (D)')
        end if
     end if

     nframes = nframes + 1
     write(*,'(a,i6,2x,a,i10,2x,a,f12.5,2x,a,i8)') 'Frame:', nframes,&
          & 'step:', step, 'time:', time_d, 'natoms:', natoms

     ! --- Optional: write XYZ per frame (coordinates only) ---
     if (write_xyz .and. have_x) then
        call append_xyz(xyz_out, nframes, natoms, time_d, rbytes, xS,&
             & xD)
     end if

     ! Cleanup arrays for next loop
     if (allocated(boxS)) deallocate(boxS)
     if (allocated(boxD)) deallocate(boxD)
     if (allocated(xS))   deallocate(xS)
     if (allocated(xD))   deallocate(xD)
     if (allocated(vS))   deallocate(vS)
     if (allocated(vD))   deallocate(vD)
     if (allocated(fS))   deallocate(fS)
     if (allocated(fD))   deallocate(fD)
  end do

  close(unit)
  write(*,*) 'Done. Total frames read:', nframes

END SUBROUTINE TRR_READER

!---------------------------------------------------------------------

SUBROUTINE PARSE_ARGS(fname, write_xyz, xyz_out)

  IMPLICIT NONE
  
  character(len=:), allocatable, intent(out) :: fname, xyz_out
  logical, intent(out) :: write_xyz
  integer :: n
  character(len=1024) :: arg
  
  n = command_argument_count()
  if (n < 1) then
     write(*,*) 'Usage: trr_reader <traj.trr> [--xyz out.xyz]'
     stop 1
  end if
  call get_command_argument(1, arg)
  fname = trim(arg)
  write_xyz = .false.
  xyz_out = ''
  if (n >= 3) then
     call get_command_argument(2, arg)
     if (trim(arg) == '--xyz') then
        call get_command_argument(3, arg)
        xyz_out = trim(arg)
        write_xyz = .true.
     end if
  end if
END SUBROUTINE PARSE_ARGS

!---------------------------------------------------------------------

SUBROUTINE OPEN_STREAM(fname, unit, fsize)

  IMPLICIT NONE
  
  character(len=*), intent(in) :: fname
  integer, intent(out) :: unit
  integer(int64), intent(out) :: fsize
  integer :: ios
  inquire(file=fname, exist=ios)
  if (ios == 0) then
     write(*,*) 'ERROR: file not found: ', trim(fname)
     stop 1
  end if
  open(newunit=unit, file=fname, access='stream', form&
       &='unformatted', status='old', action='read', iostat=ios)
  if (ios /= 0) then
     write(*,*) 'ERROR: cannot open file: ', trim(fname)
     stop 1
  end if
  inquire(unit=unit, size=fsize)
END SUBROUTINE OPEN_STREAM
!---------------------------------------------------------------------
LOGICAL FUNCTION CAN_READ_MORE(unit)
  integer, intent(in) :: unit
  integer(int64) :: pos, size
  inquire(unit=unit, pos=pos, size=size)
  can_read_more = (pos < size)
END FUNCTION CAN_READ_MORE
!---------------------------------------------------------------------
  
SUBROUTINE FAIL(msg)
  IMPLICIT NONE
  character(len=*), intent(in) :: msg
  write(*,*) 'ERROR: ', msg
  stop 1
END SUBROUTINE FAIL

!---------------------------------------------------------------------
! ---------- XDR helpers (big-endian) ----------
PURE FUNCTION IS_LITTLE_ENDIAN() result(lil)

  IMPLICIT NONE
  logical :: lil
  integer(int32) :: one
  character(len=4) :: b
  one = 1_int32
  b = transfer(one, b)
  ! On little-endian, lowest-address byte is 1
  lil = (iachar(b(1:1)) == 1)
END FUNCTION IS_LITTLE_ENDIAN
!---------------------------------------------------------------------  
SUBROUTINE BSWAP4(buf)
  IMPLICIT NONE
  integer(int8), intent(inout) :: buf(4)
  integer(int8) :: t
  t=buf(1); buf(1)=buf(4); buf(4)=t
  t=buf(2); buf(2)=buf(3); buf(3)=t
end subroutine bswap4
!---------------------------------------------------------------------  
SUBROUTINE BSWAP8(buf)
  IMPLICIT NONE
  integer(int8), intent(inout) :: buf(8)
  integer(int8) :: t
  t=buf(1); buf(1)=buf(8); buf(8)=t
  t=buf(2); buf(2)=buf(7); buf(7)=t
  t=buf(3); buf(3)=buf(6); buf(6)=t
  t=buf(4); buf(4)=buf(5); buf(5)=t
END SUBROUTINE BSWAP8
!---------------------------------------------------------------------    
LOGICAL FUNCTION READ_i32(unit,val)
  IMPLICIT NONE
  integer, intent(in) :: unit
  integer(int32), intent(out) :: val
  integer(int8) :: b(4)
  integer :: ios
  read_i32 = .false.
  read(unit, iostat=ios) b
  if (ios /= 0) return
  if (little) call bswap4(b)
  val = transfer(b, val)
  read_i32 = .true.
END FUNCTION READ_i32
!---------------------------------------------------------------------
LOGICAL FUNCTION READ_r32(unit,val)

  IMPLICIT NONE
  integer, intent(in) :: unit
  real(real32), intent(out) :: val
  integer(int8) :: b(4)
  integer :: ios
  read_r32 = .false.
  read(unit, iostat=ios) b
  if (ios /= 0) return
  if (little) call bswap4(b)
  val = transfer(b, val)
  read_r32 = .true.
END FUNCTION READ_r32
!---------------------------------------------------------------------
LOGICAL FUNCTION READ_r64(unit,val)

  IMPLICIT NONE
  integer, intent(in) :: unit
  real(real64), intent(out) :: val
  integer(int8) :: b(8)
  integer :: ios
  read_r64 = .false.
  read(unit, iostat=ios) b
  if (ios /= 0) return
  if (little) call bswap8(b)
  val = transfer(b, val)
  read_r64 = .true.
END FUNCTION READ_r64
!---------------------------------------------------------------------
LOGICAL FUNCTiON READ_r32_ARRAY(unit, arr, n) result(ok)
  IMPLICIT NONE
  integer, intent(in) :: unit
  real(real32), intent(out) :: arr(:,:)
  integer, intent(in) :: n
  integer :: ios
  integer(int8), allocatable :: b(:)
  integer :: i
  ok = .false.
  allocate(b(4*n))
  read(unit, iostat=ios) b
  if (ios /= 0) then
     deallocate(b); return
  end if
  if (little) then
     do i=1, n
        call bswap4(b(4*(i-1)+1:4*i))
     end do
  end if
  call bytes_to_r32(arr, b, n)
  deallocate(b)
  ok = .true.
END FUNCTiON READ_r32_ARRAY
!---------------------------------------------------------------------
LOGICAL FUNCTION READ_r64_ARRAY(unit, arr, n) result(ok)
  IMPLICIT NONE
  integer, intent(in) :: unit
  real(real64), intent(out) :: arr(:,:)
  integer, intent(in) :: n
  integer :: ios
  integer(int8), allocatable :: b(:)
  integer :: i
  ok = .false.
  allocate(b(8*n))
  read(unit, iostat=ios) b
  if (ios /= 0) then
     deallocate(b); return
  end if
  if (little) then
     do i=1, n
        call bswap8(b(8*(i-1)+1:8*i))
     end do
  end if
  call bytes_to_r64(arr, b, n)
  deallocate(b)
  ok = .true.
END FUNCTION READ_r64_ARRAY
!---------------------------------------------------------------------
SUBROUTINE BYTES_to_r32(arr, b, n)
  IMPLICIT NONE
  real(real32), intent(out) :: arr(:,:)
  integer(int8), intent(in) :: b(:)
  integer, intent(in) :: n
  real(real32), allocatable :: tmp(:)
  integer :: i, j, k, na, nb
  allocate(tmp(n))
  tmp = transfer(b, tmp)
  ! arr is (3, natoms) or (3,3) depending on caller; fill column-major
  na = size(arr,1); nb = size(arr,2)
  if (na*nb /= n) call fail('bytes_to_r32 size mismatch')
  k = 0
  do j=1, nb
     do i=1, na
        k = k + 1
        arr(i,j) = tmp(k)
     end do
  end do
  deallocate(tmp)
END SUBROUTINE BYTES_to_r32
!---------------------------------------------------------------------
SUBROUTINE BYTES_TO_r64(arr, b, n)

  IMPLICIT NONE
  real(real64), intent(out) :: arr(:,:)
  integer(int8), intent(in) :: b(:)
  integer, intent(in) :: n
  real(real64), allocatable :: tmp(:)
  integer :: i, j, k, na, nb
  allocate(tmp(n))
  tmp = transfer(b, tmp)
  na = size(arr,1); nb = size(arr,2)
  if (na*nb /= n) call fail('bytes_to_r64 size mismatch')
  k = 0
  do j=1, nb
     do i=1, na
        k = k + 1
        arr(i,j) = tmp(k)
     end do
  end do
  deallocate(tmp)
  
END SUBROUTINE BYTES_TO_r64
!---------------------------------------------------------------------
SUBROUTINE SKIP_BYTES(unit, nbytes)
  IMPLICIT NONE
  integer, intent(in) :: unit
  integer(int32), intent(in) :: nbytes
  integer(int8), allocatable :: dummy(:)
  integer :: ios
  if (nbytes <= 0) return
  allocate(dummy(nbytes))
  read(unit, iostat=ios) dummy
  if (ios /= 0) call fail('Failed skipping bytes')
  deallocate(dummy)
END SUBROUTINE SKIP_BYTES
!---------------------------------------------------------------------
INTEGER FUNCTION DETERMINE_REAL_BYTES(natoms, x_size, v_size, f_size, box_size) result(n)

  IMPLICIT NONE
  integer(int32), intent(in) :: natoms, x_size, v_size, f_size, box_size
  if (x_size > 0) then
     if (mod(x_size, 3*natoms) == 0) then
        n = x_size / (3*natoms); return
     end if
  end if
  if (v_size > 0) then
     if (mod(v_size, 3*natoms) == 0) then
        n = v_size / (3*natoms); return
     end if
  end if
  if (f_size > 0) then
     if (mod(f_size, 3*natoms) == 0) then
        n = f_size / (3*natoms); return
     end if
  end if
  if (box_size > 0) then
     if (mod(box_size, 9) == 0) then
        n = box_size / 9; return
     end if
  end if
  n = 4  ! fallback
  
END FUNCTION DETERMINE_REAL_BYTES
!---------------------------------------------------------------------
SUBROUTINE APPEND_XYZ(xyz_out, frame_id, natoms, time_d, rbytes, xS, xD)

  IMPLICIT NONE
  character(len=*), intent(in) :: xyz_out
  integer, intent(in) :: frame_id, natoms, rbytes
  real(real64), intent(in) :: time_d
  real(real32), intent(in), optional :: xS(:,:)
  real(real64), intent(in), optional :: xD(:,:)
  integer :: u, i
  integer :: ios
  open(newunit=u, file=xyz_out, position='append', action='write',&
       & status='unknown', iostat=ios)
  if (ios /= 0) then
     write(*,*) 'WARNING: could not open XYZ for append: ',&
          & trim(xyz_out)
     return
  end if
  write(u,'(i0)') natoms
  write(u,'(a,i0,a,f12.5)') 'frame=', frame_id, '  time=', time_d
  if (rbytes == 4 .and. present(xS)) then
     do i=1, natoms
        write(u,'(a,1x,3f16.8)') 'X', real(xS(1,i),real64),&
             & real(xS(2,i),real64), real(xS(3,i),real64)
     end do
  else if (rbytes == 8 .and. present(xD)) then
     do i=1, natoms
        write(u,'(a,1x,3f16.8)') 'X', xD(1,i), xD(2,i), xD(3,i)
     end do
  else
     write(*,*) 'WARNING: no coordinates to write for frame ',&
          & frame_id
  end if
  close(u)
END SUBROUTINE APPEND_XYZ
!---------------------------------------------------------------------


