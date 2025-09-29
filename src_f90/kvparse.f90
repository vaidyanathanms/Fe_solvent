!Using ChatGPT to parse key and value
!Date: Sept-23-2025
!Use to parse the header line to get the time value
MODULE KEY_VALUE_PARSE

  USE OMP_LIB
  IMPLICIT NONE
!!$  interface get_value
!!$     module procedure get_value_real
!!$     module procedure get_value_int
!!$     module procedure get_value_char
!!$  end interface get_value

contains
  ! Return the substring value for key=... ; success=.true. if found & parsed
  pure subroutine get_value_real(line, key, val, success)
    character(*), intent(in)  :: line, key
    real(8),      intent(out) :: val
    logical,      intent(out) :: success
    character(len=:), allocatable :: s
    integer :: ios
    call extract_valstring(line, key, s, success)
    if (.not. success) return
    read(s, *, iostat=ios) val
    success = (ios == 0)
  end subroutine get_value_real

  pure subroutine get_value_int(line, key, val, success)
    character(*), intent(in)  :: line, key
    integer,      intent(out) :: val
    logical,      intent(out) :: success
    character(len=:), allocatable :: s
    integer :: ios
    call extract_valstring(line, key, s, success)
    if (.not. success) return
    read(s, *, iostat=ios) val
    success = (ios == 0)
  end subroutine get_value_int

  pure subroutine get_value_char(line, key, val, success)
    character(*), intent(in)  :: line, key
    character(*), intent(out) :: val
    logical,      intent(out) :: success
    character(len=:), allocatable :: s
    call extract_valstring(line, key, s, success)
    if (.not. success) return
    val = s(1:min(len(val), len_trim(s)))
    success = .true.
  end subroutine get_value_char

  pure subroutine extract_valstring(line, key, out, success)
    character(*), intent(in)  :: line, key
    character(len=:), allocatable, intent(out) :: out
    logical, intent(out) :: success
    integer :: p, i, j, n
    character(len=*), parameter :: delims = ' ,;'//achar(9)  ! space,
    ! comma, semicolon, tab
    
    success = .false.; out = ''
    n = len_trim(line)
    if (n == 0) return
    
    ! Find "key="
    p = index(line(1:n), trim(key)//'=')
    if (p == 0) return
    
    ! Start after "key="
    i = p + len_trim(key) + 1
    
    ! Skip leading spaces before the value
    do while (i <= n .and. line(i:i) == ' ')
       i = i + 1
    end do
    if (i > n) return
    
    ! Find next delimiter (or end)
    j = scan(line(i:n), delims)
    if (j == 0) then
       out = adjustl(line(i:n))
    else
       out = adjustl(line(i:i+j-2))
    end if
    success = (len_trim(out) > 0)
  end subroutine extract_valstring
  
END MODULE KEY_VALUE_PARSE
