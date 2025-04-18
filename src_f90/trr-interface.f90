! 2017 Modefied by Si-Han Chen <chen.3262@osu.edu>
! Original XDR Fortran Interface Example Program is developed by James W. Barnett
! See https://github.com/wesbarnett/libgmxfort

MODULE TRR

  USE, INTRINSIC :: iso_c_binding, only: C_PTR, C_CHAR, C_FLOAT, C_INT
  IMPLICIT NONE
  private
  public xdrfile_open,xdrfile_close, read_trr_natoms, read_trr

  ! the data type located in libxdrfile
  type, public, bind(C) :: xdrfile
     type(C_PTR) :: fp, xdr
     character(kind=C_CHAR) :: mode
     integer(C_INT) :: buf1, buf1size, buf2, buf2size
  end type xdrfile

  ! interface with libxdrfile
  interface 
     
     integer(C_INT) function read_trr_natoms(filename,NATOMS) bind(C, name='read_trr_natoms')
       import
       character(kind=C_CHAR), intent(in) :: filename
       integer(C_INT), intent(out) :: NATOMS
     end function read_trr_natoms

     type(C_PTR) function xdrfile_open(filename,mode) bind(C, name='xdrfile_open')
       import
       character(kind=C_CHAR), intent(in) :: filename(*), mode(*)
     end function xdrfile_open

     integer(C_INT) function read_trr(xd,NATOMS,STEP,time,lambda,box,x,v,f) bind(C, name='read_trr')
       import
       type(xdrfile), intent(in) :: xd
       integer(C_INT), intent(out) :: NATOMS, STEP
       real(C_FLOAT), intent(out) :: lambda
       real(C_FLOAT), intent(out) :: time, box(*), x(*), v(*), f(*)
     end function read_trr

     integer(C_INT) function xdrfile_close(xd) bind(C,name='xdrfile_close')
       import
       type(xdrfile), intent(in) :: xd
     end function xdrfile_close

  end interface

END MODULE TRR


