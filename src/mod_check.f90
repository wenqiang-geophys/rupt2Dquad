module mod_check

  use mod_para
  use mod_types
  use netcdf

  implicit none

contains

subroutine check2(status,msg)
  implicit none
  integer, intent ( in) :: status
  character(len=*) :: msg

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))!,' in file ',__FILE__,' line ',__LINE__
    print *, trim(msg)
    stop 110
  end if
end subroutine

end module
