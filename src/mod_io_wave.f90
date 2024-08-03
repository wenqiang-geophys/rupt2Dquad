module mod_io_wave

  use mod_para
  use mod_types
  use netcdf

  implicit none

  integer :: wave_ncid
  integer :: wave_varid(20)

  !real(kind=rkind),allocatable,dimension(:,:) :: x,y
  real(kind=rkind),allocatable,dimension(:,:,:) :: Vx,Vy

contains

subroutine check(status)
  implicit none
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status)),' in file ',__FILE__,' line ',__LINE__
    stop "Stopped"
  end if
end subroutine

subroutine check2(status,msg)
  implicit none
  integer, intent ( in) :: status
  character(len=*) :: msg

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    print *, trim(msg)
    stop 110
  end if
end subroutine

subroutine wave_io_init(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: myrank
  integer :: i, j, ie, is, ierr
  character(len=128) :: filename
  real(kind=rkind),allocatable,dimension(:,:,:) :: x,y!,PGVx,PGVy,PGVz
  !real(kind=rkind),allocatable,dimension(:,:,:) :: vars
  integer :: dimid(20),varid(20)
  real(rkind) :: xface(Nfp),yface(Nfp)

  myrank = mesh%rank

  write(filename,'(a,a,i6.6,a)') trim(data_dir),'/wave_mpi',myrank,'.nc'
  !print*,filename
  allocate( x(Ngrid,Ngrid,mesh%nelem))
  allocate( y(Ngrid,Ngrid,mesh%nelem))
  allocate(Vx(Ngrid,Ngrid,mesh%nelem))
  allocate(Vy(Ngrid,Ngrid,mesh%nelem))
  do ie = 1,mesh%nelem
    do i = 1,Ngrid
      do j = 1,Ngrid
         x(i,j,ie) = mesh%vx(i,j,ie)
         y(i,j,ie) = mesh%vy(i,j,ie)
      enddo
    enddo
  enddo

  ! start to create
  ierr = nf90_create(filename, NF90_CLOBBER, wave_ncid)
  call check2(ierr,'nf90_create wave')

  ! define dimensions
  ierr = nf90_def_dim(wave_ncid, "Nx", Ngrid, dimid(1))
  call check2(ierr,'def_dim Ngrid')
  ierr = nf90_def_dim(wave_ncid, "Ny", Ngrid, dimid(2))
  call check2(ierr,'def_dim Ngrid')
  ierr = nf90_def_dim(wave_ncid, "Nelem", mesh%nelem, dimid(3))
  call check2(ierr,'def_dim Nelem')
  ierr = nf90_def_dim(wave_ncid, "Nt", nf90_unlimited, dimid(4))
  call check2(ierr,'def_dim Nt')

  ! define variables
  ierr = nf90_def_var(wave_ncid, "x" , NF90_DOUBLE, dimid(1:3), varid(1))
  call check2(ierr,'def_var x')
  ierr = nf90_def_var(wave_ncid, "y" , NF90_DOUBLE, dimid(1:3), varid(2))
  call check2(ierr,'def_var y')

  ierr = nf90_def_var(wave_ncid, "time", NF90_DOUBLE, dimid(4), wave_varid(4))
  call check2(ierr,'def_var time')

  ierr = nf90_def_var(wave_ncid, "Vx", NF90_FLOAT, dimid(1:4), wave_varid(1))
  call check2(ierr,'def_var Vx')
  ierr = nf90_def_var(wave_ncid, "Vy", NF90_FLOAT, dimid(1:4), wave_varid(2))
  call check2(ierr,'def_var Vy')

  ! end of define
  ierr = nf90_enddef(wave_ncid)
  call check2(ierr,'enddef')

  ! put variables
  ierr = nf90_put_var(wave_ncid, varid(1), x)
  call check2(ierr,'put_var x')
  ierr = nf90_put_var(wave_ncid, varid(2), y)
  call check2(ierr,'put_var y')

  !!!ierr = nf90_close(wave_ncid)
  deallocate(x,y)
  !deallocate(vars)
end subroutine

subroutine wave_io_save(mesh,u,it)
  implicit none
  type(meshvar) :: mesh
  real(kind=rkind) :: u(:,:,:)
  integer :: i, j, ie, is, it, ierr
  integer,dimension(4) :: start,cnt,stride
  do ie = 1,mesh%nelem
    do j = 1,Ngrid
      do i = 1,Ngrid
        Vx(i,j,ie) = sngl( u(i+(j-1)*Ngrid,ie,1)/mesh%rho(ie) )
        Vy(i,j,ie) = sngl( u(i+(j-1)*Ngrid,ie,2)/mesh%rho(ie) )
      enddo
    enddo
  enddo

  start=(/1,1,1,it+0/); cnt=(/Ngrid,Ngrid,mesh%nelem,1/); stride=(/1,1,1,1/)
  ierr = nf90_put_var(wave_ncid,wave_varid(1),Vx,start,cnt,stride)
  call check2(ierr,'put_var Vx')
  ierr = nf90_put_var(wave_ncid,wave_varid(2),Vy,start,cnt,stride)
  call check2(ierr,'put_var Vy')

  ierr = nf90_put_var(wave_ncid,wave_varid(4), &
      (/mesh%current_time/),(/it/),(/1/),(/1/))
  call check2(ierr,'put_var time')

  ierr = nf90_sync(wave_ncid)

  !deallocate(Vx,Vy,Vz)
end subroutine

subroutine wave_io_end(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: ierr
  ierr = nf90_close(wave_ncid)
  call check2(ierr,'nf90_close wave')

  deallocate(Vx,Vy)
end subroutine

end module
