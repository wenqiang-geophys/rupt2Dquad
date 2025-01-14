module mod_io_wave

  use mod_para
  use mod_types
  use netcdf
  use mod_check

  implicit none

  integer :: wave_ncid
  integer :: wave_varid(20)

  !real(kind=rkind),allocatable,dimension(:,:) :: x,y
  real(kind=rkind),allocatable,dimension(:,:,:) :: Vx,Vy,exx,eyy,exy,eta

contains

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
  allocate(exx(Ngrid,Ngrid,mesh%nelem))
  allocate(eyy(Ngrid,Ngrid,mesh%nelem))
  allocate(exy(Ngrid,Ngrid,mesh%nelem))
  allocate(eta(Ngrid,Ngrid,mesh%nelem))
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
  ierr = nf90_def_var(wave_ncid, "exx", NF90_FLOAT, dimid(1:4), wave_varid(6))
  call check2(ierr,'def_var exx')
  ierr = nf90_def_var(wave_ncid, "eyy", NF90_FLOAT, dimid(1:4), wave_varid(7))
  call check2(ierr,'def_var eyy')
  ierr = nf90_def_var(wave_ncid, "exy", NF90_FLOAT, dimid(1:4), wave_varid(8))
  call check2(ierr,'def_var exy')

  ierr = nf90_def_var(wave_ncid, "damage", NF90_FLOAT, dimid(1:4), wave_varid(5))
  call check2(ierr,'def_var eta')

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

subroutine wave_io_save(mesh,wave,it)
  implicit none
  type(meshvar) :: mesh
  type(wavevar) :: wave
  !real(kind=rkind) :: u(:,:,:)
  integer :: i, j, ie, is, it, ierr
  integer,dimension(4) :: start,cnt,stride
  do ie = 1,mesh%nelem
    do j = 1,Ngrid
      do i = 1,Ngrid
        Vx(i,j,ie) = sngl( wave%u(i+(j-1)*Ngrid,ie,1)/mesh%rho(ie) )
        Vy(i,j,ie) = sngl( wave%u(i+(j-1)*Ngrid,ie,2)/mesh%rho(ie) )
        exx(i,j,ie) = sngl( wave%u(i+(j-1)*Ngrid,ie,3) )
        eyy(i,j,ie) = sngl( wave%u(i+(j-1)*Ngrid,ie,4) )
        exy(i,j,ie) = sngl( wave%u(i+(j-1)*Ngrid,ie,5) )
        eta(i,j,ie) = sngl( wave%damage(i,j,ie) )
      enddo
    enddo
  enddo

  start=(/1,1,1,it+0/); cnt=(/Ngrid,Ngrid,mesh%nelem,1/); stride=(/1,1,1,1/)
  ierr = nf90_put_var(wave_ncid,wave_varid(1),Vx,start,cnt,stride)
  call check2(ierr,'put_var Vx')
  ierr = nf90_put_var(wave_ncid,wave_varid(2),Vy,start,cnt,stride)
  call check2(ierr,'put_var Vy')
  ierr = nf90_put_var(wave_ncid,wave_varid(6),exx,start,cnt,stride)
  call check2(ierr,'put_var exx')
  ierr = nf90_put_var(wave_ncid,wave_varid(7),eyy,start,cnt,stride)
  call check2(ierr,'put_var eyy')
  ierr = nf90_put_var(wave_ncid,wave_varid(8),exy,start,cnt,stride)
  call check2(ierr,'put_var exy')
  ierr = nf90_put_var(wave_ncid,wave_varid(5),eta,start,cnt,stride)
  call check2(ierr,'put_var eta')

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
