module mod_io_free

  use mod_para
  use mod_types
  use netcdf
  use mod_check

  implicit none

  integer :: free_ncid
  integer :: free_varid(50)

contains

subroutine free_io_init(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: myrank
  integer :: i, ie, ief, is, ierr
  character(len=128) :: filename
  real(kind=rkind),allocatable,dimension(:,:) :: x,y,z,nx,ny,nz
  real(kind=rkind),allocatable,dimension(:,:) :: mx,my,mz,lx,ly,lz
  integer,allocatable,dimension(:) :: free_id
  real(kind=rkind),allocatable,dimension(:,:,:) :: fltvars
  integer :: dimid(50),varid(50)
  real(rkind) :: xface(Nfp),yface(Nfp)

  if (mesh%nface_free==0) return

  myrank = mesh%rank

  write(filename,'(a,a,i6.6,a)') trim(data_dir),'/free_mpi',myrank,'.nc'
  !print*,filename
  allocate(x(Nfp,mesh%nface_free))
  allocate(y(Nfp,mesh%nface_free))
  allocate(nx(Nfp,mesh%nface_free))
  allocate(ny(Nfp,mesh%nface_free))
  !allocate(fltvars(Nfp,mesh%nface_free,7))
  i = 0
  do ie = 1,mesh%nelem
  !do ief = 1,mesh%nelem_free
    !ie = mesh%free2wave(ief)
    do is = 1,Nfaces
      if ( mesh%bctype(is,ie) == BC_FREE) then
        i = i + 1
        if (is == 1) then
            xface = mesh%vx(:,1,ie)
            yface = mesh%vy(:,1,ie)
        else if (is == 2) then
            xface = mesh%vx(Ngrid,:,ie)
            yface = mesh%vy(Ngrid,:,ie)
        else if (is == 3) then
            xface = mesh%vx(:,Ngrid,ie)
            yface = mesh%vy(:,Ngrid,ie)
        else if (is == 4) then
            xface = mesh%vx(1,:,ie)
            yface = mesh%vy(1,:,ie)
        end if
        x(1:Nfp,i) = xface
        y(1:Nfp,i) = yface
        nx(1:Nfp,i) = mesh%nx(:,is,ie)
        ny(1:Nfp,i) = mesh%ny(:,is,ie)
      endif
    enddo
  enddo

  ! start to create
  ierr = nf90_create(filename, NF90_CLOBBER, free_ncid)
  call check2(ierr,'nf90_create free')

  ! define dimensions
  ierr = nf90_def_dim(free_ncid, "Nfp", Nfp, dimid(1))
  call check2(ierr,'def_dim Nfp')
  ierr = nf90_def_dim(free_ncid, "Nfree", mesh%nface_free, dimid(2))
  call check2(ierr,'def_dim Nflt')
  ierr = nf90_def_dim(free_ncid, "Nt", nf90_unlimited, dimid(3))
  call check2(ierr,'def_dim Nt')

  ! define variables
  ierr = nf90_def_var(free_ncid, "x" , NF90_DOUBLE, dimid(1:2), varid(1))
  call check2(ierr,'def_var x')
  ierr = nf90_def_var(free_ncid, "y" , NF90_DOUBLE, dimid(1:2), varid(2))
  call check2(ierr,'def_var y')
  ierr = nf90_def_var(free_ncid, "nx", NF90_DOUBLE, dimid(1:2), varid(3))
  call check2(ierr,'def_var nx')
  ierr = nf90_def_var(free_ncid, "ny", NF90_DOUBLE, dimid(1:2), varid(4))
  call check2(ierr,'def_var ny')

  ierr = nf90_def_var(free_ncid, "time" , NF90_DOUBLE, dimid(3), free_varid(1))
  call check2(ierr,'def_var time')

  ierr = nf90_def_var(free_ncid, "Vx",   NF90_FLOAT, dimid(1:3), free_varid(2))
  call check2(ierr,'def_var Vx')
  ierr = nf90_def_var(free_ncid, "Vy",  NF90_FLOAT, dimid(1:3), free_varid(3))
  call check2(ierr,'def_var Vy')
  ! end of define
  ierr = nf90_enddef(free_ncid)
  call check2(ierr,'enddef')

  ! put variables
  ierr = nf90_put_var(free_ncid, varid(1), x)
  call check2(ierr,'put_var x')
  ierr = nf90_put_var(free_ncid, varid(2), y)
  call check2(ierr,'put_var y')
  ierr = nf90_put_var(free_ncid, varid(3), nx)
  call check2(ierr,'put_var nx')
  ierr = nf90_put_var(free_ncid, varid(4), ny)
  call check2(ierr,'put_var ny')

  deallocate(x,y,nx,ny)
end subroutine

subroutine free_io_save(mesh,u,it)
  implicit none
  type(meshvar) :: mesh
  real(kind=rkind) :: u(:,:,:)
  real :: Vx1(Ngrid,Ngrid)
  real :: Vy1(Ngrid,Ngrid)
  integer :: i, ie, ief, is, it, ierr
  integer :: j
  real,allocatable,dimension(:,:) :: Vx,Vy
  integer,dimension(3) :: start,cnt,stride

  if (mesh%nface_free==0) return
  allocate(Vx(Nfp,mesh%nface_free))
  allocate(Vy(Nfp,mesh%nface_free))

  i = 0
  do ie = 1,mesh%nelem
  !do ief = 1,mesh%nelem_free
    !ie = mesh%free2wave(ief)
    do is = 1,Nfaces
      if ( mesh%bctype(is,ie) == BC_FREE) then
        i = i + 1
        Vx1 = sngl(reshape(u(:,ie,1),(/Ngrid,Ngrid/)))
        Vy1 = sngl(reshape(u(:,ie,2),(/Ngrid,Ngrid/)))
        if (is==1) then
            Vx(1:Nfp,i) = Vx1(:,1)
            Vy(1:Nfp,i) = Vy1(:,1)
        elseif (is==2) then
            Vx(1:Nfp,i) = Vx1(Ngrid,:)
            Vy(1:Nfp,i) = Vy1(Ngrid,:)
        elseif (is==3) then
            Vx(1:Nfp,i) = Vx1(:,Ngrid)
            Vy(1:Nfp,i) = Vy1(:,Ngrid)
        else if (is==4) then
            Vx(1:Nfp,i) = Vx1(1,:)
            Vy(1:Nfp,i) = Vy1(1,:)
        endif
        !Vx(i,j,ie) = sngl( u(i+(j-1)*Ngrid,ie,1)/mesh%rho(ie) )
        !Vy(i,j,ie) = sngl( u(i+(j-1)*Ngrid,ie,2)/mesh%rho(ie) )
        do j=1,Nfp
            Vx(j,i)=Vx(j,i)/mesh%rho(ie)
            Vy(j,i)=Vy(j,i)/mesh%rho(ie)
        enddo

      endif
    enddo
  enddo
  start=(/1,1,it+0/); cnt=(/Nfp,mesh%nface_free,1/); stride=(/1,1,1/)
  ierr = nf90_put_var(free_ncid,free_varid(2),Vx,start,cnt,stride)
  call check2(ierr,'put_var Vx')
  ierr = nf90_put_var(free_ncid,free_varid(3),Vy,start,cnt,stride)
  call check2(ierr,'put_var Vy')

  ierr = nf90_put_var(free_ncid,free_varid(1), &
      (/mesh%current_time/),(/it/),(/1/),(/1/))
  call check2(ierr,'put_var time')

  ierr = nf90_sync(free_ncid)

end subroutine

subroutine free_io_end(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: ierr
  if (mesh%nface_free==0) return
  ierr = nf90_close(free_ncid)
  call check2(ierr,'nf90_close free')
end subroutine

end module
