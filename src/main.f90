
program main

use mod_para
use mod_gll
!use mod_nodes
!use mod_jacobi
use mod_types
use mod_mesh
use mod_geometry
use mod_wave
use mod_fault
!use mod_source
!use CMD_Progress
use mod_damp
use mod_pml
use mod_smooth
use mod_source
use mod_mpi
use mpi
use mod_exchange
use mod_io_fault
use mod_io_free
use mod_io_inter_s
use mod_io_inter_f
use mod_io_wave

implicit none

integer :: it,nt,irk
type(meshvar) :: mesh
type(buffvar) :: buff
!real(kind=rkind),allocatable,dimension(:,:,:,:) :: u,k1,k2,k3,k4
real(kind=rkind),allocatable,dimension(:,:,:) :: u,hu,mu,tu
real(kind=rkind),allocatable,dimension(:,:,:) :: au,hau,mau,tau
real(kind=rkind),allocatable,dimension(:,:,:) :: s,hs,ms,ts
real(kind=rkind) :: dt,tmax,minGLL,hmin,hmax,vpmax,vpmin,vsmax,vsmin,rhomax,rhomin
integer :: myrank,nproc
!character(len=80) :: filename
!real(kind=rkind), dimension(:,:), allocatable :: q_send,q_rec
!real(kind=rkind), dimension(:,:,:,:), allocatable :: qi
!integer,parameter :: dimens = 8
integer :: i,j,k,ie,tag,dest,req,req_r,c,je,face
integer :: is
integer :: irecv
real(kind=rkind) :: tmpu(Ngrid,Ngrid),tmpuface(Nfp)
real(kind=rkind),dimension(Nfp) :: xface,yface
integer :: sizes,ierr
!integer :: istat
integer :: tskip
character(len=128) :: filename
real(kind=rkind) :: li,lj,tmpu1,tmpu2,tmpu3
real(kind=rkind) :: amp,rise_time,t,r,rmin,srcx,srcy,svf,fc
real(kind=rkind) :: CFL,damp_s
integer :: src_idx,src_ie,src_i,src_j
real(kind=rkind) :: src_dist

integer :: fault_snap_skip
integer :: wave_snap_skip

fault_snap_skip = 1
wave_snap_skip = 1

!call init_mpi()
!call comm_rank(myrank)
call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)

mesh%nproc = nproc
mesh%rank = myrank
!print*,'myrank = ',myrank
!integer,parameter :: N = 4
!integer :: p1,p2
!real*8,allocatable :: x1(:),x2(:),wts1(:),wts2(:),G(:,:),D(:,:)
!
!p1 = N
!p2 = N+0
!
!allocate(x1(p1+1))
!allocate(x2(p2+1))
!allocate(wts1(p1+1))
!allocate(wts2(p2+1))
!allocate(G(p2+1,p1+1))
!allocate(D(p1+1,p1+1))
!
!call gll_nodes(p1,x1,wts1)
!call gll_nodes(p2,x2,wts2)
!call lagint(x1,x2,G,D)


!write(filename,'(a,i6.6)') 'mesh/meshVar',myrank
!call readMeshVar(mesh,'mesh/meshVar000000')
if (nproc == 1) then
    !call read_mesh(mesh,'mesh1')
    call readMeshVar(mesh,myrank)
else
    call readMeshVar(mesh,myrank)
end if
call build_geometry(mesh)

mesh_dir = 'data'
data_dir = 'data'

!    if (myrank==0) then
!        print*,'Filter='
!        do i = 1,NGLL
!            print*,(sngl(mesh%Filter(i,j)),j=1,NGLL)
!        enddo
!    end if
!stop 2

!if(myrank==0) then
!do i = 1,NGLL
!print*,sngl(mesh%D(i,:))
!end do
!end if
!stop 2

call write_mesh(mesh,myrank)

call init_damp(mesh,myrank)
call init_pml(mesh,myrank)

allocate( u(Np,mesh%nelem,8))
allocate(hu(Np,mesh%nelem,8))
allocate(mu(Np,mesh%nelem,8))
allocate(tu(Np,mesh%nelem,8))

allocate( au(Np,mesh%nelem,8))
allocate(hau(Np,mesh%nelem,8))
allocate(mau(Np,mesh%nelem,8))
allocate(tau(Np,mesh%nelem,8))

allocate( s(Nfp,nsurface,mesh%nelem))
allocate(hs(Nfp,nsurface,mesh%nelem))
allocate(ts(Nfp,nsurface,mesh%nelem))
allocate(ms(Nfp,nsurface,mesh%nelem))

call init_wave(mesh,u)
call fault_init(mesh)
!call write_wave(u,myrank,0)

call inter_s_io_init(mesh)
call inter_f_io_init(mesh)
call free_io_init(mesh)
call fault_io_init(mesh)
call wave_io_init(mesh)

! init mpi buffer
call init_buff(mesh,buff)
! mpi
!allocate(q_send(Nfp*mesh%mpi_ne*dimens,mesh%mpi_nn))
!allocate(q_rec(Nfp*mesh%mpi_ne*dimens,mesh%mpi_nn))
!allocate(qi(Nfp,dimens,mesh%mpi_ne,mesh%mpi_nn))

tmax = 500
!dt = 3.0*mesh%dtfactor
minGLL = mesh%xnode(2)-mesh%xnode(1)
!dt = mesh%dtfactor * minGLL / maxval(mesh%vp)
!print*,'dt=',dt
call minval_real_all(mesh%min_grid,hmin,CUSTOM_REAL)
call maxval_real_all(mesh%max_grid,hmax,CUSTOM_REAL)
call minval_real_all(mesh%min_rho,rhomin,CUSTOM_REAL)
call maxval_real_all(mesh%max_rho,rhomax,CUSTOM_REAL)
call minval_real_all(mesh%min_vp,vpmin,CUSTOM_REAL)
call maxval_real_all(mesh%max_vp,vpmax,CUSTOM_REAL)
call minval_real_all(mesh%min_vs,vsmin,CUSTOM_REAL)
call maxval_real_all(mesh%max_vs,vsmax,CUSTOM_REAL)
if (myrank==0) print*,'grid = ',sngl(hmin  ),'~',sngl(hmax  )
if (myrank==0) print*,'rho = ', sngl(rhomin),'~',sngl(rhomax)
if (myrank==0) print*,'vp = ',  sngl(vpmin ),'~',sngl(vpmax )
if (myrank==0) print*,'vs = ',  sngl(vsmin ),'~',sngl(vsmax )
!dt = 0.7*hmin/(2*order+1)/maxval(mesh%vp)
!dt = 0.7*hmin/(order**2)/maxval(mesh%vp)
!dt = 0.5*hmin/(2*order+1)/vpmax
!dt = 0.5*hmin*minGLL/vpmax

CFL = 0.8 ! O=4,Nsub=4
if (Order == 2) then
CFL = 2.0
end if
if (Order == 3) then
CFL = 2.5
CFL = 1.2
end if
if (Order == 4) then
CFL = 1.3
CFL = 0.5
CFL = 1.0/1d0
CFL = 1
end if
if (Order == 5) then
CFL = 2. * 0.5
end if
if (Order == 6) then
CFL = 1.
end if
if (order > 8) then
CFL = 0.3
end if

damp_s = 0.0
CFL = CFL / sqrt(1.0+2*damp_s)

if(myrank==0) print*,'Order=',Order
if(myrank==0) print*,'CFL=',sngl(CFL)

dt = CFL*hmin/Nsub/(2*order+1)/vpmax
#ifdef FD
dt = CFL*hmin/Ngrid/2.0/vpmax
#endif
mesh%dtfactor = dt
!#ifdef FD
!    dt = 0.2*hmin/order/2.0/maxval(mesh%vp)
!#endif
if (myrank==0) print*,'dt = ',sngl(dt),'Tmax = ',sngl(tmax)
!dt = 1.6d-3
nt = int(tmax/dt)+1
if (myrank==0) print*,'nt = ',nt
!print*,'rmin=',minval(mesh%xnode(2:order+1)-mesh%xnode(1:order))
!print*,'rmax=',maxval(mesh%xnode(2:order+1)-mesh%xnode(1:order))
!print*,'dt = ',dt
!print*,'dt = ',dt
tskip = int(0.5/dt)
if(myrank==0) print*,'tskip = ',tskip

wave_snap_skip = tskip

do i = 1,mesh%nrecv
    write(filename,'(a,i6.6,a)') 'data/recv_id',mesh%recv_fid(i),'.txt'
    open(100+i,file=trim(filename))
enddo
do i = 1,mesh%body_nrecv
    write(filename,'(a,i6.6,a)') 'data/body_recv_id',mesh%body_recv_fid(i),'.txt'
    open(900+i,file=trim(filename))
enddo

!do i = 1,mesh%nelem
!    do is = 1,nsurface
!        if (mesh%bctype(is,i)>=BC_FAULT) then
!write(filename,'(a,i6.6,a,i6.6,a)') 'data/fault_mpi',myrank,'_',i,'.txt'
!            open(2000+myrank*1000+i,file=trim(filename))
!            if (is == 1) then
!                xface = mesh%vx(:,1,i)
!                yface = mesh%vy(:,1,i)
!            else if (is == 2) then
!                xface = mesh%vx(NGLL,:,i)
!                yface = mesh%vy(NGLL,:,i)
!            else if (is == 3) then
!                xface = mesh%vx(:,NGLL,i)
!                yface = mesh%vy(:,NGLL,i)
!            else if (is == 4) then
!                xface = mesh%vx(1,:,i)
!                yface = mesh%vy(1,:,i)
!            end if
!            write(2000+myrank*1000+i,*) xface
!            write(2000+myrank*1000+i,*) yface
!            write(2000+myrank*1000+i,*) mesh%nx(:,is,i)
!            write(2000+myrank*1000+i,*) mesh%ny(:,is,i)
!        end if
!    end do
!enddo

write(filename,'(a,i6.6,a)') 'data/fault_mpi',myrank,'.txt'
open(1000+myrank,file=trim(filename))
write(1000+myrank,*) mesh%fault_buffer(:,:,1)
write(1000+myrank,*) mesh%fault_buffer(:,:,2)
write(1000+myrank,*) mesh%fault_buffer(:,:,3)
write(1000+myrank,*) mesh%fault_buffer(:,:,4)

write(filename,'(a,i6.6,a)') 'data/stress_mpi',myrank,'.txt'
open (3000+myrank,file=trim(filename))
write(3000+myrank,*) mesh%fault_buffer(:,:,1)
write(3000+myrank,*) mesh%fault_buffer(:,:,2)
write(3000+myrank,*) mesh%fault_buffer(:,:,3)
write(3000+myrank,*) mesh%fault_buffer(:,:,4)
!close(1000+myrank)

write(filename,'(a,i6.6,a)') 'data/surface_mpi',myrank,'.txt'
open (2000+myrank,file=trim(filename))
write(2000+myrank,*) mesh%surface_buffer(:,:,1)
write(2000+myrank,*) mesh%surface_buffer(:,:,2)
write(2000+myrank,*) mesh%surface_buffer(:,:,3)
write(2000+myrank,*) mesh%surface_buffer(:,:,4)

call MPI_Barrier(MPI_COMM_WORLD,ierr)

!do i = 1,mesh%body_nrecv
!    print*,mesh%body_recv_fid(i),mesh%body_recv_refx(i),mesh%body_recv_refy(i)
!end do
!stop 2


!src_idx  = minloc(sqrt((mesh%vx(:)-0d0)**2+(mesh%vy(:)-0d0)**2))
!src_dist = minval(sqrt((mesh%vx(:)-0d0)**2+(mesh%vy(:)-0d0)**2))


rmin = 1d30
do ie = 1,mesh%nelem
    do j = 1,NGLL
        do i = 1,NGLL
            r = sqrt((mesh%vx(i,j,ie)-0d0)**2+(mesh%vy(i,j,ie)-0d0)**2)
            if (r < rmin) then
                rmin = r
                src_ie = ie
                src_i = i
                src_j = j
                src_dist = r
            end if
        end do
    end do
end do


mesh%current_timestep = 0
mesh%current_time = 0d0
s = 0
tu = 0
ts = 0
mu = 0
ms = 0

u = 0
! init wave
!if (.true.) then
if (.true.) then
do ie = 1,mesh%nelem
    !if (mesh%elemtype(ie) == ELEM_SOLID) then
    !print*,"adding source"
    do j = 1,Ngrid
        do i = 1,Ngrid
            srcx = -0e3
            srcy = -3e3
            amp = 1.0
            r = sqrt((mesh%vx(i,j,ie)-srcx)**2 + (mesh%vy(i,j,ie)-srcy)**2)
            if ( r < 3e3 ) then
                u(i+(j-1)*Ngrid,ie,3) = amp * exp(-3.0*(r/3e3)**2)
                u(i+(j-1)*Ngrid,ie,4) = amp * exp(-3.0*(r/3e3)**2)
            end if
        end do
    end do
    !end if
end do
end if

do it = 1,nt
    mesh%current_timestep = it
    mesh%current_time = mesh%current_timestep * dt
    t = it*dt

    ! add source
    !!!rise_time = 0.1
    !!!if (t< rise_time) then
    !!!    amp = 0.5*(1-cos( t/rise_time * PI * 2.0))
    !!!else
    !!!    amp = 0
    !!!end if

    !do ie = 1,mesh%nelem
    !    do j = 1,NGLL
    !        do i = 1,NGLL

    !            r = sqrt((mesh%vx(i,j,ie)-0)**2 + (mesh%vy(i,j,ie)-(-2.5e3))**2)
    !            if ( r < 2e3 ) then
    !                u(i+j*NGLL,ie,3) = u(i+j*NGLL,ie,3) + amp !* exp(-3.0*(r/3e3)**2)
    !                u(i+j*NGLL,ie,4) = u(i+j*NGLL,ie,4) + amp !* exp(-3.0*(r/3e3)**2)
    !            end if
    !        end do
    !    end do
    !end do
    ! end add source

    ! add point source
    !!!amp = 1.0e9
    !!!fc = 0.4
    !!!svf = svf_ricker(it*dt,fc,1.2d0/fc)
    !!!if (.false.) then
    !!!    mesh%bctype = BC_IN
    !!!    if( src_dist < hmin) then
    !!!        i = src_i
    !!!        j = src_j
    !!!        ie = src_ie
    !!!        u(i+(j-1)*NGLL,ie,3) = u(i+(j-1)*NGLL,ie,3) + amp * svf
    !!!        u(i+(j-1)*NGLL,ie,4) = u(i+(j-1)*NGLL,ie,4) + amp * svf
    !!!    end if
    !!!end if

    mu = u
    ms = s
    tu = 0
    ts = 0
    do irk = 1,5

        call exchange_data(mesh,u,buff)

        !mu = u
        !mu(:,:,3:5) = mu(:,:,3:5) + 0.1*dt*hu(:,:,3:5)
        call rhs(mesh,u,buff%qi,hu)
        !hs = mesh%sliprate

        !do i = 1,5
        !    hu(:,:,i) = hu(:,:,i) * reshape(mesh%damp,(/Np,mesh%Nelem/))
        !end do
        !if (.true.) then
        !do ie = 1,mesh%nelem
        !    if (mesh%elemtype(ie) == ELEM_SOLID) then
        !    !print*,"adding source"
        !    do j = 1,Ngrid
        !        do i = 1,Ngrid
        !            srcx = -50e3
        !            srcy = -50e3
        !            r = sqrt((mesh%vx(i,j,ie)-srcx)**2 + (mesh%vy(i,j,ie)-srcy)**2)
        !            if ( r < 2e3 ) then
        !                hu(i+(j-1)*Ngrid,ie,3) = hu(i+(j-1)*Ngrid,ie,3) + amp * exp(-3.0*(r/2e3)**2)
        !                hu(i+(j-1)*Ngrid,ie,4) = hu(i+(j-1)*Ngrid,ie,4) + amp * exp(-3.0*(r/2e3)**2)
        !            end if
        !        end do
        !    end do
        !    end if
        !end do
        !end if

        tu = rk4a(irk)*tu + dt*hu
        !ts = rk4a(irk)*ts + dt*hs
        mesh%tslip = rk4a(irk)*mesh%tslip + dt*mesh%sliprate

        u = u + rk4b(irk)*tu
        !s = s + rk4b(irk)*ts
        mesh%slip = mesh%slip + rk4b(irk)*mesh%tslip

        ! artificial damping
        !u(:,:,1:2) = u(:,:,1:2) + 0.3*dt*hu(:,:,1:2)
        u(:,:,6:8) = damp_s*dt*hu(:,:,3:5) ! save strain rate

        !mesh%slip = s

        !call smooth_fault(mesh,u)


    end do ! irk

    do i = 1,5
        u(:,:,i) = u(:,:,i) * reshape(mesh%damp,(/Np,mesh%Nelem/))
    end do

    !call smooth_fault(mesh,u)

    mesh%eta = mesh%eta + mesh%deta * dt

    do irecv = 1,mesh%nrecv
        tmpu1 = 0
        tmpu2 = 0
        tmpu3 = 0
        do i = 1,Ngrid
            li=lagrange_basis(Ngrid,i,mesh%recv_refx(irecv),mesh%xnode)
            tmpu1 = tmpu1 + mesh%recv_buffer(irecv,i,1)*li 
            tmpu2 = tmpu2 + mesh%recv_buffer(irecv,i,2)*li 
            tmpu3 = tmpu3 + mesh%recv_buffer(irecv,i,3)*li 
        end do

        write(100+irecv,*) mesh%current_time,tmpu1,tmpu2,tmpu3
        !write(100+irecv,*) mesh%current_time,mesh%recv_buffer(irecv,1:3)

    end do
    do irecv = 1,mesh%body_nrecv
        tmpu1 = 0
        tmpu2 = 0
        do i = 1,Ngrid
            do j = 1,Ngrid
                li=lagrange_basis(Ngrid,i,mesh%body_recv_refx(irecv),mesh%xnode)
                lj=lagrange_basis(Ngrid,j,mesh%body_recv_refy(irecv),mesh%xnode)
                tmpu1 = tmpu1 + u(i+(j-1)*Ngrid,mesh%body_recv_ie(irecv),1)*li*lj 
                tmpu2 = tmpu2 + u(i+(j-1)*Ngrid,mesh%body_recv_ie(irecv),2)*li*lj 
            end do
        end do
        !tmpu1 = tmpu1 / mesh%rho(mesh%body_recv_ie(irecv))
        !tmpu2 = tmpu2 / mesh%rho(mesh%body_recv_ie(irecv))
        write(900+irecv,*) mesh%current_time,&
            tmpu1,tmpu2
            !u(mesh%body_recv_i(irecv),mesh%body_recv_j(irecv),mesh%body_recv_ie(irecv),1:2)
    end do

    !do i = 1,mesh%nelem
    !    do is = 1,nsurface
    !        if (mesh%bctype(is,i)>=BC_FAULT) then
    !            if (is == 1) then
    !                tmpuface = u(:,1,i,1)
    !            else if (is == 2) then
    !                tmpuface = u(NGLL,:,i,1)
    !            else if (is == 3) then
    !                tmpuface = u(:,NGLL,i,1)
    !            else if (is == 4) then
    !                tmpuface = u(1,:,i,1)
    !            end if
    !            write(2000+myrank*1000+i,*) tmpuface
    !        end if
    !    end do
    !enddo
    if (mod(it-1,5) == 0) then
    k = 0
    do i = 1,mesh%nelem
        do is = 1,nsurface
            if (mesh%bctype(is,i)>=BC_FAULT) then
                k = k + 1
                mesh%fault_buffer(:,k,1) = mesh%sliprate(:,is,i)
                !mesh%fault_buffer(:,k,1) = (s(:,is,i)-ms(:,is,i))/dt
            end if
        end do
    end do
    write(1000+myrank,*) mesh%fault_buffer(:,:,1)

    k = 0
    do i = 1,mesh%nelem
        do is = 1,nsurface
            if (mesh%bctype(is,i)>=BC_FAULT) then
                k = k + 1
                mesh%fault_buffer(:,k,1) = mesh%stress(:,is,i)
                !mesh%fault_buffer(:,k,1) = (s(:,is,i)-ms(:,is,i))/dt
            end if
        end do
    end do
    write(3000+myrank,*) mesh%fault_buffer(:,:,1)

    k = 0
    do i = 1,mesh%nelem
        do is = 1,nsurface
            if (mesh%bctype(is,i)==BC_FREE .or. mesh%bctype(is,i)==BC_FREE_G) then
                k = k + 1
                mesh%surface_buffer(:,k,1) = mesh%eta(:,is,i)
            end if
        end do
    end do
    write(2000+myrank,*) mesh%surface_buffer(:,:,1)
    end if

    if (mod(it-1,tskip)==0) then
        !if (myrank==0) print*,it,'/',nt
    !    !if (myrank==0) print*,nt

        !if(myrank==0) print*,'writing data/wave @ it = ',it
        !call write_wave(u,myrank,1)
    end if


    if (mod(it-1,fault_snap_skip) == 0) then
        call inter_s_io_save(mesh,u,(it-1)/fault_snap_skip+1)
        call inter_f_io_save(mesh,u,(it-1)/fault_snap_skip+1)
    end if
    if (mod(it-1,fault_snap_skip) == 0) then
        call free_io_save(mesh,u,(it-1)/fault_snap_skip+1)
    end if
    if (mod(it-1,fault_snap_skip) == 0) then
        call fault_io_save(mesh,(it-1)/fault_snap_skip+1)
    end if
    if (mod(it-1,wave_snap_skip) == 0) then
        if(myrank==0) print*,'writing data/wave_mpi @ it = ',it
        call wave_io_save(mesh,u,(it-1)/wave_snap_skip+1)
    end if
enddo ! it

call inter_s_io_end(mesh)
call inter_f_io_end(mesh)
call free_io_end(mesh)
call fault_io_end(mesh)
call wave_io_end(mesh)

call MPI_Barrier(MPI_COMM_WORLD,ierr)
do i = 1,mesh%nrecv
    close(100+i)
end do
do i = 1,mesh%body_nrecv
    close(900+i)
end do

call MPI_Finalize(ierr)


end program
