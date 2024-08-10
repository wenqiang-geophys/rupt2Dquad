module mod_types

use mod_para

implicit none

public :: meshvar

type :: meshvar
    integer :: current_timestep
    integer :: nelem ! number of elements
    integer :: ncoord ! number of points
    integer :: nglob
    integer :: rank,nproc
    integer :: nfree
    integer :: nfault
    integer :: nelem_fault
    integer :: nface_fault
    integer :: nface_free
    integer :: nface_inter_s
    integer :: nface_inter_f
    integer :: mpi_nn
    integer :: mpi_ne
    integer :: mpi_nemax
    integer :: pinterfaces
    integer :: nrecv
    integer :: body_nrecv
    real(kind=rkind) :: dtfactor
    real(kind=rkind) :: current_time
    real(kind=rkind) :: max_grid
    real(kind=rkind) :: min_grid
    real(kind=rkind) :: max_vp
    real(kind=rkind) :: min_vp
    real(kind=rkind) :: max_vs
    real(kind=rkind) :: min_vs
    real(kind=rkind) :: max_rho
    real(kind=rkind) :: min_rho

    ! coord, elem, and connectivity
    real(kind=rkind),dimension(:,:),pointer :: coord => null()
    integer,dimension(:,:),pointer :: elem => null()
    integer,dimension(:,:),pointer :: neigh => null()
    integer,dimension(:,:),pointer :: face => null() ! is always 0
    !integer,dimension(:,:),pointer :: face_flag => null() ! save when fault face flag < 0
    integer,dimension(:,:),pointer :: direction => null()
    integer,dimension(:,:),pointer :: bctype => null()
    integer,dimension(:,:),pointer :: fluxtype => null()
    integer,dimension(:),pointer :: elemtype => null()
    integer,dimension(:),pointer :: ispml => null()
    !integer,dimension(:,:,:),pointer :: vmapM => null()
    !integer,dimension(:,:,:),pointer :: vmapP => null()
    integer :: vmapM(Nfp,Nfaces)
    integer,dimension(:,:,:),pointer :: vmapP => null()

    integer,dimension(:),pointer :: mpi_neighbor => null()
    integer,dimension(:,:,:),pointer :: mpi_connection => null()
    integer,dimension(:,:),pointer :: mpi_ibool => null()
    integer,dimension(:,:,:),pointer :: mpi_interface => null()
    ! mpi media
    real(kind=rkind),dimension(:,:),pointer :: mpi_rho => null()
    real(kind=rkind),dimension(:,:),pointer :: mpi_vp => null()
    real(kind=rkind),dimension(:,:),pointer :: mpi_vs => null()

    integer,dimension(:),pointer :: free_node => null()
    integer,dimension(:),pointer :: fault_node => null()
    integer,dimension(:),pointer :: fault_elem => null()
    integer,dimension(:),pointer :: fault_face => null()
    integer,dimension(:),pointer :: fault2wave => null()
    integer,dimension(:),pointer :: wave2fault => null()
    real(kind=rkind),dimension(:,:,:),pointer :: tau_n => null()
    real(kind=rkind),dimension(:,:,:),pointer :: tau_0 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: stress => null()
    real(kind=rkind),dimension(:,:,:),pointer :: Slip => null()
    real(kind=rkind),dimension(:,:,:),pointer :: mSlip => null()
    real(kind=rkind),dimension(:,:,:),pointer :: tSlip => null()
    real(kind=rkind),dimension(:,:,:),pointer :: Sliprate => null()
    real(kind=rkind),dimension(:,:,:),pointer :: mu_s => null()
    real(kind=rkind),dimension(:,:,:),pointer :: mu_d => null()
    real(kind=rkind),dimension(:,:,:),pointer :: Dc => null()
    real(kind=rkind),dimension(:,:,:),pointer :: C0 => null()
    logical,dimension(:,:,:),pointer :: isslip => null()

    ! surface height eta
    real(kind=rkind),dimension(:,:,:),pointer :: eta => null()
    real(kind=rkind),dimension(:,:,:),pointer :: deta => null()

    ! matrices
    !real(kind=rkind),dimension(Np,Np) :: Dr,Ds,Drw,Dsw
    real(kind=rkind) :: xnode(Ngrid)
    real(kind=rkind) :: fd_D(Ngrid,Ngrid),fd_Dt(Ngrid,Ngrid)
    real(kind=rkind) :: G(order2+1,order+1),D(order2+1,order+1),W(order2+1,order2+1),Dt(order+1,order2+1)
    real(kind=rkind) :: Dp(NGLL,NGLL),Dm(NGLL,NGLL)
    real(kind=rkind),dimension(order+1,order+1) :: Mass,invMass,Vdm,invVdm,Filter
    !integer :: Fmask(NGLL,nsurface)
    !real(kind=rkind) :: LIFT(Np,3*NGLL)
    integer :: iglob(NGLL,Nsub)
    real(kind=rkind) :: globInvM(Ngrid)

    !real(kind=rkind),dimension(:),pointer :: r => null()
    !real(kind=rkind),dimension(:),pointer :: s => null()
    real(kind=rkind),dimension(:,:,:),pointer :: vx => null()
    real(kind=rkind),dimension(:,:,:),pointer :: vy => null()
    real(kind=rkind),dimension(:,:,:),pointer :: rx => null()
    real(kind=rkind),dimension(:,:,:),pointer :: ry => null()
    real(kind=rkind),dimension(:,:,:),pointer :: sx => null()
    real(kind=rkind),dimension(:,:,:),pointer :: sy => null()
    real(kind=rkind),dimension(:,:,:),pointer :: detJ => null()
    real(kind=rkind),dimension(:,:,:),pointer :: nx => null()
    real(kind=rkind),dimension(:,:,:),pointer :: ny => null()
    real(kind=rkind),dimension(:,:,:),pointer :: norm_n => null()
    real(kind=rkind),dimension(:,:,:),pointer :: damp => null()

    ! PML
    real(kind=rkind),dimension(:,:,:),pointer :: pax => null()
    real(kind=rkind),dimension(:,:,:),pointer :: pbx => null()
    real(kind=rkind),dimension(:,:,:),pointer :: pdx => null()
    real(kind=rkind),dimension(:,:,:),pointer :: pay => null()
    real(kind=rkind),dimension(:,:,:),pointer :: pby => null()
    real(kind=rkind),dimension(:,:,:),pointer :: pdy => null()
    real(kind=rkind),dimension(:,:,:),pointer :: pax2 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: pbx2 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: pdx2 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: pay2 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: pby2 => null()
    real(kind=rkind),dimension(:,:,:),pointer :: pdy2 => null()

    ! media
    real(kind=rkind),dimension(:),pointer :: vp => null()
    real(kind=rkind),dimension(:),pointer :: vs => null()
    real(kind=rkind),dimension(:),pointer :: rho => null()
    real(kind=rkind),dimension(:),pointer :: zp => null()
    real(kind=rkind),dimension(:),pointer :: zs => null()
    real(kind=rkind),dimension(:),pointer :: lam => null()
    real(kind=rkind),dimension(:),pointer :: mu => null()


    ! fault recvs
    integer,dimension(:),pointer :: recv_fid => null()
    integer,dimension(:),pointer :: recv_i => null()
    integer,dimension(:),pointer :: recv_ie => null()
    real(kind=rkind),dimension(:),pointer :: recv_refx => null()
    real(kind=rkind),dimension(:,:,:),pointer :: recv_buffer => null()
    real(kind=rkind),dimension(:,:,:),pointer :: fault_buffer => null()

    real(kind=rkind),dimension(:,:,:),pointer :: surface_buffer => null()

    ! body recvs
    integer,dimension(:),pointer :: body_recv_fid => null()
    integer,dimension(:),pointer :: body_recv_i => null()
    integer,dimension(:),pointer :: body_recv_j => null()
    integer,dimension(:),pointer :: body_recv_ie => null()
    real(kind=rkind),dimension(:),pointer :: body_recv_refx => null()
    real(kind=rkind),dimension(:),pointer :: body_recv_refy => null()

end type

type buffvar
  real(kind=rkind), dimension(:,:), allocatable :: q_send,q_rec
  real(kind=rkind), dimension(:,:,:,:), allocatable :: qi
end type

type mytimer
  character(len=8) :: date
  character(len=10) :: time
  character(len=5) :: zone
  integer,dimension(8) :: values
end type

type wavevar
    real(kind=rkind),dimension(:,:,:),pointer :: u,hu,mu,tu
end type

contains

end module
