module mod_para

implicit none

character(len=256) :: mesh_dir
character(len=256) :: data_dir

integer,parameter :: rkind = 8
integer,parameter :: SIZE_REAL = 4
integer,parameter :: SIZE_DOUBLE = 8
integer,parameter :: CUSTOM_REAL = SIZE_DOUBLE

integer,parameter :: nsurface = 4
integer,parameter :: Nfaces = 4
integer,parameter :: FDorder = 8
integer,parameter :: Order = pOrder
integer,parameter :: Nsub = 1
integer,parameter :: Ngrid = Nsub*Order+1
integer,parameter :: Order2 = order
integer,parameter :: NGLL = Order+1
integer,parameter :: NGLL2 = Order2+1
!integer,parameter :: Np = (Order+1)*(Order+2)/2
!integer,parameter :: Nfp = NGLL
integer,parameter :: Nfp = Ngrid
!integer,parameter :: Np = NGLL**2
integer,parameter :: Np = Ngrid**2
#ifdef PSV
integer,parameter :: Nvar = 5 !P-SV:5, SH:3
#endif
#ifdef SH
integer,parameter :: Nvar = 3 !P-SV:5, SH:3
#endif
integer,parameter :: dimens = 8

!real(kind=rkind),parameter :: steph = 250.0d0
real(kind=rkind),parameter :: steph = 1.0d0

real(kind=rkind),parameter :: PI = 3.141592653589793238463
real(kind=rkind),parameter :: EPS = 1e-5
real(kind=rkind),parameter :: GRAVITY = 9.8

!real(kind=rkind),parameter :: damp_s = 0.3

! Boundary type
integer,parameter :: BC_IN             = 0
integer,parameter :: BC_SOLID_SOLID    = 0
integer,parameter :: BC_FLUID_FLUID    = 1
integer,parameter :: BC_SOLID_FLUID    = 2
integer,parameter :: BC_FREE     = -1
integer,parameter :: BC_OUT      = -2
integer,parameter :: BC_FREE_G   = -4
integer,parameter :: BC_FAULT    = 100

! Element type
integer,parameter :: ELEM_SOLID = 0
integer,parameter :: ELEM_FLUID = 1
!integer :: Order, NGLL, nsurface
!integer :: Np

! Low storage Runge-Kutta coefficients
real(kind=rkind),parameter :: rk4a(5) = (/ &
        0.0,                               &
        -567301805773.0/1357537059087.0,   &
        -2404267990393.0/2016746695238.0,  &
        -3550918686646.0/2091501179385.0,  &
        -1275806237668.0/842570457699.0/)
real(kind=rkind),parameter :: rk4b(5) = (/ &
        1432997174477.0/9575080441755.0,   &
        5161836677717.0/13612068292357.0,  &
        1720146321549.0/2090206949498.0,   &
        3134564353537.0/4481467310338.0,   &
        2277821191437.0/14882151754819.0/)
real(kind=rkind),parameter :: rk4c(5) = (/ &
        0.0,                               &
        1432997174477.0/9575080441755.0,   &
        2526269341429.0/6820363962896.0,   &
        2006345519317.0/3224310063776.0,   &
        2802321613138.0/2924317926251.0/)

contains


!subroutine para_init()
!end subroutine


end module
