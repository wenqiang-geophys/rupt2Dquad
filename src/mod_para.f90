module mod_para

implicit none

#ifndef pOrder
#define pOrder 2
#endif

character(len=256) :: problem
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

! mpi
logical :: masternode

logical :: nice_print

real(kind=RKIND) :: simu_time_max
real(kind=RKIND) :: default_cfl_number
real(kind=RKIND) :: cfl_number
real(kind=RKIND) :: timestep

integer :: fault_snap_skip
integer :: wave_snap_skip
integer :: grdsurf_snap_skip

integer :: flux_method

integer :: use_damp
integer :: use_pml

! Friction laws
! 0 : linear slip weakening (default)
! 1 : rate state, ageing law
! 2 : rate state, slip law
! 3 : rate state, slip law, flash heating
! 4 : time weakening
integer :: friction_law

integer :: input_stress_type
! 0 : stress tensor Sxx0 Syy0 Sxy0
! 1 : local stress sigma0, tau0
! 2 : fault traction Tx0, Ty0

real(kind=rkind) :: RS_f0
real(kind=rkind) :: RS_V0
real(kind=rkind) :: RS_fw

! off-fault plasticity
integer :: plasticity
real(kind=rkind) :: cohesion
real(kind=rkind) :: blkfric
real(kind=rkind) :: Tvisc
real(kind=rkind) :: SmaxAngle
!real(kind=rkind) :: coef_byy
!real(kind=rkind) :: coef_bxx
!real(kind=rkind) :: coef_bxy
!real(kind=rkind) :: fluidpres_profile_h1
!real(kind=rkind) :: fluidpres_profile_h2
!real(kind=rkind) :: fluidpres_profile_o1
!real(kind=rkind) :: fluidpres_profile_o2

! thermal pressurization
integer :: thermalpressure

! io
integer :: export_grdsurf_velo
integer :: export_grdsurf_displ
integer :: export_grdsurf_strain
integer :: export_media
integer :: export_wave
integer :: export_wave_component
integer :: export_wave_timestep

! smoothly loading stress perturbation: T0+coef*dT0
integer :: smooth_load
real(kind=RKIND) :: smooth_load_time

! Initial condition
integer :: initial_condition_wave
real(kind=rkind) :: src_loc(3)
real(kind=rkind) :: src_gaussian_width
real(kind=rkind) :: src_mxx
real(kind=rkind) :: src_myy
real(kind=rkind) :: src_mzz
real(kind=rkind) :: src_myz
real(kind=rkind) :: src_mxz
real(kind=rkind) :: src_mxy
real(kind=rkind) :: src_m0

! parameters for time weakening law
real(kind=rkind) :: nucleate_y0
real(kind=rkind) :: nucleate_z0
real(kind=rkind) :: nucleate_rcrit
real(kind=rkind) :: nucleate_Vrup
real(kind=rkind) :: TimeForcedRup

contains


!subroutine para_init()
!end subroutine


end module
