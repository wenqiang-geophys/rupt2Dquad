module mod_plastic

  use mod_para
  use mod_types
  use mod_funcs

  implicit none

contains

subroutine Return_Map(sxx,syy,sxy,rho,depth,dt)
  implicit none
  real(kind=RKIND),intent(inout) :: sxx, syy, sxy
  real(kind=RKIND),intent(in) :: rho, depth, dt

  real(kind=RKIND) :: sxx0, syy0, sxy0, fluidpresh
  real(kind=RKIND) :: sm, sdxx, sdyy, sdxy
  real(kind=RKIND) :: secinv, tau, taulim, decay, yldfac
  real(kind=RKIND) :: omeg, h1, h2, omeg1, omeg2
  real(kind=RKIND) :: byy, bxx, bxy

  !real(kind=RKIND) :: cohes, blkfric, angfric, tvisc!, dist
  real(kind=RKIND) :: angfric
  !integer :: i, j, k

  !real(kind=RKIND) :: lam,miu,chi

  ! For tpv27
  cohesion = 1.36e6
  blkfric = 0.1934
  tvisc = 0.03

  !!! For tpv29
  !cohes = 1.18e6
  !blkfric = 0.1680
  !tvisc = 0.05

  angfric = atan(blkfric)
  !print*, 'blkfric = ', blkfric, 'angfric = ',angfric

  !do k = nk1, nk2
  !  do j = nj1, nj2
  !    do i = ni1, ni2
        !sxx = Txx(i,j,k); syy = Tyy(i,j,k); szz = Tzz(i,j,k)
        !sxy = Txy(i,j,k); sxz = Txz(i,j,k); syz = Tyz(i,j,k)
        !depth = -z(i,j,k)
        !depth = -z

  ! For tpv27
  !@fluidpresh = 1.0*9.8*depth ! rho*g*h
  !@if(depth<=15.0)then
  !@  omeg = 1.0
  !@elseif(depth<=20.0)then
  !@  omeg = (20.0 - depth)/5.0
  !@else
  !@  omeg = 0.0
  !@endif
  !@!szz0 = min(-2670.0*9.8*steph/3.0, -2670.0*9.8*depth)
  !@!szz0 = -2670.0*9.8*depth
  !@szz0 = -2.670*9.8*depth ! in MPa
  !@syy0 = omeg*(0.926793*(szz0+fluidpresh)-fluidpresh) + (1.0-omeg)*szz0
  !@sxx0 = omeg*(1.073206*(szz0+fluidpresh)-fluidpresh) + (1.0-omeg)*szz0
  !@sxy0 = omeg*(0.169029*(szz0+fluidpresh))
  !@sxz0 = 0.0
  !@syz0 = 0.0

  h1 = 17e3
  h2 = 22e3
  omeg1 = 1.0
  omeg2 = 0.0
  byy = 1.025837
  bxx = 0.974162
  bxy = 0.158649
  !h1 = fluidpres_profile_h1
  !h2 = fluidpres_profile_h2
  !omeg1 = fluidpres_profile_o1
  !omeg2 = fluidpres_profile_o2
  !byy = coef_byy
  !bxx = coef_bxx
  !bxy = coef_bxy

  ! For tpv30
  fluidpresh = 1.0*9.8*depth
  if(depth<=h1)then
    omeg = omeg1
  elseif(depth<=h2)then
    omeg = (h2-depth)/(h2-h1)
  else
    omeg = omeg2
  endif
  syy0 = min(-1.0e-3, -rho*9.8*depth) ! Pa
  sxx0 = omeg*(bxx*(syy0+fluidpresh)-fluidpresh) + (1.0-omeg)*syy0
  sxy0 = omeg*(bxy*(syy0+fluidpresh))

  fluidpresh = 0
  syy0 = -120e6
  sxy0 = 70e6
  sxx0 = 0

  sxx = sxx0 + sxx
  syy = syy0 + syy
  sxy = sxy0 + sxy

  sm = (sxx+syy)/2.0
  sdxx = sxx - sm
  sdyy = syy - sm
  sdxy = sxy

  secinv = 0.5*(sdxx**2+sdyy**2)+sdxy**2
  tau = sqrt(secinv)
  taulim = cohesion*cos(angfric) - (sm+fluidpresh)*sin(angfric)
  taulim = max(0.0, taulim)

  if(tau .gt. taulim) then
    decay = exp(-dt/tvisc)
    yldfac = decay + (1.0-decay)*taulim/tau
    sxx = sdxx*yldfac + sm
    syy = sdyy*yldfac + sm
    sxy = sdxy*yldfac
  endif

  sxx = sxx - sxx0
  syy = syy - syy0
  sxy = sxy - sxy0
endsubroutine

subroutine update_plastic(mesh,u)
  implicit none
  type(meshvar) :: mesh
  real(kind=RKIND),dimension(:,:,:) :: u ! (Np*Nelem,Nvar)
  real(kind=RKIND) :: dt,depth,rho,cp,cs
  real(kind=RKIND) :: exx,eyy,exy
  real(kind=RKIND) :: sxx,syy,sxy
  integer :: ie,i,j

  dt = mesh%deltat

  do ie = 1,mesh%Nelem
    !do j = 1,Np
    !  iv(j) = j+(ie-1)*Np
    !end do
    cp = mesh%vp(ie)
    cs = mesh%vs(ie)
    rho = mesh%rho(ie)

    do i = 1,Ngrid
    do j = 1,Ngrid

      depth = -mesh%vy(i,j,ie)
      exx = u(i+(j-1),ie,3)
      eyy = u(i+(j-1),ie,4)
      exy = u(i+(j-1),ie,5)

      call strain2stress(exx,eyy,exy,rho,cp,cs,sxx,syy,sxy)

      call Return_Map(sxx,syy,sxy,rho,depth,dt)

      call stress2strain(sxx,syy,sxy,rho,cp,cs,exx,eyy,exy)

      u(i+(j-1)*Ngrid,ie,3) = exx
      u(i+(j-1)*Ngrid,ie,4) = eyy
      u(i+(j-1)*Ngrid,ie,5) = exy
    end do
    end do
  end do ! element

end subroutine

end module
