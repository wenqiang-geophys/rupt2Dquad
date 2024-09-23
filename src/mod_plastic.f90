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
  real(kind=RKIND) :: szz0, szz, sdzz

  real(kind=RKIND) :: angfric

  angfric = atan(blkfric)
  fluidpresh = 0
  syy0 = -120e6
  sxy0 = 70e6
  sxx0 = syy0-2*sxy0/tan(2*SmaxAngle/180.0*pi)
  szz0 = (sxx0+syy0)/2

  sxx = sxx0 + sxx
  syy = syy0 + syy
  sxy = sxy0 + sxy
  !szz = szz0 + 0

  !sm = (sxx+syy+szz)/3.0
  sm = (sxx+syy)/2.0
  sdxx = sxx - sm
  sdyy = syy - sm
  !sdzz = szz - sm
  sdxy = sxy

  !secinv = 0.5*(sdxx**2+sdyy**2+sdzz**2)+sdxy**2
  secinv = 0.5*(sdxx**2+sdyy**2)+sdxy**2
  tau = sqrt(secinv)
  taulim = cohesion*cos(angfric) - (sm+fluidpresh)*sin(angfric)
  taulim = max(0.0, taulim)

  if(tau .gt. taulim) then
    decay = exp(-dt/Tvisc)
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
    cp = mesh%vp(ie)
    cs = mesh%vs(ie)
    rho = mesh%rho(ie)

    do j = 1,Ngrid
    do i = 1,Ngrid

      depth = -mesh%vy(i,j,ie)
      exx = u(i+(j-1)*Ngrid,ie,3)
      eyy = u(i+(j-1)*Ngrid,ie,4)
      exy = u(i+(j-1)*Ngrid,ie,5)

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
