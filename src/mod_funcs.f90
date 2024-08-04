module mod_funcs

use mod_para
use mod_gll
use mod_fd
use mod_mesh

contains

subroutine Flux1(U,dU,rho,cp,cs,Fx)
    implicit none
    real(kind=rkind),intent(in) :: U(:),dU(:),rho,cp,cs
    real(kind=rkind),intent(out) :: Fx(:)
    real(kind=rkind) :: rrho,miu,lam,chi

    ! U=(rho*Vx,rho*Vy,Exx,Eyy,Exy)
    ! f=(Sxx,Sxy,Vx,0,Vy)
    ! g=(Sxy,Syy,0,Vy,Vx)
    rrho = 1d0/rho
    miu = rho*cs**2
    chi = rho*cp**2
    lam = chi-2d0*miu
    Fx(:) = 0

#ifdef STRAIN
    Fx(1) = U(3) * chi + U(4) * lam
    Fx(2) = U(5) * miu
    Fx(3) = U(1) * rrho
    Fx(4) = 0
    Fx(5) = U(2) * rrho
    Fx(1) = (U(3)+U(3+3)) * chi + (U(4)+U(4+3)) * lam
    Fx(2) = (U(5)+U(5+3)) * miu
#else
    ! U = (Vx,Vy,Sxx,Syy,Sxy)
    Fx(1) = U(3) * rrho
    Fx(2) = U(5) * rrho
    Fx(3) = U(1) * chi
    Fx(4) = U(1) * lam
    Fx(5) = U(2) * miu
#endif
#ifdef SYM
    ! U = (Vx,Vy,Sxx,Syy,Sxy)
    Fx(1) = U(3)
    Fx(2) = U(5)
    Fx(3) = U(1)
    Fx(4) = 0
    Fx(5) = U(2)
    Fx(6:8) = 0
#endif
end subroutine

subroutine Flux2(U,dU,rho,cp,cs,Fy)
    implicit none
    real(kind=rkind),intent(in) :: U(:),dU(:),rho,cp,cs
    real(kind=rkind),intent(out) :: Fy(:)
    real(kind=rkind) :: rrho,miu,lam,chi
    ! U=(rho*Vx,rho*Vy,Exx,Eyy,Exy)
    ! f=(Sxx,Sxy,Vx,0,Vy)
    ! g=(Sxy,Syy,0,Vy,Vx)
    rrho = 1d0/rho
    miu = rho*cs**2
    chi = rho*cp**2
    lam = chi-2d0*miu
    Fy(:) = 0

#ifdef STRAIN
    Fy(1) = U(5) * miu
    Fy(2) = U(3) * lam + U(4) * chi
    Fy(3) = 0
    Fy(4) = U(2) * rrho
    Fy(5) = U(1) * rrho
    Fy(1) = (U(5)+U(5+3)) * miu
    Fy(2) = (U(3)+U(3+3)) * lam + (U(4)+U(4+3)) * chi
#else
    ! U = (Vx,Vy,Sxx,Syy,Sxy)
    Fy(1) = U(5) * rrho
    Fy(2) = U(4) * rrho
    Fy(3) = U(2) * lam
    Fy(4) = U(2) * chi
    Fy(5) = U(1) * miu
#endif
#ifdef SYM
    ! U = (Vx,Vy,Sxx,Syy,Sxy)
    Fy(1) = U(5)
    Fy(2) = U(4)
    Fy(3) = 0
    Fy(4) = U(2)
    Fy(5) = U(1)
    Fy(6:8) = 0
#endif
end subroutine

subroutine Flux1_acoustic(U,rho,cp,Fx)
    implicit none
    real(kind=rkind),intent(in) :: U(:),rho,cp
    real(kind=rkind),intent(out) :: Fx(:)
    real(kind=rkind) :: chi
    ! U=(rho*Vx,rho*Vy,p/kappa)
    ! f=-(p,0,Vx)
    ! g=-(0,p,Vy)
    chi = rho*cp**2

    Fx(1) = -U(3) * chi
    Fx(2) = -0
    Fx(3) = -U(1) / rho
    Fx(4) = 0
    Fx(5) = 0
end subroutine

subroutine Flux2_acoustic(U,rho,cp,Fy)
    implicit none
    real(kind=rkind),intent(in) :: U(:),rho,cp
    real(kind=rkind),intent(out) :: Fy(:)
    real(kind=rkind) :: chi
    ! U=(rho*Vx,rho*Vy,p/kappa)
    ! f=-(p,0,Vx)
    ! g=-(0,p,Vy)
    chi = rho*cp**2

    Fy(1) = -0
    Fy(2) = -U(3) * chi
    Fy(3) = -U(2) / rho
    Fy(4) = 0
    Fy(5) = 0
end subroutine

!old subroutine get_Ar(rho,cp,cs,rx,ry,detJ,Ar)
!old     implicit none
!old     real(kind=rkind),intent(in) :: rho,cp,cs,rx,ry,detJ
!old     real(kind=rkind),intent(out) :: Ar(:,:)
!old     real(kind=rkind) :: rrho,miu,lam,chi
!old 
!old     real(kind=rkind) :: Ax(5,5),Ay(5,5)
!old 
!old     rrho = 1d0/rho
!old     miu = rho*cs**2
!old     chi = rho*cp**2
!old     lam = chi-2d0*miu
!old 
!old     Ax(:,:) = 0d0
!old     Ay(:,:) = 0d0
!old 
!old     Ax(1,3) = rrho
!old     Ax(2,5) = rrho
!old     Ax(3,1) = chi
!old     Ax(4,1) = lam
!old     Ax(5,2) = miu
!old 
!old     Ay(1,5) = rrho
!old     Ay(2,4) = rrho
!old     Ay(3,2) = lam
!old     Ay(4,2) = chi
!old     Ay(5,1) = miu
!old 
!old     Ar = detJ * (rx*Ax+ry*Ay)
!old end subroutine


!old subroutine RefFlux(U,rho,cp,cs,rx,ry,sx,sy,detJ,modFx,modFy)
!old     implicit none
!old     real(kind=rkind),intent(in) :: U(:),rho,cp,cs,rx,ry,sx,sy,detJ
!old     real(kind=rkind),intent(out) :: modFx(:),modFy(:)
!old     real(kind=rkind) :: Fx(5),Fy(5)
!old 
!old     call Flux1(U,rho,cp,cs,Fx)
!old     call Flux2(U,rho,cp,cs,Fy)
!old     modFx = detJ * (rx*Fx+ry*Fy)
!old     modFy = detJ * (sx*Fx+sy*Fy)
!old end subroutine

!subroutine numflux1(fL,fR,rho,cp,cs,fstar)
subroutine numflux1(fL,fR,zp_in,zp_out,zs_in,zs_out,fstar)
    implicit none
    real(kind=rkind),intent(in) :: fL(:),fR(:),zp_in,zs_in,zp_out,zs_out
    real(kind=rkind),intent(out) :: fstar(:)
    real(kind=rkind) :: a1,a2
    !Zp = rho*cp
    !Zs = rho*cs
    !a1 = 0.5d0*(fR(3)-fL(3)) + 0.5d0*(fR(1)-fL(1))/Zp
    !a2 = 0.5d0*(fR(5)-fL(5)) + 0.5d0*(fR(2)-fL(2))/Zs
    ! f (sxx,sxy,vx,0,vy)
    a1 = ( (fR(1)-fL(1)) + (fR(3)-fL(3))*zp_out) / (zp_in+zp_out)
    a2 = ( (fR(2)-fL(2)) + (fR(5)-fL(5))*zs_out) / (zs_in+zs_out)

    fstar(1) = 0.5d0*(fR(1)-fL(1)) + 0.5d0*(fR(3)-fL(3))*Zp_in
    fstar(2) = 0.5d0*(fR(2)-fL(2)) + 0.5d0*(fR(5)-fL(5))*Zs_in
    fstar(3) = 0.5d0*(fR(3)-fL(3)) + 0.5d0*(fR(1)-fL(1))/Zp_in
    fstar(4) = 0.0d0
    fstar(5) = 0.5d0*(fR(5)-fL(5)) + 0.5d0*(fR(2)-fL(2))/Zs_in
    !fstar(1) = a1*zp_in
    !fstar(2) = a2*zs_in
    !fstar(3) = a1
    !fstar(4) = 0.0d0
    !fstar(5) = a2
end subroutine

subroutine numflux1_mod(fL,fR,zp_in,zp_out,zs_in,zs_out,gam,fstar)
    implicit none
    real(kind=rkind),intent(in) :: fL(:),fR(:),zp_in,zs_in,zp_out,zs_out,gam
    real(kind=rkind),intent(out) :: fstar(:)
    real(kind=rkind) :: a1,a2
    !Zp = rho*cp
    !Zs = rho*cs
    !a1 = 0.5d0*(fR(3)-fL(3)) + 0.5d0*(fR(1)-fL(1))/Zp
    !a2 = 0.5d0*(fR(5)-fL(5)) + 0.5d0*(fR(2)-fL(2))/Zs
    ! f (sxx,sxy,vx,0,vy)
    a1 = ( (fR(1)-fL(1)) + (fR(3)-fL(3))*zp_out) / (zp_in+zp_out)
    a2 = ( (fR(2)-fL(2)) + (fR(5)-fL(5))*zs_out) / (zs_in+zs_out)

    fstar(1) = 0.5d0*(fR(1)-fL(1)) + gam*0.5d0*(fR(3)-fL(3))*Zp_in
    fstar(2) = 0.5d0*(fR(2)-fL(2)) + gam*0.5d0*(fR(5)-fL(5))*Zs_in
    fstar(3) = 0.5d0*(fR(3)-fL(3)) + gam*0.5d0*(fR(1)-fL(1))/Zp_in
    fstar(4) = 0.0d0
    fstar(5) = 0.5d0*(fR(5)-fL(5)) + gam*0.5d0*(fR(2)-fL(2))/Zs_in
    !fstar(1) = a1*zp_in
    !fstar(2) = a2*zs_in
    !fstar(3) = a1
    !fstar(4) = 0.0d0
    !fstar(5) = a2
end subroutine

subroutine numflux1_elastic_acoustic(fL,fR,zp_in,zp_out,zs_in,zs_out,fstar)
    implicit none
    real(kind=rkind),intent(in) :: fL(:),fR(:),zp_in,zp_out,zs_in,zs_out
    real(kind=rkind),intent(out) :: fstar(:)
    real(kind=rkind) :: a1,a2
    !Zp = rho*cp
    !Zs = rho*cs
    !a1 = 0.5d0*(-fR(3)-fL(3)) + 0.5d0*(fR(1)-fL(1))/Zp
    !a1 = 0.5d0 * ( (-fR(3)) - fL(3) ) - 0.5d0 * ( (-fR(1) ) + fL(1) )/Zp
    !a2 = -fL(2)/Zs

    a1 = ( fR(1)-fL(1) + (-fR(3)-fL(3))*zp_out ) / (zp_in+zp_out)
    a2 = -fL(2)/zs_in

    ! fL (sxx,sxy,vx,0,vy)
    ! fR -(p,0,vx)
    fstar(1) = a1*zp_in
    fstar(2) = a2*zs_in
    fstar(3) = a1
    fstar(4) = 0.0d0
    fstar(5) = a2
end subroutine

subroutine numflux1_acoustic(fL,fR,zp_in,zp_out,fstar)
    implicit none
    real(kind=rkind),intent(in) :: fL(:),fR(:),zp_in,zp_out
    real(kind=rkind),intent(out) :: fstar(:)
    real(kind=rkind) :: a1
    !Zp = rho*cp
    !a1 = 0.5d0 * ( (fR(3)-fL(3)) - (fR(1)-fL(1))/Zp )
    a1 = - ( (fR(1)-fL(1)) - (fR(3)-fL(3))*zp_out ) / (zp_in+zp_out)
    ! f -(p,0,vx)
    fstar(1) = -a1*zp_in
    fstar(2) = 0
    fstar(3) = a1
    fstar(4) = 0
    fstar(5) = 0
end subroutine

subroutine numflux1_acoustic_elastic(fL,fR,zp_in,zp_out,fstar)
    implicit none
    real(kind=rkind),intent(in) :: fL(:),fR(:),zp_in,zp_out
    real(kind=rkind),intent(out) :: fstar(:)
    real(kind=rkind) :: a1
    !Zp = rho*cp
    !a1 = -0.5d0*( (fR(3)+fL(3)) + (fR(1)-fL(1))/Zp )
    !a1 = -0.5d0*( ( fR(3) - ( -fL(3) ) ) + ( fR(1) + (-fL(1)) )/Zp )
    a1 = -( ( fR(1) + (-fL(1)) ) + ( fR(3) - ( -fL(3) ) )*zp_out ) / (zp_in+zp_out)
    ! fL -(p,0,vx)
    ! fR (sxx,sxy,vx,0,vy)
    fstar(1) = -a1*zp_in
    fstar(2) = 0
    fstar(3) = a1
    fstar(4) = 0
    fstar(5) = 0
end subroutine

function getInvT(nxe,nze)
    real(kind=rkind), intent(in) :: nxe,nze
    real(kind=rkind), dimension(:,:), allocatable :: getInvT
    real(kind=rkind) :: matD(3,3),invD(3,3)
    allocate(getInvT(5,5))
    getInvT=0.0
    getInvT(1,1) = nxe
    getInvT(1,2) = nze
    getInvT(2,1) = -nze
    getInvT(2,2) = nxe
#ifdef STRAIN
    getInvT(3,3) = nxe*nxe
    getInvT(3,4) = nze*nze
    getInvT(3,5) = nxe*nze
    getInvT(4,3) = nze*nze
    getInvT(4,4) = nxe*nxe
    getInvT(4,5) = -nxe*nze
    getInvT(5,3) = -2*nxe*nze
    getInvT(5,4) = 2*nxe*nze
    getInvT(5,5) = nxe*nxe-nze*nze
#else
    getInvT(3,3) = nxe*nxe
    getInvT(3,4) = nze*nze
    getInvT(3,5) = 2*nxe*nze
    getInvT(4,3) = nze*nze
    getInvT(4,4) = nxe*nxe
    getInvT(4,5) = -2*nxe*nze
    getInvT(5,3) = -nxe*nze
    getInvT(5,4) = nxe*nze
    getInvT(5,5) = nxe*nxe-nze*nze
#endif
!    matD = 0.0d0
!    matD(1,1) = 1.0d0
!    matD(2,2) = 1.0d0
!    matD(3,3) = 2.0d0
!    invD = 0.0d0
!    invD(1,1) = 1.0d0
!    invD(2,2) = 1.0d0
!    invD(3,3) = 0.5d0
!#ifdef STRAIN
!    !!! matmul for small matrices is quite inefficient !!!
!    getInvT(3:5,3:5) = matmul(matmul(matD,getInvT(3:5,3:5)),invD)
!#endif
end function

function getT(nxe,nze)
    real(kind=rkind), intent(in) :: nxe,nze
    real(kind=rkind), dimension(:,:), allocatable :: getT
    real(kind=rkind) :: matD(3,3),invD(3,3)
    allocate(getT(5,5))
    getT=0.0
    getT( 1, 1) = nxe
    getT( 1, 2) = -nze
    getT( 2, 1) = nze
    getT( 2, 2) = nxe
#ifdef STRAIN
    getT( 3, 3) = nxe*nxe
    getT( 3, 4) = nze*nze
    getT( 3, 5) = -nxe*nze
    getT( 4, 3) = nze*nze
    getT( 4, 4) = nxe*nxe
    getT( 4, 5) = nxe*nze
    getT( 5, 3) = 2*nxe*nze
    getT( 5, 4) = -2*nxe*nze
    getT( 5, 5) = nxe*nxe-nze*nze
#else
    getT( 3, 3) = nxe*nxe
    getT( 3, 4) = nze*nze
    getT( 3, 5) = -2*nxe*nze
    getT( 4, 3) = nze*nze
    getT( 4, 4) = nxe*nxe
    getT( 4, 5) = 2*nxe*nze
    getT( 5, 3) = nxe*nze
    getT( 5, 4) = -nxe*nze
    getT( 5, 5) = nxe*nxe-nze*nze
#endif
!    matD = 0.0d0
!    matD(1,1) = 1.0d0
!    matD(2,2) = 1.0d0
!    matD(3,3) = 2.0d0
!    invD = 0.0d0
!    invD(1,1) = 1.0d0
!    invD(2,2) = 1.0d0
!    invD(3,3) = 0.5d0
!#ifdef STRAIN
!    !!! matmul for small matrices is quite inefficient !!!
!    getT(3:5,3:5) = matmul(matmul(matD,getT(3:5,3:5)),invD)
!#endif
end function

subroutine rotate_u(Tv,Ts,u)
    implicit none
    real*8,intent(inout) :: u(8)
    real*8,intent(in) :: Ts(3,3),Tv(2,2)
    real*8 :: u1(8)

    u1(:) = 0
    !u(1:2) = matmul(Tv,u(1:2))
    u1(1) = Tv(1,1)*u(1)+Tv(1,2)*u(2)
    u1(2) = Tv(2,1)*u(1)+Tv(2,2)*u(2)
    !u(3:5) = matmul(Ts,u(3:5))
    u1(3) = Ts(1,1)*u(3)+Ts(1,2)*u(4)+Ts(1,3)*u(5)
    u1(4) = Ts(2,1)*u(3)+Ts(2,2)*u(4)+Ts(2,3)*u(5)
    u1(5) = Ts(3,1)*u(3)+Ts(3,2)*u(4)+Ts(3,3)*u(5)

    u = u1
end subroutine

subroutine rotate_u_acoustic(Tv,u)
    implicit none
    real*8,intent(inout) :: u(8)
    real*8,intent(in) :: Tv(2,2)
    real*8 :: u1(8)

    u1(:) = 0
    !u(1:2) = matmul(Tv,u(1:2))
    u1(1) = Tv(1,1)*u(1)+Tv(1,2)*u(2)
    u1(2) = Tv(2,1)*u(1)+Tv(2,2)*u(2)
    u1(3) = u(3)
    u1(4) = 0
    u1(5) = 0

    u = u1
end subroutine

subroutine rotate_xy2nm(n,m,Tx,Ty,Tn,Tm)
    implicit none
    real(kind=rkind),intent(in) :: n(2),m(2),Tx,Ty
    real(kind=rkind),intent(out) :: Tn,Tm

    Tn = Tx*n(1)+Ty*n(2)
    Tm = Tx*m(1)+Ty*m(2)
end subroutine

subroutine rotate_nm2xy(n,m,Tn,Tm,Tx,Ty)
    implicit none
    real(kind=rkind),intent(in) :: n(2),m(2),Tn,Tm
    real(kind=rkind),intent(out) :: Tx,Ty

    Tx = Tn*n(1)+Tm*m(1)
    Ty = Tn*n(2)+Tm*m(2)
end subroutine

subroutine extract_traction_velocity(u,n,vx,vy,Tx,Ty)
    implicit none
    !real(kind=rkind),intent(in) :: u(5),n(2)
    real(kind=rkind),intent(in) :: u(5),n(2)
    real(kind=rkind),intent(out) :: vx,vy,Tx,Ty
    real(kind=rkind) :: sxx,syy,sxy

    vx = u(1)
    vy = u(2)
    sxx = u(3)
    syy = u(4)
    sxy = u(5)

    Tx = sxx*n(1) + sxy*n(2)
    Ty = sxy*n(1) + syy*n(2)
end subroutine

!subroutine extract_traction_velocity_strain(u,n,rho,lam,miu,vx,vy,Tx,Ty)
subroutine extract_traction_velocity_strain(u,n,rho,cp,cs,vx,vy,Tx,Ty)
    implicit none
    real(kind=rkind),intent(in) :: u(5),n(2),rho,cp,cs
    real(kind=rkind),intent(out) :: vx,vy,Tx,Ty
    real(kind=rkind) :: sxx,syy,sxy
    real(kind=rkind) :: exx,eyy,exy
    real(kind=rkind) :: chi, lam, miu

    miu = rho*cs**2
    chi = rho*cp**2
    lam = chi-2d0*miu
    !chi = lam+2.0*miu

    vx = u(1)/rho
    vy = u(2)/rho
    exx = u(3)
    eyy = u(4)
    exy = u(5)

    ! sxx sxy syy
    !Fx(1) = U(3) * chi + U(4) * lam
    !Fx(2) = U(5) * miu
    !Fy(1) = U(5) * miu
    !Fy(2) = U(3) * lam + U(4) * chi
    sxx = exx * chi + eyy * lam
    syy = exx * lam + eyy * chi
    sxy = exy * miu

    Tx = sxx*n(1) + sxy*n(2)
    Ty = sxy*n(1) + syy*n(2)
end subroutine

subroutine riemannSolver_continuous(v_p,v_m,sigma_p,sigma_m,z_p,z_m,v_hat_p,v_hat_m,sigma_hat_p,sigma_hat_m)
    implicit none
    real(kind=rkind),intent(in) :: v_p,v_m,sigma_p,sigma_m,z_p,z_m
    real(kind=rkind),intent(out) :: v_hat_p,v_hat_m,sigma_hat_p,sigma_hat_m
    real(kind=rkind) :: p
    real(kind=rkind) :: q
    real(kind=rkind) :: phi
    !real(kind=rkind) :: v_hat
    real(kind=rkind) :: eta

    p=z_m*v_p + sigma_p
    q=z_p*v_m - sigma_m
    ! fix bug
    p=z_p*v_p + sigma_p
    q=z_m*v_m - sigma_m

    !eta=(z_p*z_m)/(z_p+z_m)
    !phi= eta*(p/z_p - q/z_m)
    !sigma_hat_p=phi
    !sigma_hat_m=phi
    !v_hat_p=(q+phi)/z_m
    !v_hat_m=(p-phi)/z_p

    if(z_p > 0 .and. z_m > 0) then
      eta=(z_p*z_m)/(z_p+z_m);
      phi= eta*(p/z_p - q/z_m);
      sigma_hat_p=phi;
      sigma_hat_m=phi;
      v_hat_p=(q+phi)/z_m;
      v_hat_m=(p-phi)/z_p;
    else if(z_p > 0) then
      sigma_hat_p=0;
      sigma_hat_m=sigma_m;
      v_hat_p=v_p;
      v_hat_m=v_m;
    else if(z_m > 0) then
      sigma_hat_p=sigma_p;
      sigma_hat_m=0;
      v_hat_p=v_p;
      v_hat_m=v_m;
    else
      sigma_hat_p=sigma_p;
      sigma_hat_m=sigma_m;
      v_hat_p=v_p;
      v_hat_m=v_m;
    end if

end subroutine

subroutine generate_fluctuations_left(z,T,T_hat,v,v_hat,F)
    implicit none
    real(kind=rkind),intent(in) :: z,T,T_hat,v,v_hat
    real(kind=rkind),intent(out) :: F
    F = -0.5d0*(z*(v-v_hat) + (T-T_hat))
end subroutine

subroutine generate_fluctuations_right(z,T,T_hat,v,v_hat,F)
    implicit none
    real(kind=rkind),intent(in) :: z,T,T_hat,v,v_hat
    real(kind=rkind),intent(out) :: F
    F = -0.5d0*(z*(v-v_hat) - (T-T_hat))
end subroutine

subroutine SlipWeakeningFriction(vn_p, vn_m, Tn_p, Tn_m, zn_p, zn_m, &
                                 vn_hat_p, vn_hat_m, Tn_hat_p, Tn_hat_m, &
                                 vm_p,vm_m, Tm_p, Tm_m, zm_p , zm_m, &
                                 vm_hat_p, vm_hat_m, Tm_hat_p, Tm_hat_m, &
                                 T0_n, T0_m, S, Vel, Dc, mu_s, mu_d)
    implicit none
    real(kind=rkind),intent(in) :: vn_p, vn_m, Tn_p, Tn_m, zn_p, zn_m
    real(kind=rkind),intent(in) :: vm_p,vm_m, Tm_p, Tm_m, zm_p , zm_m
    real(kind=rkind),intent(in) :: T0_n, T0_m, S, Dc, mu_s, mu_d
    real(kind=rkind),intent(out) :: vn_hat_p,vn_hat_m,Tn_hat_p,Tn_hat_m
    real(kind=rkind),intent(out) :: vm_hat_p,vm_hat_m,Tm_hat_p,Tm_hat_m
    real(kind=rkind),intent(out) :: Vel

    real(kind=rkind) :: p_m,p_n,q_m,q_n
    real(kind=rkind) :: eta_s,eta_n
    real(kind=rkind) :: phi_m,phi_n
    real(kind=rkind) :: Tau_lock,sigma_n,tau_str
    real(kind=rkind) :: Tl, Tm, vv_l, vv_m
    ! compute characteristics
    p_m = zm_p*vm_p + Tm_p
    p_n = zn_p*vn_p + Tn_p

    q_m = zm_m*vm_m - Tm_m
    q_n = zn_m*vn_m - Tn_m

    ! half of the harmonic mean of Z1_s, Z2_s
    eta_s=(zm_p*zm_m)/(zm_p+zm_m)
    eta_n=(zn_p*zn_m)/(zn_p+zn_m)

    phi_m=eta_s*(p_m/zm_p - q_m/zm_m)
    phi_n=eta_n*(p_n/zn_p - q_n/zn_m)

    !double T0_m=0
    !double T0_n=0

    ! get prestress (where normal traction is effective normal traction)
    !prestress(T0_n, T0_m, x, 0.0, m, n)

    vn_hat_m = (p_n - phi_n)/zn_p   ! continuity of normal velocity
    Tn_hat_m = phi_n                ! continuity of normal stress
    vn_hat_p = (q_n + phi_n)/zn_m   ! continuity of normal velocity
    Tn_hat_p = phi_n

    tau_lock = dsqrt((T0_m + phi_m)**2)
    sigma_n = max(0.0, -(T0_n + phi_n))   ! including prestress

    call TauStrength(tau_str, sigma_n, S, Dc, mu_s, mu_d)

    !Vel = 0
    if (tau_lock >= tau_str) then
        call slip_weakening(vv_m, Vel, Tm, phi_m+T0_m, eta_s, tau_str, sigma_n)
        Tm_hat_m = Tm - T0_m
        Tm_hat_p = Tm - T0_m
    else
        Tm_hat_m = phi_m
        Tm_hat_p = phi_m
        vv_m = 0.0
        Vel = 0.0
    end if

    vm_hat_p = (Tm_hat_m + q_m)/zm_m + vv_m
    vm_hat_m = (p_m - Tm_hat_p)/zm_p - vv_m
end subroutine

!! solve for slip-rate (vv):
subroutine slip_weakening(v1,Vel,tau1,phi_1,eta,tau_str,sigma_n)
    implicit none
    real(kind=rkind),intent(in) :: phi_1,eta,tau_str,sigma_n
    real(kind=rkind),intent(out) :: v1,Vel,tau1
    real(kind=rkind) :: Phi

    Phi = abs(phi_1)            ! stress-transfer functional
    Vel = (Phi - tau_str)/eta   ! slip-rate

    !compute slip velocities
    v1 = phi_1/(eta+tau_str/Vel)
    !compute shear stress on the fault
    tau1 = phi_1 - eta*v1
end subroutine

subroutine TauStrength(tau_str, sigma_n, S, S_c, mu_s, mu_d)
    implicit none
    real(kind=rkind),intent(in) :: sigma_n,S,S_c,mu_s,mu_d
    real(kind=rkind),intent(out) :: tau_str
    real(kind=rkind) :: fric_coeff

    !mu_s = 0.677                        ! stastic friction
    !mu_d = 0.525                        ! dynamic friction
    !sigma0 = 120.0e6                    ! normal stress
    !S_c = 0.4                           ! critical slip

    fric_coeff = mu_s - (mu_s-mu_d) * min(abs(S),S_c)/S_c
    tau_str = fric_coeff*sigma_n   
end subroutine

subroutine get_face1(u,face,u_face)
    implicit none
    ! u(Ngrid,Ngrid)
    real(kind=rkind),intent(in) :: u(:,:)
    integer,intent(in) :: face
    real(kind=rkind),intent(out) :: u_face(:)

    if (face == 1) then
        u_face = u(:,1)
    elseif (face == 2) then
        u_face = u(Ngrid,:)
    elseif (face == 3) then
        u_face = u(:,Ngrid)
    elseif (face == 4) then
        u_face = u(1,:)
    end if
end subroutine

!old subroutine get_face(u,i,face,ie,u_face)
!old     implicit none
!old     ! u(NGLL,NGLL,nelem,5)
!old     real(kind=rkind),intent(in) :: u(:,:,:,:)!,xi(:)
!old     integer,intent(in) :: i,ie,face
!old     real(kind=rkind),intent(out) :: u_face(:)
!old     !integer :: ivar
!old
!old     if (face == 1) then
!old         u_face = u(i,1,ie,:)
!old         !do ivar = 1,5
!old         !    u_face(ivar) = lagrange_interpol(order+1,0d0,xi,u(i,:,ie,ivar))
!old         !end do
!old     elseif (face == 2) then
!old         u_face = u(NGLL,i,ie,:)
!old         !do ivar = 1,5
!old         !    u_face(ivar) = lagrange_interpol(order+1,1d0,xi,u(:,i,ie,ivar))
!old         !end do
!old     elseif (face == 3) then
!old         u_face = u(i,NGLL,ie,:)
!old         !do ivar = 1,5
!old         !    u_face(ivar) = lagrange_interpol(order+1,1d0,xi,u(i,:,ie,ivar))
!old         !end do
!old     elseif (face == 4) then
!old         u_face = u(1,i,ie,:)
!old         !do ivar = 1,5
!old         !    u_face(ivar) = lagrange_interpol(order+1,0d0,xi,u(:,i,ie,ivar))
!old         !end do
!old     end if
!old end subroutine


end module
