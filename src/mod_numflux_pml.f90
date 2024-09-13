module mod_numflux_pml

use mod_para
use mod_gll
use mod_fd
use mod_mesh
use mod_funcs

contains

subroutine get_flux_pml(mesh,u,uax,uay,i,is,ie,qi,fstar)
    implicit none
    ! u(Ngrid,Ngrid,nelem,5)
    real(kind=rkind),intent(in) :: u(:,:,:),qi(:,:,:,:)
    real(kind=rkind),intent(in) :: uax(:,:,:), uay(:,:,:)
    !real(kind=rkind),intent(in) :: du(:,:,:)
    integer,intent(in) :: i,is,ie
    integer :: j
    integer :: ief
    type(meshvar),intent(in) :: mesh
    !real(kind=rkind),intent(in) :: rho,cp,cs
    real(kind=rkind) :: rho,cp,cs
    real(kind=rkind) :: rho_out,cp_out,cs_out
    real(kind=rkind),intent(out) :: fstar(:)
    ! local
    integer :: ine,flag,face,ir
    real(kind=rkind),dimension(8) :: uL,uR,uLx,uRx,uLy,uRy,fL,fR,duL,duR,fL1,fL2,fR1,fR2,dQ
    real(kind=rkind) :: nx,ny,norm_vec_n
    real(kind=rkind) :: matT(5,5),invT(5,5),matTv(2,2),invTv(2,2),matTs(3,3),invTs(3,3)
    !real(kind=rkind) :: mu_s,mu_d,Dc,fault_mu,Phi,tau_0,tau_lock,Tau_n,Tn_m,vv,Vel,Tau_str
    real(kind=rkind) :: mu_s,mu_d,Dc,fault_mu,Phi,tau_0,tau_lock,Tau_n,vv,Vel,Tau_str
    real(kind=rkind) :: eta,slip,sliprate
    !real(kind=rkind) :: slip,sliprate
    real(kind=rkind) :: zs_m,zs_p,zp_m,zp_p
    real(kind=rkind) :: C0
    !real(kind=rkind) :: xi(order+1)
    real(kind=rkind) :: x,y
    integer :: mpi_e,mpi_n,n
    real(kind=rkind) :: vec_n(2),vec_m(2)
    real(kind=rkind) :: vx_m,vy_m,Tx_m,Ty_m
    real(kind=rkind) :: vx_p,vy_p,Tx_p,Ty_p
    real(kind=rkind) :: vn_m,vm_m,Tn_m,Tm_m
    real(kind=rkind) :: vn_p,vm_p,Tn_p,Tm_p
    real(kind=rkind) :: dVx,dVy,dTx,dTy
    real(kind=rkind) :: vn_hat_p, vn_hat_m, Tn_hat_p, Tn_hat_m
    real(kind=rkind) :: vm_hat_p, vm_hat_m, Tm_hat_p, Tm_hat_m
    real(kind=rkind) :: FLn,FLm,FL_n,FL_m,FLx,FLy,FL_x,FL_y
    real(kind=rkind) :: FRn,FRm,FR_n,FR_m,FRx,FRy,FR_x,FR_y
    real(kind=rkind) :: miu,lam,chi
    real(kind=rkind) :: damp
    real(kind=rkind) :: alphaS,alphaV
    real(kind=rkind) :: fc(8),fp(8)
    real(kind=rkind),dimension(8) :: Fx,Fy,Fx_out,Fy_out
    real(kind=rkind) :: pbx,pby
    real(kind=rkind) :: pbx_out,pby_out
    real(kind=rkind) :: vecNp(Np)

    damp = 0

    rho = mesh%rho(ie)
    cp  = mesh%vp(ie)
    cs  = mesh%vs(ie)

    ! correct them later
    rho_out = rho
    cp_out = cp
    cs_out = cs

    vecNp = reshape(mesh%pbx(:,:,ie),(/Np/))
    pbx = vecNp(i)
    vecNp = reshape(mesh%pby(:,:,ie),(/Np/))
    pby = vecNp(i)

    pbx_out = pbx;
    pby_out = pby;

    !rrho = 1d0/rho
    miu = rho*cs**2
    chi = rho*cp**2
    lam = chi-2d0*miu

    !xi = mesh%xnode

    zs_m = rho*cs
    zp_m = rho*cp
    zs_p = rho*cs ! is not correct now
    zp_p = rho*cp

    !zp_in  = rho*cp
    !zs_in  = rho*cs

    ! default (for BC_OUT,BC_FREE,...)
    !zp_out = zp_in
    !zs_out = zs_in

    !call get_face(u,i,is,ie,uL)
    uL = u(mesh%vmapM(i,is),ie,:)
    uLx = uax(mesh%vmapM(i,is),ie,:)
    uLy = uay(mesh%vmapM(i,is),ie,:)
    ine = mesh%neigh(is,ie)
    face = mesh%face(is,ie)
    flag = mesh%direction(is,ie)
    !flag = 0
    uRx = 0
    uRy = 0
    ir = i; if (flag > 0) ir = Ngrid-i+1 ! reverse
    !call get_face(u,ir,face,ine,xi,uR)
    !if (ine == ie) then; uR(:) = 0; end if ! absorbing
    if (ine > 0) then
        !call get_face(reshape(u,(/Ngrid,Ngrid,mesh%Nelem,5/)),ir,face,ine,uR)
        if (mesh%ispml(ie)==1) then
            vecNp = reshape(mesh%pbx(:,:,ine),(/Np/))
            pbx_out = vecNp(i)
            vecNp = reshape(mesh%pby(:,:,ine),(/Np/))
            pby_out = vecNp(i)
        else
            pbx_out = pbx
            pby_out = pby
        end if
        uR = u(mesh%vmapP(i,is,ie),ine,:)
        uRx = uax(mesh%vmapP(i,is,ie),ine,:)
        uRy = uay(mesh%vmapP(i,is,ie),ine,:)
        rho_out = mesh%rho(ine)
        cp_out = mesh%vp(ine)
        cs_out = mesh%vs(ine)
        zp_p = mesh%rho(ine) * mesh%vp(ine)
        zs_p = mesh%rho(ine) * mesh%vs(ine)
    elseif (ine == 0) then
        uR(:) = 0 ! absorbing
        duR(:) = 0
        !uR(1:2) =  uL(1:2)
        !uR(3:5) = -uL(3:5)
    elseif (ine == -1) then
        !!! set neighbor media !!!
        !print*,'mpi'
        mpi_e = mesh%mpi_ibool(is,ie)
        mpi_n = mesh%mpi_interface(4,is,ie)
        !print*,'mpi_e=',mpi_e
        !print*,'mpi_n=',mpi_n
        !uR(:) = qi(ir,1:5,mpi_e,mpi_n)  ! mpi
        uR(:) = qi(ir,:,mpi_e,mpi_n)  ! mpi
        duR(:) = 0
        !print*,sngl(uR)
        do j = 1, mesh%pinterfaces
          !Check if the rank and the element in that rank are in the mpi-impedance array (they should)
          if ((abs(mesh%mpi_vp(j,2) - mesh%mpi_interface(1,is,ie)) < epsilon(mesh%mpi_vp(j,2))) .and. &
              (abs(mesh%mpi_vp(j,3) - mesh%mpi_interface(2,is,ie)) < epsilon(mesh%mpi_vp(j,3)))) then
          !if ( int(mesh%mpi_vp(i,2)) .eq. mesh%mpi_interface(1,is,ie) .and. &
          !     int(mesh%mpi_vp(i,3)) .eq. mesh%mpi_interface(2,is,ie) ) then
            !associate corresponding impedance value to zout
            rho_out = mesh%mpi_rho(j,1)
            cp_out = mesh%mpi_vp(j,1)
            cs_out = mesh%mpi_vs(j,1)
            zp_p = rho_out*cp_out
            zs_p = rho_out*cs_out
          end if
        end do
    end if
    if (mesh%bctype(is,ie) == BC_FREE) then ! free surface
        uR(:) = 0
        uR(1:2) =  uL(1:2)
        uR(3) = -0.5*2.0*rho*GRAVITY*mesh%eta(i,is,ie)/chi - uL(3)
        uR(4) = -0.5*2.0*rho*GRAVITY*mesh%eta(i,is,ie)/chi - uL(4)
        !uR(3) = -uL(3)
        !uR(4) = -uL(4)
        uR(5) = -uL(5)
    end if
    !if (mesh%bctype(is,ie) == BC_FREE_G) then ! free surface with gravity
    !    ! u = (rho*Vx,rho*Vy,p/kappa)
    !    uR(:) = 0
    !    uR(1:2) =  uL(1:2) ! u(2,4,5) is not used
    !    uR(3) = 2.0*rho*9.8**mesh%eta(i,is,ie)/chi - uL(3)
    !end if

    duL(:) = 0
    duR(:) = 0

    nx = mesh%nx(i,is,ie)
    ny = mesh%ny(i,is,ie)
    vec_n = (/nx,ny/)
    vec_m = (/-ny,nx/)
    norm_vec_n = mesh%norm_n(i,is,ie)

    !call extract_traction_velocity_strain(uL(1:5),vec_n,rho,    cp,    cs,    vx_m,vy_m,Tx_m,Ty_m)
    !call extract_traction_velocity_strain(uR(1:5),vec_n,rho_out,cp_out,cs_out,vx_p,vy_p,Tx_p,Ty_p)
    call Flux1((uL+uLx)/pbx,duL,rho,cp,cs,Fx)
    call Flux2((uL+uLy)/pby,duL,rho,cp,cs,Fy)
    call Flux1((uR+uRx)/pbx_out,duL,rho_out,cp_out,cs_out,Fx_out)
    call Flux2((uR+uRy)/pby_out,duL,rho_out,cp_out,cs_out,Fy_out)

    ! penalty flux
    !if (mesh%bctype(is,ie) == BC_IN) then
    if (mesh%ispml(ie)==1) then
        fc = nx*(Fx_out-Fx)+ny*(Fy_out-Fy)
        !dQ = uR-uL
        !fc(1) = dTx
        !fc(2) = dTy
        !fc(3) = nx*dVx
        !fc(4) = ny*dVy
        !fc(5) = nx*dVy+ny*dVx

        !fc(1) = nx*dQ(3)+ny*dQ(5)
        !fc(2) = nx*dQ(5)+ny*dQ(4)
        !fc(3) = nx*dQ(1)
        !fc(4) =          ny*dQ(2)
        !fc(5) = nx*dQ(2)+ny*dQ(1)

        fc = 0.5*fc

        ! vx vy sxx syy sxy
        fp(1) = nx*fc(3)+ny*fc(5)
        fp(2) = nx*fc(5)+ny*fc(4)
        fp(3) = nx*fc(1)
        fp(4) =          ny*fc(2)
        fp(5) = nx*fc(2)+ny*fc(1)

        ! stress penalties
        fp(1:2) = fc(1:2) + 0.1*fp(1:2)
        ! velocity penalties
        fp(3:5) = fc(3:5) + 0.1*fp(3:5)

        fp = fc + 0.0 * 0.5 * cp * (uR-uL) !! LLF

        fp = fp * norm_vec_n

        fstar = fp
        fstar(6:8) = 0.0

        return
    end if

end subroutine

end module
