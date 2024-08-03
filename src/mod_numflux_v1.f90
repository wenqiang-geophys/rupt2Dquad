module mod_numflux

use mod_para
use mod_gll
use mod_fd
use mod_mesh
use mod_funcs

contains

subroutine get_flux(mesh,u,du,i,is,ie,qi,fstar)
    implicit none
    ! u(Ngrid,Ngrid,nelem,5)
    real(kind=rkind),intent(in) :: u(:,:,:),qi(:,:,:,:)
    real(kind=rkind),intent(in) :: du(:,:,:)
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
    real(kind=rkind),dimension(8) :: uL,uR,fL,fR,duL,duR,fL1,fL2,fR1,fR2,dQ
    real(kind=rkind) :: nx,ny,norm_vec_n
    real(kind=rkind) :: matT(5,5),invT(5,5),matTv(2,2),invTv(2,2),matTs(3,3),invTs(3,3)
    !real(kind=rkind) :: mu_s,mu_d,Dc,fault_mu,Phi,tau_0,tau_lock,Tau_n,Tn_m,vv,Vel,Tau_str
    real(kind=rkind) :: mu_s,mu_d,Dc,fault_mu,Phi,tau_0,tau_lock,Tau_n,vv,Vel,Tau_str
    real(kind=rkind) :: zp_in,zs_in,zp_out,zs_out,eta,slip,sliprate
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
    real(kind=rkind) :: rrho,miu,lam,chi
    real(kind=rkind) :: damp
    real(kind=rkind) :: alphaS,alphaV
    real(kind=rkind) :: fc(8),fp(8)

    damp = 0

    rho = mesh%rho(ie)
    cp  = mesh%vp(ie)
    cs  = mesh%vs(ie)

    ! correct them later
    rho_out = rho
    cp_out = cp
    cs_out = cs

    rrho = 1d0/rho
    miu = rho*cs**2
    chi = rho*cp**2
    lam = chi-2d0*miu

    !xi = mesh%xnode

    zs_m = rho*cs
    zs_p = rho*cs ! is not correct now
    zp_m = rho*cp
    zp_p = rho*cp

    zp_in  = rho*cp
    zs_in  = rho*cs

    ! default (for BC_OUT,BC_FREE,...)
    zp_out = zp_in
    zs_out = zs_in

    !call get_face(u,i,is,ie,uL)
    uL = u(mesh%vmapM(i,is),ie,:)
    duL = du(mesh%vmapM(i,is),ie,:) * damp * mesh%dtfactor
    ine = mesh%neigh(is,ie)
    face = mesh%face(is,ie)
    flag = mesh%direction(is,ie)
    !flag = 0
    uR = 0
    duR = 0
    ir = i; if (flag > 0) ir = Ngrid-i+1 ! reverse
    !call get_face(u,ir,face,ine,xi,uR)
    !if (ine == ie) then; uR(:) = 0; end if ! absorbing
    if (ine > 0) then
        !call get_face(reshape(u,(/Ngrid,Ngrid,mesh%Nelem,5/)),ir,face,ine,uR)
        uR = u(mesh%vmapP(i,is,ie),ine,:)
        duR = du(mesh%vmapP(i,is,ie),ine,:) * damp * mesh%dtfactor
        rho_out = mesh%rho(ine)
        cp_out = mesh%vp(ine)
        cs_out = mesh%vs(ine)
        zp_out = mesh%rho(ine) * mesh%vp(ine)
        zs_out = mesh%rho(ine) * mesh%vs(ine)
    elseif (ine == 0) then
        uR(:) = 0 ! absorbing
        duR(:) = 0
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
            !rrho_out = 1d0/rho_out
          end if
        end do
    end if
    if (mesh%bctype(is,ie) == BC_FREE) then ! free surface
        uR(1:2) =  uL(1:2)
        uR(3:5) = -uL(3:5)
    end if
    if (mesh%elemtype(ie) == ELEM_FLUID .and. mesh%bctype(is,ie) == BC_FREE) then ! free surface
        ! v+ = v-
        ! p+ = 2*rho*g*eta - p-
        ! p+/kappa = 2*rho*g*eta/kappa - p-/kappa
        uR(1) =  uL(1)
        uR(3) = 2.0*rho*GRAVITY*mesh%eta(i,is,ie)/chi - uL(3)
    end if
    if (mesh%bctype(is,ie) == BC_FREE_G) then ! free surface with gravity
        ! u = (rho*Vx,rho*Vy,p/kappa)
        uR(:) = 0
        uR(1:2) =  uL(1:2) ! u(2,4,5) is not used
        uR(3) = 2.0*rho*9.8**mesh%eta(i,is,ie)/chi - uL(3)
    end if

    duL(:) = 0
    duR(:) = 0

    nx = mesh%nx(i,is,ie)
    ny = mesh%ny(i,is,ie)
    norm_vec_n = mesh%norm_n(i,is,ie)

    !if (mesh%bctype(is,ie) == BC_IN .and. ine .ne. 0 ) then
    !    !fstar = 0.5*fstar + 0.5*(fR-fL)/2.0 ! central flux
    !    call Flux1(uL,duL,rho,cp,cs,fL1)
    !    call Flux1(uR,duR,rho,cp,cs,fR1)
    !    call Flux2(uL,duL,rho,cp,cs,fL2)
    !    call Flux2(uR,duR,rho,cp,cs,fR2)
    !    fstar = 0.5*(fR1*nx-fL1*nx+fR2*ny-fL2*ny)  ! LLF flux
    !    fstar = fstar + 0.5*0.5*cp*(uR-uL)  ! LLF flux
    !    fstar = fstar * norm_vec_n
    !    return
    !end if

    invT = getInvT(nx,ny)
    matT = getT   (nx,ny)
    invTv = invT(1:2,1:2)
    invTs = invT(3:5,3:5)
    matTv = matT(1:2,1:2)
    matTs = matT(3:5,3:5)
    !uL = matmul(invT,uL)
    !uR = matmul(invT,uR)
    if (mesh%elemtype(ie) == ELEM_SOLID) then
        !print*,'solid elem'
        call rotate_u(invTv,invTs,uL)
        call rotate_u(invTv,invTs,duL)
        call Flux1(uL,duL,rho,cp,cs,fL)
        ! deal with neighbors
        if (mesh%bctype(is,ie) == BC_SOLID_SOLID .or. &
                mesh%bctype(is,ie) == BC_FREE .or. &
                mesh%bctype(is,ie) == BC_FAULT  &
                 ) then
            call rotate_u(invTv,invTs,uR)
            call rotate_u(invTv,invTs,duR)
            call Flux1(uR,duR,rho_out,cp_out,cs_out,fR)
            call numflux1(fL,fR,zp_in,zp_out,zs_in,zs_out,fstar)
        else if (mesh%bctype(is,ie) == BC_SOLID_FLUID) then
            call rotate_u_acoustic(invTv,uR)
            call Flux1_acoustic(uR,rho_out,cp_out,fR)
            call numflux1_elastic_acoustic(fL,fR,zp_in,zp_out,zs_in,zs_out,fstar)
        else
            print*,'error of bctype (ELEM SOLID)'
        end if
        !fstar = matmul(matT, fstar)
        !call rotate_u(matTv,matTs,fstar)
        !fstar = fstar * norm_vec_n
    elseif (mesh%elemtype(ie) == ELEM_FLUID) then
        call rotate_u_acoustic(invTv,uL)
        call Flux1_acoustic(uL,rho,cp,fL)
        ! deal with neighbors
        if ( &
                mesh%bctype(is,ie) == BC_FLUID_FLUID .or. &
                mesh%bctype(is,ie) == BC_FREE .or. &
                mesh%bctype(is,ie) == BC_FREE_G  &
                ) then
            call rotate_u_acoustic(invTv,uR)
            call Flux1_acoustic(uR,rho_out,cp_out,fR)
            call numflux1_acoustic(fL,fR,zp_in,zp_out,fstar)
            if (mesh%bctype(is,ie) == BC_FREE) then
                ! f = -(p,0,vx)  u = (rho*vx,rho*vy,p/kappa)
                mesh%deta(i,is,ie) = -(fL(3)+fstar(3))*1.0 ! v*
            else
                mesh%deta(i,is,ie) = 0
            end if
        else if ( mesh%bctype(is,ie) == BC_SOLID_FLUID) then
            call rotate_u(invTv,invTs,uR)
            call rotate_u(invTv,invTs,duR)
            call Flux1(uR,duR,rho_out,cp_out,cs_out,fR)
            call numflux1_acoustic_elastic(fL,fR,zp_in,zp_out,fstar)
        else
            print*,'error of bctype (ELEM FLUID)'!,'bctype=',mesh%bctype(is,ie)
        end if
    else
        print*,'error of elemtype'
    end if

    !if (mesh%bctype(is,ie) == BC_IN .and. ine .ne. 0 ) then
    !    !fstar = 0.5*fstar + 0.5*(fR-fL)/2.0 ! central flux
         !call numflux1_mod(fL,fR,zp_in,zp_out,zs_in,zs_out,0.1d0,fstar)
    !    fstar = 0.5*cp*(uR-uL) + 1.0*(fR+fL)/2.0 - fL ! LLF flux
    !end if
    !if (mesh%bctype(is,ie) == BC_FAULT ) then
    !    alphaV = 1.0
    !    alphaS = 1.0
    !else
    !    alphaV = 0.5
    !    alphaS = 0.5
    !end if
    !fstar(1) = 0.5d0*(fR(1)-fL(1)) + alphaV*0.5d0*(fR(3)-fL(3))*Zp_in
    !fstar(2) = 0.5d0*(fR(2)-fL(2)) + alphaV*0.5d0*(fR(5)-fL(5))*Zs_in
    !fstar(3) = 0.5d0*(fR(3)-fL(3)) + alphaS*0.5d0*(fR(1)-fL(1))/Zp_in
    !fstar(4) = 0.0d0
    !fstar(5) = 0.5d0*(fR(5)-fL(5)) + alphaS*0.5d0*(fR(2)-fL(2))/Zs_in

    if (mesh%bctype(is,ie) == BC_FAULT) then
        ief = mesh%wave2fault(ie)
        !if (is == 1) then
        !    x = mesh%vx(i,1,ie)
        !    y = mesh%vy(i,1,ie)
        !else if (is == 2) then
        !    x = mesh%vx(NGLL,i,ie)
        !    y = mesh%vy(NGLL,i,ie)
        !else if (is == 3) then
        !    x = mesh%vx(i,NGLL,ie)
        !    y = mesh%vy(i,NGLL,ie)
        !else if (is == 4) then
        !    x = mesh%vx(1,i,ie)
        !    y = mesh%vy(1,i,ie)
        !else
        !    y = 1e30
        !end if

        !if (abs(x+0.e3) < 1.5e-3) then
        !    print*,sngl(x),sngl(y),sngl(fstar)
        !    stop 2
        !end if

        !Tau_n = +120e6
        !Tau_n = y/sqrt(3.0)*2.0*7378.0
        zs_in  = rho*cs
        zs_out = rho*cs
        eta = zs_in*zs_out/(zs_in+zs_out)
        !print*,eta
        Tau_n    = mesh%Tau_n   (i,is,ief)
        Tau_0    = mesh%Tau_0   (i,is,ief)
        slip     = mesh%Slip    (i,is,ief)
        sliprate = mesh%Sliprate(i,is,ief)
        mu_s     = mesh%mu_s    (i,is,ief)
        mu_d     = mesh%mu_d    (i,is,ief)
        Dc       = mesh%Dc      (i,is,ief)
        C0       = mesh%C0      (i,is,ief)

        ! f (sxx,sxy,vx,0,vy)
        !Phi = fstar(2) ! sxy stress transfer function
        Tau_lock = fL(2) + fstar(2) + Tau_0
        Tn_m = 0.0
        Tn_m = fL(1) + fstar(1)  ! sxx
        !Tn_m = sign(1d0,Tau_lock) * Tn_m


        fault_mu = mu_s - (mu_s-mu_d) * (min(abs(slip),Dc)/Dc)
        ! instantanenous healing
        !if (.false.) then
        !if (abs(sliprate)<1e-10) then
        !    fault_mu = mu_s
        !    slip = 0
        !    mesh%slip(i,is,ie) = 0
        !end if
        !end if


        !Tau_str = fault_mu * max(0.0,+Tn_m-Tau_n) + C0
        Tau_str = fault_mu * max(0.0,-(Tn_m+Tau_n)) + C0
        !Tau_str = fault_mu * max(0.0,-Tau_n) + C0
        if (abs(Tau_lock) .gt. abs(Tau_str)) then ! fault is slipping
        !if (Tau_lock .ge. Tau_str) then ! fault is slipping
            ! solve for slip rate
            !vv = (Tau_lock-Tau_str)/eta
            !vv = (abs(Tau_lock)-Tau_str)/eta * sign(1d0,Tau_lock)
            Vel = (abs(Tau_lock)-abs(Tau_str))/eta
            vv = Tau_lock*Vel/(eta*Vel+Tau_str) ! same sign with Tau_lock, parallel condition
            !vv = ((Phi+Tau_0)-Tau_str)/eta
            !vv = sqrt(vv**2+1e-6**2)-1e-6
        else
            vv = 0.0
        endif
        !print*,vv
        !mesh%sliprate(i,is,ie) = abs(vv)!abs(V_p-V_m)
        mesh%sliprate(i,is,ief) = abs(vv) !abs(V_p-V_m)
        !mesh%stress(i,is,ief) = Tau_lock - eta*vv
        mesh%stress(i,is,ief) = Tm_hat_m+Tau_0
        !mesh%recv_buffer(n,i,2) = Tm_hat_m+Tau_0
        !mesh%recv_buffer(n,i,3) = Tn_hat_m+Tau_n


        ! f (sxx,sxy,vx,0,vy)
        fstar(2) = fstar(2) - eta*vv       ! sxy
        fstar(5) = fstar(5) - eta*vv/zs_in ! vy

        !if ( nx>0 )then
        !if ( &
        !        abs(y+0.0e3*dsqrt(3d0)/2d0) < 1e1 .or. &
        !        abs(y+1.5e3*dsqrt(3d0)/2d0) < 1e1 .or. &
        !        abs(y+3.0e3*dsqrt(3d0)/2d0) < 1e1 .or. &
        !        abs(y+4.5e3*dsqrt(3d0)/2d0) < 1e1 .or. &
        !        abs(y+7.5e3*dsqrt(3d0)/2d0) < 1e1 .or. &
        !        abs(y+12.e3*dsqrt(3d0)/2d0) < 1e1 &
        !        ) then
        !    write(*,'(i3,i8,f10.2,5e16.4)') is,ie,-sngl(y)/dsqrt(3d0)*2d0*1d-3,sngl(nx),&
        !            sngl(mesh%current_time),sngl(vv),&
        !            sngl((fstar(2)+Tau_0)),&
        !            sngl((fstar(1)+Tau_n))
        !end if
        !end if

        ! save current time fault recvs
        do n = 1,mesh%nrecv
            !if (ie == mesh%recv_ie(n) .and. i == mesh%recv_i(n) )  then
            if (ie == mesh%recv_ie(n))  then
                !print*,ie,i,sngl(mesh%current_time),sngl(vv)
                mesh%recv_buffer(n,i,1) = abs(vv)
                mesh%recv_buffer(n,i,2) = fL(2) + fstar(2) + Tau_0
                mesh%recv_buffer(n,i,3) = fL(1) + fstar(1) + Tau_n
                !print*,n,ie,i,mesh%recv_fid(n),sngl(mesh%recv_buffer(n,1:3))
            end if
        enddo
    end if ! bc = fault

   !fault     !print*,fstar(5),nx,ny,norm_vec_n
   !fault end if
    !fstar = matmul(matT, fstar) * norm_vec_n
    !fstar = matmul(matT, fstar)
    if (mesh%elemtype(ie) == ELEM_SOLID) then
        call rotate_u(matTv,matTs,fstar)
    elseif (mesh%elemtype(ie) == ELEM_FLUID) then
        call rotate_u_acoustic(matTv,fstar)
    else
        print*,'error of elemtype'
    end if
    fstar = fstar * norm_vec_n
    !if (mesh%bctype(is,ie) == BC_FAULT) then
    !print*,fstar(5)
end subroutine

end module
