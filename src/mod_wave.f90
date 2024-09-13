module mod_wave

use mod_para
use mod_gll
use mod_fd
use mod_mesh
use mod_funcs
use mod_numflux

contains

subroutine RHS(mesh,u,uax,uay,uax2,uay2,qi,ru,ruax,ruay,ruax2,ruay2)
    implicit none
    type(meshvar),intent(inout) :: mesh
    !type(meshvar),intent(in) :: mesh
    real(kind=rkind),intent(in) :: u(:,:,:)
    real(kind=rkind),intent(in) :: uax(:,:,:)
    real(kind=rkind),intent(in) :: uay(:,:,:)
    real(kind=rkind),intent(in) :: uax2(:,:,:)
    real(kind=rkind),intent(in) :: uay2(:,:,:)
    real(kind=rkind),intent(in) :: qi(:,:,:,:)
    real(kind=rkind),intent(inout) :: ru(:,:,:)
    real(kind=rkind),intent(inout) :: ruax(:,:,:)
    real(kind=rkind),intent(inout) :: ruay(:,:,:)
    real(kind=rkind),intent(inout) :: ruax2(:,:,:)
    real(kind=rkind),intent(inout) :: ruay2(:,:,:)
    integer :: ie,j,i,iv,isub
    real(kind=rkind),dimension(Ngrid,Ngrid,8) :: q,r
    !real(kind=rkind),dimension(5) :: u1,Fx,Fy,Fr,Fs
    real(kind=rkind),dimension(8) :: Fx,Fy
    !real(kind=rkind),dimension(Ngrid) :: rx1,ry1,sx1,sy1,detJ1
    real(kind=rkind),dimension(Ngrid,Ngrid) :: rx,ry,sx,sy,detJ
    real(kind=rkind),dimension(Ngrid,Ngrid,8) :: u2d,uax2d,uay2d,uax22d,uay22d,ru2d,ru2d2
    real(kind=rkind),dimension(Ngrid,8) :: Fline
    !real(kind=rkind),dimension(Ngrid,5,5) :: Ar,As
    !real(kind=rkind),dimension(Ngrid,5,5) :: DA
    !real(kind=rkind),dimension(Ngrid,5) :: Du
    real(kind=rkind),dimension(Ngrid,8) :: temp
    !real(kind=rkind),dimension(5) :: uL,uR,fL,fR,fstar
    real(kind=rkind),dimension(8) :: fstar
    real(kind=rkind) :: rho,cp,cs,chi,lam,miu
    !real(kind=rkind) :: rx,ry,sx,sy,detJ
    !real(kind=rkind) :: vec_n(2),norm_vec_n
    !real(kind=rkind) :: nx,ny,norm_n
    !real(kind=rkind) :: matT(5,5),invT(5,5)
    !real(kind=rkind) :: bs
    integer :: ix(NGLL)
    real(kind=rkind) :: invm
    real(kind=rkind),dimension(NGLL) :: temp1
    real(kind=rkind),dimension(Np) :: pax,pay,pbx,pby,pdx,pdy
    real(kind=rkind),dimension(Np) :: pax2,pay2,pbx2,pby2,pdx2,pdy2
    real(kind=rkind) :: alpha

#ifdef FD
    invm = mesh%invMass(1,1) 
#else
    invm = mesh%invMass(1,1) * dble(Nsub)
    !invm = mesh%invMass(1,1) * 1d0!dble(Nsub)
#endif
    !print*,'invm = ',invm

    !fstar(:) = 1.2d30
    Fline(:,:) = 0

    do ie = 1,mesh%nelem
        ! This loop computes the q = x derivative
        ! sweeps through the rows of each element

        cp = mesh%vp(ie)
        cs = mesh%vs(ie)
        rho = mesh%rho(ie)
        chi = rho*cp*cp
        miu = rho*cs*cs
        lam = chi - 2d0*miu

        !       3
        !    4------3
        !    |      | 2
        !  4 |      | 
        !    1------2
        !       1
        u2d = reshape(u(:,ie,1:8),(/Ngrid,Ngrid,8/))
        uax2d = reshape(uax(:,ie,1:8),(/Ngrid,Ngrid,8/))
        uay2d = reshape(uay(:,ie,1:8),(/Ngrid,Ngrid,8/))
        uax22d = reshape(uax2(:,ie,1:8),(/Ngrid,Ngrid,8/))
        uay22d = reshape(uay2(:,ie,1:8),(/Ngrid,Ngrid,8/))
        if (mesh%ispml(ie)==1) then
        do i=1,8
            !uax2d(:,:,i) = (u2d(:,:,i)+uax2d(:,:,i)+uax22d(:,:,i))/(mesh%pbx(:,:,ie)*mesh%pbx2(:,:,ie))
            !uay2d(:,:,i) = (u2d(:,:,i)+uay2d(:,:,i)+uay22d(:,:,i))/(mesh%pby(:,:,ie)*mesh%pby2(:,:,ie))
            uax2d(:,:,i) = (u2d(:,:,i)+uax2d(:,:,i))/(mesh%pbx(:,:,ie))
            uay2d(:,:,i) = (u2d(:,:,i)+uay2d(:,:,i))/(mesh%pby(:,:,ie))
            !uax2d(:,:,i) = (u2d(:,:,i)+uax2d(:,:,i))/mesh%pbx(:,:,ie)+1.0*uayx2d(:,:,i)
            !uay2d(:,:,i) = (u2d(:,:,i)+uay2d(:,:,i))/mesh%pby(:,:,ie)+1.0*uaxy2d(:,:,i)
        end do
        else
            uax2d = u2d
            uay2d = u2d
        end if
        pax = reshape(mesh%pax(:,:,ie),(/Np/))
        pay = reshape(mesh%pay(:,:,ie),(/Np/))
        pbx = reshape(mesh%pbx(:,:,ie),(/Np/))
        pby = reshape(mesh%pby(:,:,ie),(/Np/))
        pdx = reshape(mesh%pdx(:,:,ie),(/Np/))
        pdy = reshape(mesh%pdy(:,:,ie),(/Np/))
        pax2 = reshape(mesh%pax2(:,:,ie),(/Np/))
        pay2 = reshape(mesh%pay2(:,:,ie),(/Np/))
        pbx2 = reshape(mesh%pbx2(:,:,ie),(/Np/))
        pby2 = reshape(mesh%pby2(:,:,ie),(/Np/))
        pdx2 = reshape(mesh%pdx2(:,:,ie),(/Np/))
        pdy2 = reshape(mesh%pdy2(:,:,ie),(/Np/))
        ru2d = reshape(ru(:,ie,1:8),(/Ngrid,Ngrid,8/)) * 0 * mesh%dtfactor
        ! vel + eta*dt*acc
        !u2d(:,:,1:2) = u2d(:,:,1:2) + 0.3 * mesh%dtfactor * ru2d(:,:,1:2)
        rx = mesh%rx(:,:,ie)
        ry = mesh%ry(:,:,ie)
        sx = mesh%sx(:,:,ie)
        sy = mesh%sy(:,:,ie)
        detJ = mesh%detJ(:,:,ie)


        ! Line is in the xi direction
        do j = 1,Ngrid
            do i = 1,Ngrid
#ifdef VERSION1
                if (mesh%elemtype(ie) == ELEM_SOLID) then
#endif
                    call Flux1(uax2d(i,j,:),ru2d(i,j,:),rho,cp,cs,Fx)
                    call Flux2(uay2d(i,j,:),ru2d(i,j,:),rho,cp,cs,Fy)
#ifdef VERSION1
                else if (mesh%elemtype(ie) == ELEM_FLUID) then
                    call Flux1_acoustic(u2d(i,j,:),rho,cp,Fx)
                    call Flux2_acoustic(u2d(i,j,:),rho,cp,Fy)
                else
                    print*,'error of elemtype'
                end if
#endif
                !Fx = (Fx+uax2d(i,j,:))/mesh%pbx(i,j,ie)
                !Fy = (Fy+uay2d(i,j,:))/mesh%pby(i,j,ie)
                Fline(i,:) = ( rx(i,j)*Fx(:) + ry(i,j)*Fy(:) ) * detJ(i,j)
                !Fline = Fline + uax2d(i,j,:) * detJ(i,j)
            end do
#ifdef FD
            do iv = 1,5
                call fd_d1(mesh%fd_D,FDorder,Fline(:,iv),temp(:,iv)) ! strong form
            enddo
#else
            !temp = matmul(mesh%D,Fline) ! strong form
            temp = 0d0
            do isub = 1,Nsub
                ix = mesh%iglob(:,isub)
                !temp(ix,:) = temp(ix,:) + matmul(mesh%D,Fline(ix,:)) * dble(Nsub)
                do iv = 1,5
                    call mxm(mesh%Dt,Fline(ix,iv),temp1)
                    temp(ix,iv) = temp(ix,iv) + temp1 * dble(Nsub)
                end do
                !temp(ix,:) = temp(ix,:) + matmul(mesh%Filter,matmul(mesh%D,Fline(ix,:))) * dble(Nsub)
            end do
            do i = 1,Ngrid
                temp(i,:) = temp(i,:) * mesh%globInvM(i)
            end do
#endif

            call get_flux(mesh,u,ru,uax,uay,j,4,ie,qi,fstar)
            temp(1,:) = temp(1,:) + invm * fstar
            call get_flux(mesh,u,ru,uax,uay,j,2,ie,qi,fstar)
            temp(Ngrid,:) = temp(Ngrid,:) + invm * fstar

            if (.false.) then
            do isub = 1,Nsub
                ix = mesh%iglob(:,isub)
                temp(ix,:) = matmul(mesh%Filter,temp(ix,:))
            end do
            end if

            q(:,j,1:8) = temp ! (n,n) x (n,5)
        end do

        ! Line is in the eta direction
        do i = 1,Ngrid
            do j = 1,Ngrid
#ifdef VERSION1
                if (mesh%elemtype(ie) == ELEM_SOLID) then
#endif
                    call Flux1(uax2d(i,j,:),ru2d(i,j,:),rho,cp,cs,Fx)
                    call Flux2(uay2d(i,j,:),ru2d(i,j,:),rho,cp,cs,Fy)
#ifdef VERSION1
                else if (mesh%elemtype(ie) == ELEM_FLUID) then
                    call Flux1_acoustic(u2d(i,j,:),rho,cp,Fx)
                    call Flux2_acoustic(u2d(i,j,:),rho,cp,Fy)
                else
                    print*,'error of elemtype'
                end if
#endif
                !Fx = (Fx+uax2d(i,j,:))/mesh%pbx(i,j,ie)
                !Fy = (Fy+uay2d(i,j,:))/mesh%pby(i,j,ie)
                Fline(j,:) = ( sx(i,j)*Fx(:) + sy(i,j)*Fy(:) ) * detJ(i,j)
                !Fline = Fline + uay2d(i,j,:) * detJ(i,j)
            end do
#ifdef FD
            do iv = 1,5
                call fd_d1(mesh%fd_D,FDorder,Fline(:,iv),temp(:,iv)) ! strong form
            enddo
#else
            !temp = matmul(mesh%D,Fline) ! (n,5)
            temp = 0d0
            do isub = 1,Nsub
                ix = mesh%iglob(:,isub)
                !temp(ix,:) = temp(ix,:) + matmul(mesh%D,Fline(ix,:)) * dble(Nsub)
                do iv = 1,5
                    call mxm(mesh%Dt,Fline(ix,iv),temp1)
                    temp(ix,iv) = temp(ix,iv) + temp1 * dble(Nsub)
                end do
                !temp(ix,:) = temp(ix,:) + matmul(mesh%Filter,matmul(mesh%D,Fline(ix,:))) * dble(Nsub)
            end do
            do j = 1,Ngrid
                temp(j,:) = temp(j,:) * mesh%globInvM(j)
            end do
#endif

            call get_flux(mesh,u,ru,uax,uay,i,1,ie,qi,fstar)
            temp(1,:) = temp(1,:) + invm * fstar
            call get_flux(mesh,u,ru,uax,uay,i,3,ie,qi,fstar)
            temp(Ngrid,:) = temp(Ngrid,:) + invm * fstar

            if (.false.) then
            do isub = 1,Nsub
                ix = mesh%iglob(:,isub)
                temp(ix,:) = matmul(mesh%Filter,temp(ix,:))
            end do
            end if

            r(i,:,1:8) = temp ! (n,n) x (n,5)
        end do

        do iv = 1,5
            ru2d(:,:,iv) = (q(:,:,iv)+r(:,:,iv))/mesh%detJ(:,:,ie) 
        end do

#ifdef SYM
        ru2d2(:,:,1) = ru2d(:,:,1)/rho
        ru2d2(:,:,2) = ru2d(:,:,2)/rho
        ru2d2(:,:,3) = chi*ru2d(:,:,3) + lam*ru2d(:,:,4)
        ru2d2(:,:,4) = lam*ru2d(:,:,3) + chi*ru2d(:,:,4)
        ru2d2(:,:,5) = miu*ru2d(:,:,5)
        ru2d2(:,:,6:8) = 0.

        ru2d = ru2d2
#endif

        ru(:,ie,1:8) = reshape(ru2d,(/Np,8/))
        !ru(:,ie,6:8) = 20*mesh%dtfactor * (u(:,ie,6:8)-u(:,ie,3:5))
        ru(:,ie,6:8) = 0.1*mesh%dtfactor * ru(:,ie,3:5)
        ru(:,ie,6:8) = 0

        alpha = 0.3

        if (mesh%ispml(ie)==1) then
        do i = 1,Np
            !call Flux1(u(i,ie,:),0*u(i,ie,:),rho,cp,cs,Fx)
            !call Flux2(u(i,ie,:),0*u(i,ie,:),rho,cp,cs,Fy)
        do iv = 1,5
            !alpha = 1-pdx2(i)*pbx(i)/(pbx(i)*(pdx2(i)+pax2(i)*pbx2(i))-pbx2(i)*(pdx(i)+pax(i)*pbx(i)))
            alpha = 1.
            if (alpha/=alpha) print*,'alpha=',alpha
            ruax(i,ie,iv) = -(pax(i)+pdx(i)/pbx(i))*uax(i,ie,iv)-(pdx(i)/pbx(i))*u(i,ie,iv)*alpha
            !alpha = 1-pdy2(i)*pby(i)/(pby(i)*(pdy2(i)+pay2(i)*pby2(i))-pby2(i)*(pdy(i)+pay(i)*pby(i)))
            alpha = 1.
            if (alpha/=alpha) print*,'alpha=',alpha
            ruay(i,ie,iv) = -(pay(i)+pdy(i)/pby(i))*uay(i,ie,iv)-(pdy(i)/pby(i))*u(i,ie,iv)*alpha

            !alpha = 1+pdx(i)*pbx2(i)/(pbx(i)*(pdx2(i)+pax2(i)*pbx2(i))-pbx2(i)*(pdx(i)+pax(i)*pbx(i)))
            !if (alpha/=alpha) print*,'alpha=',alpha
            !ruax2(i,ie,iv) = -(pax2(i)+pdx2(i)/pbx2(i))*uax2(i,ie,iv)-(pdx2(i)/pbx2(i))*u(i,ie,iv)*alpha
            !alpha = 1+pdy(i)*pby2(i)/(pby(i)*(pdy2(i)+pay2(i)*pby2(i))-pby2(i)*(pdy(i)+pay(i)*pby(i)))
            !if (alpha/=alpha) print*,'alpha=',alpha
            !ruay2(i,ie,iv) = -(pay2(i)+pdy2(i)/pby2(i))*uay2(i,ie,iv)-(pdy2(i)/pby2(i))*u(i,ie,iv)*alpha
            !ruax(i,ie,iv) = -(pax(i)+pdx(i)/pbx(i))*uax(i,ie,iv)-(pdx(i)/pbx(i))*Fx(iv)
            !ruay(i,ie,iv) = -(pay(i)+pdy(i)/pby(i))*uay(i,ie,iv)-(pdy(i)/pby(i))*Fy(iv)
            !ruaxy(:,ie,iv) = -alpha*pdy*uaxy(:,ie,iv)-alpha*pdy*u(:,ie,iv)
            !ruayx(:,ie,iv) = -alpha*pdx*uayx(:,ie,iv)-alpha*pdx*u(:,ie,iv)
        end do
        end do
        end if

        ! these terms can be ignored
!!!        ! add gravity
!!!        if (mesh%elemtype(ie) == ELEM_FLUID) then
!!!            !print*,'add body force ...'
!!!            ! u = (rho*vx rho*vy p/kappa)
!!!            ! u = (rho*vx rho*vy exx eyy exy)
!!!#ifdef VERSION1
!!!            ru(:,ie,2) = ru(:,ie,2) - GRAVITY*rho*(u(:,ie,3))
!!!            ru(:,ie,3) = ru(:,ie,3) + GRAVITY/chi*u(:,ie,2)
!!!#else
!!!            ru(:,ie,2) = ru(:,ie,2) - GRAVITY*rho*(u(:,ie,3)+u(:,ie,4))*(-1)
!!!            ru(:,ie,3) = ru(:,ie,3) + GRAVITY/chi*u(:,ie,2)
!!!            ru(:,ie,4) = ru(:,ie,4) + GRAVITY/chi*u(:,ie,2)
!!!#endif
!!!        end if

    end do ! elem

end subroutine

subroutine write_wave(u,rank,flag)
    implicit none
    real(kind=rkind) :: u(:,:,:)
    integer,intent(in) :: rank,flag
    character(len=128) :: filename
    !write(filename,'(a,i6.6,a)') 'data/wave_it',it,'.bin'
    !open(100,file=trim(filename),access='direct',form='unformatted',recl=sizeof(1.0)*size(u),status='replace')
    !write(100,rec=1) sngl(u)

    write(filename,'(a,i6.6)') 'data/wave',rank
    if (flag == 0) then
        open(100,file=trim(filename),access='stream',form='unformatted',status='replace',action='write')
    else
        open(100,file=trim(filename),access='stream',form='unformatted',status='unknown',position='append')
    end if
    write(100) sngl(u)
    close(100)
end subroutine write_wave

subroutine init_wave(mesh,u)
    implicit none
    type(meshvar),intent(in) :: mesh
    real(kind=rkind),intent(out) :: u(:,:,:)
    real(kind=rkind) :: xmin,xmax,ymin,ymax
    real(kind=rkind) :: srcx,srcy,a
    real(kind=rkind),allocatable,dimension(:,:) :: r,amp
    integer :: i,j,ie,Nelem

    Nelem = mesh%Nelem

    allocate(r(Np,mesh%nelem))
    allocate(amp(Np,mesh%nelem))

    xmin = minval(mesh%vx)
    xmax = maxval(mesh%vx)
    ymin = minval(mesh%vy)
    ymax = maxval(mesh%vy)

    !print*,'Physical domian:'
    !print*,'x = ',sngl(xmin),'~',sngl(xmax)
    !print*,'y = ',sngl(ymin),'~',sngl(ymax)
    !write(*,'(a,f8.2,a,f8.2,a,f8.2,a,f8.2)') 'x=(',xmin,',',xmax,'); y=(',ymin,',',ymax,')'
    print*, 'x=(',sngl(xmin),',',sngl(xmax),'); y=(',sngl(ymin),',',sngl(ymax),')'

    !print*,sngl(mesh%vx(:,1))
    !print*,sngl(mesh%vy(:,1))
    srcx = (xmin+xmax)/2.0;
    srcy = (ymin+ymax)/2.0;
    srcx = -5*1e3
    srcy = -5*1e3
    srcx = -20e3
    srcx = 0
    srcy = 0
    !srcx = 5e3
    !srcy = 5e3
    r = (reshape(mesh%vx,(/Np,Nelem/))-srcx)**2+(reshape(mesh%vy,(/Np,Nelem/))-srcy)**2
    a = min(xmax-xmin,ymax-ymin)/20.0
    a = 2e3
    amp = 0
    !amp = 1*exp(-r/a**2)
    !mesh%bctype = 0

    !do ie = 1,mesh%nelem
    !do j = 1,NGLL
    !do i = 1,NGLL
    !if (mesh%vy(i,j,ie)<0) then
    !amp(i,j,ie) = 1*exp(-r(i,j,ie)/a**2)
    !end if
    !enddo
    !enddo
    !enddo

    !mesh%bctype = 0

    ! rhoVx rhoVy Exx Eyy Exy
    !u(:,:,1) = 0
    !u(:,:,2) = 0
    !u(:,:,3) = amp*1
    !u(:,:,4) = amp*1
    !u(:,:,5) = amp*0

    u = 0
    do ie = 1,mesh%nelem
    if (mesh%elemtype(ie) == ELEM_SOLID) then
    u(:,ie,1) = 0
    u(:,ie,2) = 0
    u(:,ie,3) = amp(:,ie)*1
    u(:,ie,4) = amp(:,ie)*1
    u(:,ie,5) = amp(:,ie)*0
    end if
    !if (mesh%elemtype(ie) == ELEM_FLUID) then
    !u(:,ie,1) = 0
    !u(:,ie,2) = 0
    !u(:,ie,3) = amp(:,ie)*1
    !u(:,ie,4) = 0
    !u(:,ie,5) = 0
    !end if

    !if (mesh%elemtype(ie) == ELEM_FLUID) then
    !u(:,ie,1) = 0
    !u(:,ie,2) = 0
    !u(:,ie,3) = amp(:,ie)*1
    !u(:,ie,4) = 0
    !u(:,ie,5) = 0
    !end if
    end do
end subroutine

end module
