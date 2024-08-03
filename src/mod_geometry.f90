module mod_geometry
use mod_para
use mod_gll
use mod_mesh
use mod_fd
!use mod_wave
!use mod_jacobi

contains

subroutine build_geometry(mesh)
    implicit none
    type(meshvar) :: mesh
    !real(kind=rkind),dimension(order+1) :: r1,s1
    real(kind=rkind),dimension(Ngrid,Ngrid) :: r,s
    real(kind=rkind),dimension(Ngrid,Ngrid) :: rx,sx,ry,sy,detJ
    real(kind=rkind),dimension(Ngrid,Ngrid) :: xr,xs,yr,ys
    real(kind=rkind),dimension(Ngrid,nsurface) :: nx,ny,norm_n
    !real(kind=rkind),dimension(Np,Np) :: Vdm,Dr,Ds
    !real(kind=rkind),dimension(Np,Np) :: Drw,Dsw,Vdmr,Vdms,tempM1,tempM2
    integer :: i,j,N,is,ie,face,direction
    !integer :: iv(Np)
    integer :: nelem
    integer :: p1,p2
    real(kind=rkind),allocatable :: x1(:),x2(:),wts1(:),wts2(:)
    real(kind=rkind),allocatable :: G(:,:),D(:,:),G1(:,:),D1(:,:),W(:,:),Mass(:,:),invMass(:,:),V(:,:),iV(:,:)
    real(kind=rkind),allocatable :: Dp(:,:),Dm(:,:)
    real(kind=rkind) :: Filter(NGLL,NGLL),Ccut(NGLL,NGLL),u0(NGLL),u1(NGLL),Gd(NGLL-1,NGLL),Gu(NGLL,NGLL-1),xd(NGLL-1),wtsd(NGLL-1)
    real(kind=rkind) :: VX(4),VY(4)
    !real(kind=rkind) :: nx,ny,norm_n
    integer :: iv1d(Np),iv2d(Ngrid,Ngrid)
    real(kind=rkind) :: v2d(Ngrid,Ngrid)
    real(kind=rkind) :: eps,alpha
    real(kind=rkind) :: h
    integer :: sp,Nc
    !real(kind=rkind) :: tmpface(NGLL)
    integer :: ix(NGLL)

    Nelem = mesh%Nelem

    !       3
    !    4------3
    !    |      | 2
    !  4 |      | 
    !    1------2
    !       1
    N = Ngrid
    do i = 1,Np; iv1d(i) = i; end do
    iv2d = reshape(iv1d,(/N,N/))
    do i = 1,Ngrid
        do j = 1,Ngrid
            v2d(i,j) = iv2d(i,j)
        end do
    end do
    mesh%vmapM(:,1) = iv2d(:,1)
    mesh%vmapM(:,2) = iv2d(N,:)
    mesh%vmapM(:,3) = iv2d(:,N)
    mesh%vmapM(:,4) = iv2d(1,:)

    allocate(mesh%vmapP(Nfp,Nfaces,Nelem))
    do ie = 1,Nelem
        do is = 1,Nfaces
            face = mesh%face(is,ie)
            direction = mesh%direction(is,ie)

            if (direction == 0) then
                mesh%vmapP(:,is,ie) = mesh%vmapM(:,face)
            else if (direction == 1) then
                mesh%vmapP(:,is,ie) = mesh%vmapM(N:1:-1,face)
            else
                print*,'Error at direction'
            end if
            !print*,ie,is,mesh%vmapM(:,is)!,'vmapP=',mesh%vmapP(:,is,ie)
            !call get_face1(v2d,is,tmpface)!,'vmapP=',mesh%vmapP(:,is,ie)
            !if (direction == 1) tmpface = tmpface(N:1:-1)
            !print*,ie,is,sngl(tmpface)-mesh%vmapM(:,is)!,'vmapP=',mesh%vmapP(:,is,ie)

            !print*,ie,is,mesh%vmapP(:,is,ie)!,'vmapP=',mesh%vmapP(:,is,ie)
            !call get_face1(v2d,face,tmpface)!,'vmapP=',mesh%vmapP(:,is,ie)
            !if (direction == 1) tmpface = tmpface(N:1:-1)
            !print*,ie,is,sngl(tmpface)-mesh%vmapP(:,is,ie)!,'vmapP=',mesh%vmapP(:,is,ie)
            !print*,ie,is,direction
        end do
    end do

    mesh%globInvM(:) = 0
    do i = 1,Nsub
        do j = 1,NGLL
            ix(j) = (i-1)*(NGLL-1)+j
        end do

        mesh%globInvM(ix) = mesh%globInvM(ix) + 1.0
        mesh%iglob(:,i) = ix
    end do

    do i = 1,Ngrid
        mesh%globInvM(i) = 1.0/mesh%globInvM(i)
    end do

    !print*,mesh%iglob
    !print*,mesh%globInvM
    !stop 2

    
    p1 = order
    p2 = order2
    
#ifdef FD
    allocate(x1(ngrid))
    allocate(x2(ngrid))
    allocate(wts1(ngrid))
    allocate(wts2(ngrid))
    allocate(G(ngrid,ngrid))
    allocate(Dp(ngrid,ngrid))
    allocate(Dm(ngrid,ngrid))
    allocate(D(ngrid,ngrid))
    allocate(D1(ngrid,ngrid))
    allocate(G1(ngrid,ngrid))
    allocate(W(ngrid,ngrid))
    allocate(V(ngrid,ngrid))
    allocate(iV(ngrid,ngrid))
    allocate(Mass(ngrid,ngrid))
    allocate(invMass(ngrid,ngrid))
#else
    allocate(x1(p1+1))
    allocate(x2(p2+1))
    allocate(wts1(p1+1))
    allocate(wts2(p2+1))
    allocate(G(p2+1,p1+1))
    allocate(D(p2+1,p1+1))
    allocate(D1(p1+1,p1+1))
    allocate(Mass(p1+1,p1+1))
    allocate(invMass(p1+1,p1+1))
    allocate(G1(p1+1,p1+1))
    allocate(W(p2+1,p2+1))
    allocate(V(p1+1,p1+1))
    allocate(iV(p1+1,p1+1))
#endif

    call gll_nodes(p1,x1,wts1)
    call gll_nodes(p2,x2,wts2)
#ifndef FD
    call gll_nodes(p1-1,xd,wtsd)
    call lagint(x1,x1,G1,D1) ! square matrix
    call lagint(x1,x2,G,D) ! rectangle matrix

    call igllm(x1,xd,Gd)
    call igllm(xd,x1,Gu)
#endif

#ifdef FD
    call D1_o4(ngrid,1d0,D,invMass,x1,h)
    !call fd_coeff_upwind_o6(ngrid,1d0,Dp,Dm,invMass,x1,h)
    mesh%Dp = Dp
    mesh%Dm = Dm
    D1 = (Dp+Dm)
    Mass = invMass
    call invert(Mass)
    mesh%fd_D = D
    mesh%fd_Dt = transpose(D)
    !print*,'h=',h
    !print*,sngl(D); stop 2
#endif

    !D = transpose(D)
    !D1 = transpose(D1)

#ifdef FD
    mesh%xnode = x1
#else
    mesh%xnode = xnode_sub(x1,NGLL,Nsub)
#endif
    !print*,mesh%xnode; stop 2

    W(:,:) = 0d0
    do i = 1,p2+1
        W(i,i) = wts2(i)
    end do

    do j = 1,Ngrid
        do i = 1,Ngrid
            !r(i,j) = x1(i)
            !s(i,j) = x1(j)
            r(i,j) = mesh%xnode(i)
            s(i,j) = mesh%xnode(j)
        end do
    end do

#ifndef FD
    Mass = matmul(matmul(transpose(G),W),G)
    invMass = Mass
    call invert(invMass)
#endif

    
    ! Filter
    if (.false.) then
    V(:,:) = 0d0
    do i = 0,order
        do j = 0,order
            V(i+1,j+1) = legendreP(j,x1(i+1)*2d0-1d0)
        end do
    end do
    iV = V
    call invert(iV)
    end if


    !!!if (.false.) then
    !!!    ! modal matrices D,M^-1,V,V^-1
    !!!    invMass(:,:) = 0d0
    !!!    do i = 1,order+1
    !!!        invMass(i,i) = i-0.5d0
    !!!    enddo
    !!!    D(:,:) = 0d0
    !!!    do j = 1,order+1
    !!!        do i = j+1,order+1,2
    !!!            D(i,j) = 2d0
    !!!        end do
    !!!    end do
    !!!end if
#ifndef FD
    eps = 1d-16
    alpha = -dlog(eps)
    sp = 2
    Nc = 3
    Ccut(:,:) = 0d0
    do i = 0,order
        if (i>Nc) then
        Ccut(i+1,i+1) = dexp(-alpha*(dble(i-Nc)/dble(order-Nc))**sp)
        else
        Ccut(i+1,i+1) = 1d0
        end if
    enddo
    !print*,'Cut=',(sngl(Ccut(i,i)),i=1,order+1)

    Filter = matmul(V,matmul(Ccut,iV))
    !Filter = matmul(V,iV)
    !Filter = Ccut
    alpha = 0.1
    Filter = alpha * matmul(Gu , Gd)
    do i = 1,NGLL
        Filter(i,i) = Filter(i,i) + (1.0-alpha) * 1.0
    end do
    !print*,'Filter=',sngl(Filter); stop 2
    !u0 = (/1.0,2.0,1.5,0.3,2.0/)
    !u1 = matmul(Filter,u0)
    !print*,'u0=',u0
    !print*,'u1=',u1
    !stop 2
#endif
    

    mesh%G = G
    mesh%D = D
    mesh%Dt = transpose(D)
    mesh%W = W
    mesh%Mass = Mass
    mesh%invMass= invMass
    mesh%Vdm = V
    mesh%invVdm = iV
    mesh%Filter = Filter

    !print*,invMass
    !print*,sngl(D); stop 2
    !print*,shape(D); stop 2

    nelem = mesh%nelem
    allocate(mesh%vx(Ngrid,Ngrid,nelem))
    allocate(mesh%vy(Ngrid,Ngrid,nelem))
    allocate(mesh%rx(Ngrid,Ngrid,nelem))
    allocate(mesh%ry(Ngrid,Ngrid,nelem))
    allocate(mesh%sx(Ngrid,Ngrid,nelem))
    allocate(mesh%sy(Ngrid,Ngrid,nelem))
    allocate(mesh%detJ(Ngrid,Ngrid,nelem))

    ! surface
    allocate(mesh%nx(Ngrid,nsurface,nelem))
    allocate(mesh%ny(Ngrid,nsurface,nelem))
    allocate(mesh%norm_n(Ngrid,nsurface,nelem))
    ! media
    !allocate(mesh%vp(nelem))
    !allocate(mesh%vs(nelem))
    !allocate(mesh%rho(nelem))
    !allocate(mesh%zp(nelem))
    !allocate(mesh%zs(nelem))
    !allocate(mesh%mu(nelem))
    !allocate(mesh%lam(nelem))

    do i = 1,nelem
        VX = mesh%coord(1,mesh%elem(:,i))
        VY = mesh%coord(2,mesh%elem(:,i))
        !mesh%vx(:,:,i) = (1d0-r)*(1d0-s)*VX(1)+r*(1d0-s)*VX(2)+r*s*(VX(3))+(1d0-r)*s*VX(4)
        !mesh%vy(:,:,i) = (1d0-r)*(1d0-s)*VY(1)+r*(1d0-s)*VY(2)+r*s*(VY(3))+(1d0-r)*s*VY(4)
        mesh%vx(:,:,i) = (1d0-r/steph)*(1d0-s/steph)*VX(1)+r/steph*(1d0-s/steph)*VX(2) &
                       + r/steph*s/steph*(VX(3))+(1d0-r/steph)*s/steph*VX(4)
        mesh%vy(:,:,i) = (1d0-r/steph)*(1d0-s/steph)*VY(1)+r/steph*(1d0-s/steph)*VY(2) &
                       + r/steph*s/steph*(VY(3))+(1d0-r/steph)*s/steph*VY(4)

        ! compute metric analytically ?
        xr = s/steph/steph*vx(3) - s/steph/steph*vx(4) + vx(1)*(s/steph - 1d0)/steph - vx(2)*(s/steph - 1d0)/steph
        xs = r/steph/steph*vx(3) - r/steph/steph*vx(2) + vx(1)*(r/steph - 1d0)/steph - vx(4)*(r/steph - 1d0)/steph
        yr = s/steph/steph*vy(3) - s/steph/steph*vy(4) + vy(1)*(s/steph - 1d0)/steph - vy(2)*(s/steph - 1d0)/steph
        ys = r/steph/steph*vy(3) - r/steph/steph*vy(2) + vy(1)*(r/steph - 1d0)/steph - vy(4)*(r/steph - 1d0)/steph

        detJ = xr*ys-xs*yr

        rx =  ys/detJ
        sx = -yr/detJ
        ry = -xs/detJ
        sy =  xr/detJ

        mesh%rx(:,:,i) = rx
        mesh%sx(:,:,i) = sx
        mesh%ry(:,:,i) = ry
        mesh%sy(:,:,i) = sy
        mesh%detJ(:,:,i) = detJ


#ifdef FD
        !call geometricFactors2d(rx,sx,ry,sy,detJ,mesh%vx(:,:,i),mesh%vy(:,:,i),mesh%D)
#endif
        !print*,shape(rx)
        !print*,shape(mesh%vx(:,:,i))
        !print*,shape(mesh%fd_D)
        mesh%rx(:,:,i) = rx
        mesh%sx(:,:,i) = sx
        mesh%ry(:,:,i) = ry
        mesh%sy(:,:,i) = sy
        mesh%detJ(:,:,i) = detJ

        nx(:,4) = -detJ(1,:)      *rx(1,:)  
        ny(:,4) = -detJ(1,:)      *ry(1,:)  
        nx(:,2) =  detJ(Ngrid,:)  *rx(Ngrid,:)
        ny(:,2) =  detJ(Ngrid,:)  *ry(Ngrid,:)
        nx(:,1) = -detJ(:,1)      *sx(:,1)  
        ny(:,1) = -detJ(:,1)      *sy(:,1)  
        nx(:,3) =  detJ(:,Ngrid)  *sx(:,Ngrid)
        ny(:,3) =  detJ(:,Ngrid)  *sy(:,Ngrid)

        norm_n = dsqrt(nx*nx + ny*ny)

        mesh%nx    (:,:,i) = nx / norm_n
        mesh%ny    (:,:,i) = ny / norm_n
        mesh%norm_n(:,:,i) = norm_n

    end do

    !print*,'r=',mesh%xnode
    !print*,'rx=(',sngl(minval(mesh%rx)),'~',sngl(maxval(mesh%rx)),')'
    !print*,'ry=(',sngl(minval(mesh%ry)),'~',sngl(maxval(mesh%ry)),')'
    !print*,'sx=(',sngl(minval(mesh%sx)),'~',sngl(maxval(mesh%sx)),')'
    !print*,'sy=(',sngl(minval(mesh%sy)),'~',sngl(maxval(mesh%sy)),')'
    !print*,'detJ=(',sngl(minval(mesh%detJ)),'~',sngl(maxval(mesh%detJ)),')'
    !print*,'rx=',sngl(mesh%detJ(:,:,11))
    !print*,'sx=',sngl(mesh%sx(1,:))
    !print*,'ry=',sngl(mesh%ry(1,:))
    !print*,'sy=',sngl(mesh%sy(1,:))
    !print*,'Jac=',sngl(mesh%jac(1,:))

    !call write_mesh(mesh)

    !do i = 1,mesh%nelem
    !    do j = 1,Np
    !        iv(j) = (i-1)*Np + j
    !    enddo

    !    call normals2d(mesh%nx(:,i),mesh%ny(:,i),mesh%sj(:,i),&
    !            mesh%Dr,mesh%Ds,mesh%vx(iv),mesh%vy(iv),mesh%fmask)
    !enddo
    !print*,'nx=',sngl(mesh%nx(:,1))
    !print*,'ny=',sngl(mesh%ny(:,1))
    !print*,'sJ=',sngl(mesh%sJ(:,1))


!    mesh%vp (:) = 6000e-0
!    mesh%vs (:) = 3464.1e-0
!    mesh%rho(:) = 2670e-0
!
!#ifdef TPV10
!    mesh%vp (:) = 5716.0
!    mesh%vs (:) = 3300.0
!    mesh%rho(:) = 2700.0
!#endif
   ! mesh%zp (:) = mesh%rho*mesh%vp
   ! mesh%zs (:) = mesh%rho*mesh%vs
   ! mesh%mu (:) = mesh%rho*(mesh%vs**2)
   ! mesh%lam(:) = mesh%rho*(mesh%vp**2)-2.0*mesh%mu

    call setDT(mesh)

    !print*,'dt =', sngl(mesh%dtfactor)

    !call build_maps2d(mesh)

end subroutine build_geometry

subroutine geometricFactors2d(rx,sx,ry,sy,detJ,x,y,Dr)
    !Input
    real(kind=rkind), intent(in) :: x(:,:),y(:,:),Dr(:,:)
    !Output
    real(kind=rkind), dimension(:,:), intent(out) :: rx,sx,ry,sy,detJ
    !Local
    real(kind=rkind), allocatable, dimension(:,:) :: xr,xs,yr,ys
    !integer :: i

    integer :: s(2),m,n
    s = shape(x)
    m = s(1)
    n = s(2)
    allocate(xr(m,n))
    allocate(xs(m,n))
    allocate(yr(m,n))
    allocate(ys(m,n))

    xr = matmul(Dr,x)
    xs = transpose(matmul(Dr,transpose(x)))
    yr = matmul(Dr,y)
    ys = transpose(matmul(Dr,transpose(y)))

    detJ = -xs*yr + xr*ys

!    do i=1,size(J)
!        if (J(i) .le. 0.0) then
!            call add(errmsg, 2, "Jacobian negative or zero!", myname)
!            if (.level.errmsg == 2) then; call print(errmsg, .true.); stop; endif
!        end if
!    end do

    rx = ys/detJ
    sx = -yr/detJ
    ry = -xs/detJ
    sy = xr/detJ
end subroutine geometricFactors2d

!subroutine normals2d(nx,nz,sj,Dr,Ds,x,z,fmask)
!    !input
!    real(kind=rkind), dimension(:), intent(in) :: x,z
!    real(kind=rkind), dimension(:,:):: Dr,Ds
!    integer, dimension(:,:), intent(in) :: fmask
!    !output
!    real(kind=rkind), dimension(3*NGLL), intent(out) :: nx,nz,sj
!    !local
!    real(kind=rkind), dimension(Np) :: xr,zr,xs,zs
!    real(kind=rkind), dimension(3*NGLL) :: fxr,fzr,fxs,fzs
!    integer, dimension(3*NGLL) :: fmaskv
!    integer, dimension(NGLL) :: fid1,fid2,fid3
!    integer :: i,k,l
!
!    fid1 = 0
!    fid2 = 0
!    fid3 = 0
!    k=1
!    do l=1,3
!        do i=1,NGLL
!            fmaskv(k)=fmask(i,l)
!            if (l==1) then
!                fid1(i)=k
!            else if (l==2) then
!                fid2(i)=k
!            else
!                fid3(i)=k
!            end if
!        k=k+1
!        end do
!    end do
!
!    xr = matmul(Dr,x)
!    zr = matmul(Dr,z)
!    xs = matmul(Ds,x)
!    zs = matmul(Ds,z)
!
!    ! interpolate geometric factors to face nodes
!    fxr = xr(fmaskv)
!    fxs = xs(fmaskv)
!    fzr = zr(fmaskv)
!    fzs = zs(fmaskv)
!
!    ! build normals
!    ! face 1
!    nx(fid1) = fzr(fid1)
!    nz(fid1) = -fxr(fid1)
!    ! face 2
!    nx(fid2) = fzs(fid2)-fzr(fid2)
!    nz(fid2) = -fxs(fid2)+fxr(fid2)
!    ! face 3
!    nx(fid3) = -fzs(fid3)
!    nz(fid3) = fxs(fid3)
!
!    ! normalize
!    sj = sqrt(nx*nx+nz*nz)
!    nx = nx/sJ
!    nz = nz/sJ
!end subroutine normals2d

subroutine setDT(mesh)
    implicit none
    type(meshvar) :: mesh
    real(kind=rkind) :: hmin,hmax,vmax,r,x1,x2,y1,y2,x(nsurface),y(nsurface)
    integer :: i,j,j1

    hmin = 1d30
    hmax = 0

    do i = 1,mesh%nelem
        x = mesh%coord(1,mesh%elem(:,i))
        y = mesh%coord(2,mesh%elem(:,i))
        
        do j = 1,nsurface
            j1 = j+1
            if (j1 > nsurface) then
                j1 = j1-nsurface
            end if
            x1 = x(j); x2 = x(j1)
            y1 = y(j); y2 = y(j1)
            
            r = sqrt((x1-x2)**2+(y1-y2)**2)
            !print*,'i=',i,'face=',j,'r=',r
            if (hmin > r) then
                hmin = r
            end if
            if (hmax < r) then
                hmax = r
            end if
            
        end do
    end do

    vmax = 6000.
#ifdef PSV
    vmax=maxval(mesh%vp)
#endif
#ifdef SH
    vmax=maxval(mesh%vs)
#endif
    !print*,hmin,vmax
    vmax = 1.0
    mesh%dtfactor = hmin/vmax/1.0!dble(order**2)
    mesh%min_grid = hmin
    mesh%max_grid = hmax
    mesh%min_rho = minval(mesh%rho)
    mesh%min_vp  = minval(mesh%vp)
    mesh%min_vs  = minval(mesh%vs)
    mesh%max_rho = maxval(mesh%rho)
    mesh%max_vp  = maxval(mesh%vp)
    mesh%max_vs  = maxval(mesh%vs)

end subroutine setDT

subroutine write_mesh(mesh,n)
    type(meshvar),intent(in) :: mesh
    integer,intent(in) :: n
    integer :: i,j,ie
    character(len=80) :: filename
    write(filename,'(a,i6.6)') 'data/vertices',n
    !open(100,file='data/vertices.txt')
    open(100,file=trim(filename))
    write(100,*) Ngrid,Ngrid,mesh%nelem
    do ie = 1,mesh%nelem
    do j = 1,Ngrid
    do i = 1,Ngrid
        write(100,*) mesh%vx(i,j,ie),mesh%vy(i,j,ie)
    end do
    end do
    end do
    close(100)
end subroutine write_mesh

subroutine invert(V)
    !invert matrix
    real(kind=rkind), dimension(:,:), intent(inout) :: V
    real(kind=rkind), dimension(2*size(v(:,1))) :: work
    integer, dimension(size(v(:,1))) :: ipvt
    integer :: ierr, iwork
    integer :: N
    !character(len=6) :: myname = "invert"
    character(len=180) :: err_str

    N=size(v(:,1))
    iwork=2*N
    !LU trafo
    if (rkind==SIZE_REAL) then
        call sgetrf(N,N,V,N,ipvt,ierr)
    else
        call dgetrf(N,N,V,N,ipvt,ierr)
    endif
    if (ierr/=0) then
        write(err_str,*) "Error LU vdm2d ",ierr
    endif

    ! ivert Pr
    if (rkind==SIZE_REAL) then
        call sgetri(N,V,N,ipvt,work,iwork,ierr)
    else
        call dgetri(N,V,N,ipvt,work,iwork,ierr)
    endif
    if (ierr/=0) then
        write(err_str,*) "Error invert vdm2d ",ierr
    endif
end subroutine invert

end module mod_geometry
