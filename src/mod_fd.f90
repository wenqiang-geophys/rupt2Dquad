!program main
!
!!use mod_FDo4
!
!implicit none
!
!integer :: N,FDorder,i
!real*8 :: L,h
!real*8,allocatable :: D1(:,:), HI(:,:), x(:)
!real*8,allocatable :: u(:,:), du(:,:)
!
!N = 21
!FDorder = 4
!L = 1d0
!
!allocate(D1(N,N))
!allocate(HI(N,N))
!allocate(u(N,N))
!allocate(du(N,N))
!allocate(x(N))
!
!call D1_o4(N,L,D1,HI,x,h)
!
!do i = 1,N
!    u(i,i) = 1d0
!end do
!
!call fd_d1(D1,FDorder,u,du)
!
!
!
!!print*,'h=',h
!!print*,'x=',x
!!do i = 1,N
!!print*,HI(i,i)
!!end do
!
!do i = 1,FDorder
!    print*,sngl(D1(i,1:FDorder+FDorder/2))
!    print*,sngl(du(i,1:FDorder+FDorder/2))
!end do
!print*,sngl(d1(5,3:3+4))
!print*,sngl(du(5,3:3+4))
!!print*,sngl(d1)
!
!
!contains

module mod_fd

use mod_para
use mod_gll

contains

subroutine mxm(D,u,du)
implicit none
integer :: N,j,i
real*8,intent(in) :: D(:,:),u(:)
real*8,intent(out) :: du(:)

N = size(D,1)
!print*,'N=',N

do i = 1,N
    du(i) = 0
    do j = 1,N
        du(i) = du(i) + D(j,i) * u(j)
    end do
    !du(i) = dot_product(D(i,1:N),u(1:N))
end do

end subroutine

subroutine fd_d1(D,order,u,du)
implicit none
integer :: N,order,nb,m,i
!real*8 :: D(:,:),u(:),du(:)
real*8,intent(in) :: D(:,:),u(:)
real*8,intent(out) :: du(:)
!real*8 :: Db(:,:),D0(:,:)

N = size(D,1)
!print*,'N=',N

m = order/2
nb = order+m

do i = 1,order
    du(i) = dot_product(D(i,1:nb),u(1:nb))
end do
do i = order+1,N-order
    du(i) = dot_product(D(i,(i-m):(i+m)),u((i-m):(i+m)))
end do
do i = N-order+1,N
    du(i) = dot_product(D(i,(N-nb+1):N),u((N-nb+1):N))
end do

end subroutine

subroutine D1_o4(N,L,D1,HI,x,h)
implicit none

integer :: N
real*8 :: L,D1(:,:),HI(:,:),x(:),h
integer :: BP,m,order,i,j
! local vars
real*8,allocatable :: xb(:),P(:),Q(:,:),dd(:)

! BP: Number of boundary points
! m:  Number of nonequidistant spacings
! order: Accuracy of interior stencil
BP = FDorder
m = FDorder/2
order = FDorder

!!!! Non-equidistant grid points !!!!!

!xx = zeros(m+1,1);
allocate(xb(BP+1))
allocate(P(N))
allocate(Q(N,N))
allocate(dd(2*m+1))

if (FDorder == 4) then
xb(0+1) =  0.0000000000000d+00;
xb(1+1) =  6.8764546205559d-01;
xb(2+1) =  1.8022115125776d+00;
xb(3+1) =  2.8022115125776d+00;
xb(4+1) =  3.8022115125776d+00;
else if (FDorder == 8) then
xb(0+1) =  0.0000000000000e+00;
xb(1+1) =  3.8118550247622e-01;
xb(2+1) =  1.1899550868338e+00;
xb(3+1) =  2.2476300175641e+00;
xb(4+1) =  3.3192851303204e+00;
xb(5+1) =  4.3192851303204e+00;
xb(6+1) =  5.3192851303204e+00;
xb(7+1) =  6.3192851303204e+00;
xb(8+1) =  7.3192851303204e+00;
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! Compute h !!!!!!!!!!
h = L/(2d0*xb(m+1) + N-1d0-2d0*m);
!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! Define grid !!!!!!!!
!x = h*[xb; linspace(xb(m+1)+1,L/h-xb(m+1)-1,N-2*(m+1))'; L/h-flip(xb) ];
!x = zeros(N,1);
do i = 0,m
    x(1+i) = xb(1+i);
    x(N-i) = L/h-xb(1+i);
end do

x(m+2:(N-m-1)) = linspace(xb(m+1)+1d0,L/h-xb(m+1)-1d0,N-2*(m+1));
x = x*h;
!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! Norm matrix !!!!!!!!
!P = zeros(BP,1);
!#ok<*NASGU>
if (FDorder == 4) then
P(0+1) =  2.1259737557798d-01;
P(1+1) =  1.0260290400758d+00;
P(2+1) =  1.0775123588954d+00;
P(3+1) =  9.8607273802835d-01;
else if (FDorder == 8) then
P(1+0) =  1.0758368078310e-01;
P(1+1) =  6.1909685107891e-01;
P(1+2) =  9.6971176519117e-01;
P(1+3) =  1.1023441350947e+00;
P(1+4) =  1.0244688965833e+00;
P(1+5) =  9.9533550116831e-01;
P(1+6) =  1.0008236941028e+00;
P(1+7) =  9.9992060631812e-01;
end if

!H = ones(N,1);
!H(:,:) = 1;
!H(1:BP) = P;
!H(N:-1:N-BP+1) = P;

!HI = spdiags(1./H/h,0,N,N);
!H = spdiags(h*H,0,N,N);
do i = 1,N
    HI(i,i) = 1d0
end do
do i = 1,BP
    HI(i,i) = 1d0/P(i)
    HI(N+1-i,N+1-i) = 1d0/P(i)
end do
HI = HI/h

!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! Q matrix !!!!!!!!!!!

! interior stencil
if (FDorder == 4) then
dd = (/1d0/12d0,-2d0/3d0,0d0,2d0/3d0,-1d0/12d0/)
else if (FDorder == 8) then
dd = (/1d0/280d0,-4d0/105d0,1d0/5d0,-4d0/5d0,0d0,4d0/5d0,-1d0/5d0,4d0/105d0,-1d0/280d0/)
end if
!Q = zeros(N,N);
Q(:,:) = 0
do i = 1,N
    do j = -m,m
        if(i+j>0 .and. i+j<=N) then
            Q(i,i+j) = dd(j+m+1)
        end if
    end do
end do

! Boundaries
if (FDorder == 4) then
Q(0+1,0+1) = -5.0000000000000d-01;
Q(0+1,1+1) =  6.5605279837843d-01;
Q(0+1,2+1) = -1.9875859409017d-01;
Q(0+1,3+1) =  4.2705795711740d-02;
Q(0+1,4+1) =  0.0000000000000d+00;
Q(0+1,5+1) =  0.0000000000000d+00;
Q(1+1,0+1) = -6.5605279837843d-01;
Q(1+1,1+1) =  0.0000000000000d+00;
Q(1+1,2+1) =  8.1236966439895d-01;
Q(1+1,3+1) = -1.5631686602052d-01;
Q(1+1,4+1) =  0.0000000000000d+00;
Q(1+1,5+1) =  0.0000000000000d+00;
Q(2+1,0+1) =  1.9875859409017d-01;
Q(2+1,1+1) = -8.1236966439895d-01;
Q(2+1,2+1) =  0.0000000000000d+00;
Q(2+1,3+1) =  6.9694440364211d-01;
Q(2+1,4+1) = -8.3333333333333d-02;
Q(2+1,5+1) =  0.0000000000000d+00;
Q(3+1,0+1) = -4.2705795711740d-02;
Q(3+1,1+1) =  1.5631686602052d-01;
Q(3+1,2+1) = -6.9694440364211d-01;
Q(3+1,3+1) =  0.0000000000000d+00;
Q(3+1,4+1) =  6.6666666666667d-01;
Q(3+1,5+1) = -8.3333333333333d-02;

else if (FDorder == 8) then
Q(1+0,0+1) = -5.0000000000000e-01;
Q(1+0,1+1) =  6.7284756079369e-01;
Q(1+0,2+1) = -2.5969732837062e-01;
Q(1+0,3+1) =  1.3519390385721e-01;
Q(1+0,4+1) = -6.9678474730984e-02;
Q(1+0,5+1) =  2.6434024071371e-02;
Q(1+0,6+1) = -5.5992311465618e-03;
Q(1+0,7+1) =  4.9954552590464e-04;
Q(1+0,8+1) =  0.0000000000000e+00;
Q(1+0,9+1) =  0.0000000000000e+00;
Q(1+0,10+1) =  0.0000000000000e+00;
Q(1+0,11+1) =  0.0000000000000e+00;
Q(1+1,0+1) = -6.7284756079369e-01;
Q(1+1,1+1) =  0.0000000000000e+00;
Q(1+1,2+1) =  9.4074021172233e-01;
Q(1+1,3+1) = -4.0511642426516e-01;
Q(1+1,4+1) =  1.9369192209331e-01;
Q(1+1,5+1) = -6.8638079843479e-02;
Q(1+1,6+1) =  1.3146457241484e-02;
Q(1+1,7+1) = -9.7652615479254e-04;
Q(1+1,8+1) =  0.0000000000000e+00;
Q(1+1,9+1) =  0.0000000000000e+00;
Q(1+1,10+1) =  0.0000000000000e+00;
Q(1+1,11+1) =  0.0000000000000e+00;
Q(1+2,0+1) =  2.5969732837062e-01;
Q(1+2,1+1) = -9.4074021172233e-01;
Q(1+2,2+1) =  0.0000000000000e+00;
Q(1+2,3+1) =  9.4316393361096e-01;
Q(1+2,4+1) = -3.5728039257451e-01;
Q(1+2,5+1) =  1.1266686855013e-01;
Q(1+2,6+1) = -1.8334941452280e-02;
Q(1+2,7+1) =  8.2741521740941e-04;
Q(1+2,8+1) =  0.0000000000000e+00;
Q(1+2,9+1) =  0.0000000000000e+00;
Q(1+2,10+1) =  0.0000000000000e+00;
Q(1+2,11+1) =  0.0000000000000e+00;
Q(1+3,0+1) = -1.3519390385721e-01;
Q(1+3,1+1) =  4.0511642426516e-01;
Q(1+3,2+1) = -9.4316393361096e-01;
Q(1+3,3+1) =  0.0000000000000e+00;
Q(1+3,4+1) =  8.7694387866575e-01;
Q(1+3,5+1) = -2.4698058719506e-01;
Q(1+3,6+1) =  4.7291642094198e-02;
Q(1+3,7+1) = -4.0135203618880e-03;
Q(1+3,8+1) =  0.0000000000000e+00;
Q(1+3,9+1) =  0.0000000000000e+00;
Q(1+3,10+1) =  0.0000000000000e+00;
Q(1+3,11+1) =  0.0000000000000e+00;
Q(1+4,0+1) =  6.9678474730984e-02;
Q(1+4,1+1) = -1.9369192209331e-01;
Q(1+4,2+1) =  3.5728039257451e-01;
Q(1+4,3+1) = -8.7694387866575e-01;
Q(1+4,4+1) =  0.0000000000000e+00;
Q(1+4,5+1) =  8.1123946853807e-01;
Q(1+4,6+1) = -2.0267150541446e-01;
Q(1+4,7+1) =  3.8680398901392e-02;
Q(1+4,8+1) = -3.5714285714286e-03;
Q(1+4,9+1) =  0.0000000000000e+00;
Q(1+4,10+1) =  0.0000000000000e+00;
Q(1+4,11+1) =  0.0000000000000e+00;
Q(1+5,0+1) = -2.6434024071371e-02;
Q(1+5,1+1) =  6.8638079843479e-02;
Q(1+5,2+1) = -1.1266686855013e-01;
Q(1+5,3+1) =  2.4698058719506e-01;
Q(1+5,4+1) = -8.1123946853807e-01;
Q(1+5,5+1) =  0.0000000000000e+00;
Q(1+5,6+1) =  8.0108544742793e-01;
Q(1+5,7+1) = -2.0088756283071e-01;
Q(1+5,8+1) =  3.8095238095238e-02;
Q(1+5,9+1) = -3.5714285714286e-03;
Q(1+5,10+1) =  0.0000000000000e+00;
Q(1+5,11+1) =  0.0000000000000e+00;
Q(1+6,0+1) =  5.5992311465618e-03;
Q(1+6,1+1) = -1.3146457241484e-02;
Q(1+6,2+1) =  1.8334941452280e-02;
Q(1+6,3+1) = -4.7291642094198e-02;
Q(1+6,4+1) =  2.0267150541446e-01;
Q(1+6,5+1) = -8.0108544742793e-01;
Q(1+6,6+1) =  0.0000000000000e+00;
Q(1+6,7+1) =  8.0039405922650e-01;
Q(1+6,8+1) = -2.0000000000000e-01;
Q(1+6,9+1) =  3.8095238095238e-02;
Q(1+6,10+1) = -3.5714285714286e-03;
Q(1+6,11+1) =  0.0000000000000e+00;
Q(1+7,0+1) = -4.9954552590464e-04;
Q(1+7,1+1) =  9.7652615479254e-04;
Q(1+7,2+1) = -8.2741521740941e-04;
Q(1+7,3+1) =  4.0135203618880e-03;
Q(1+7,4+1) = -3.8680398901392e-02;
Q(1+7,5+1) =  2.0088756283071e-01;
Q(1+7,6+1) = -8.0039405922650e-01;
Q(1+7,7+1) =  0.0000000000000e+00;
Q(1+7,8+1) =  8.0000000000000e-01;
Q(1+7,9+1) = -2.0000000000000e-01;
Q(1+7,10+1) =  3.8095238095238e-02;
Q(1+7,11+1) = -3.5714285714286e-03;
end if

do i = 1,BP
    do j = 1,BP
        Q(N+1-i,N+1-j) = -Q(i,j)
    end do
end do

D1 = matmul(HI,Q)

end subroutine

subroutine DpO6(D,u,du)
implicit none
integer :: N,i,j
real*8 :: D(:,:),u(:),du(:)

N = size(D,1)

du(1) = dot_product(D(1,1:6),u(1:6))
du(2) = dot_product(D(2,1:6),u(1:6))
du(3) = dot_product(D(3,1:7),u(1:7))
du(4) = dot_product(D(4,1:8),u(1:8))
du(5) = dot_product(D(5,1:9),u(1:9))
du(6) = dot_product(D(6,1:10),u(1:10))

do i = 7,N-6
    du(i) = dot_product(D(i,(i-2):(i+4)),u((i-2):(i+4)))
end do

du(N  ) = dot_product(D(N  ,(N-5):N),u((N-5):N))
du(N-1) = dot_product(D(N-1,(N-5):N),u((N-5):N))
du(N-2) = dot_product(D(N-2,(N-5):N),u((N-5):N))
du(N-3) = dot_product(D(N-3,(N-5):N),u((N-5):N))
du(N-4) = dot_product(D(N-4,(N-6):N),u((N-6):N))
du(N-5) = dot_product(D(N-5,(N-7):N),u((N-7):N))

end subroutine

subroutine DmO6(D,u,du)
implicit none
integer :: N,i,j
real*8 :: D(:,:),u(:),du(:)

N = size(D,1)

du(1) = dot_product(D(1,1:6),u(1:6))
du(2) = dot_product(D(2,1:6),u(1:6))
du(3) = dot_product(D(3,1:6),u(1:6))
du(4) = dot_product(D(4,1:6),u(1:6))
du(5) = dot_product(D(5,1:7),u(1:7))
du(6) = dot_product(D(6,1:8),u(1:8))

do i = 7,N-6
    du(i) = dot_product(D(i,(i-4):(i+2)),u((i-4):(i+2)))
end do

du(N  ) = dot_product(D(N  ,(N-5):N),u((N-5):N))
du(N-1) = dot_product(D(N-1,(N-5):N),u((N-5):N))
du(N-2) = dot_product(D(N-2,(N-6):N),u((N-6):N))
du(N-3) = dot_product(D(N-3,(N-7):N),u((N-7):N))
du(N-4) = dot_product(D(N-4,(N-8):N),u((N-8):N))
du(N-5) = dot_product(D(N-5,(N-9):N),u((N-9):N))

end subroutine

subroutine fd_coeff_upwind_o6(N,L,Dp,Dm,HI,x,h)
implicit none

integer :: N
real*8 :: L,Dp(:,:),Dm(:,:),HI(:,:),x(:),h
integer :: BP,m,order,i,j
! local vars
real*8,allocatable :: P(:),Q(:,:),Qp(:,:),Qm(:,:),Qtmp(:,:),q_stencil(:)
integer,allocatable :: q_diags(:)

! BP: Number of boundary points
! m:  Number of nonequidistant spacings
! order: Accuracy of interior stencil
BP = 6;
m = 3;
order = 6;

!!!! Non-equidistant grid points !!!!!

!xx = zeros(m+1,1);
allocate(P(N))
allocate(Q(N,N))
allocate(Qtmp(N,N))
allocate(q_diags(2*m+1))
allocate(q_stencil(2*m+1))

h = 1d0/dble(N-1)

x = linspace(0d0,1d0,N)

!x = x*h;
!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! Norm matrix !!!!!!!!
P(1:6)=(/0.13613e5/0.43200e5,0.12049e5/0.8640e4,0.535e3/0.864e3,0.1079e4/0.864e3,0.7841e4/0.8640e4,0.43837e5/0.43200e5/)

do i = 1,N
    HI(i,i) = 1d0
end do
do i = 1,BP
    HI(i,i) = 1d0/P(i)
    HI(N+1-i,N+1-i) = 1d0/P(i)
end do
HI = HI/h

!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! Q matrix !!!!!!!!!!!

! interior stencil
q_diags   = (/-2,-1,0,1,2,3,4/)
q_stencil = (/1d0/30d0,-2d0/5d0,-7d0/12d0,4d0/3d0,-1d0/2d0,2d0/15d0,-1d0/60d0/)
!Q = zeros(N,N);
Q(:,:) = 0d0
do i = 1,N
    !do j = -m,m
    do j = q_diags(1),q_diags(2*m+1)
        if(i+j>0 .and. i+j<=N) then
            !Q(i,i+j) = q_stencil(j+m+1)
            Q(i,i+j) = q_stencil(j-q_diags(1)+1)
        end if
    end do
end do

! Boundarids
Q(1,1) = -2.0592440631605124e-03;
Q(1,2) =  6.5975852420554337e-01;
Q(1,3) = -8.6848291993037424e-02;
Q(1,4) = -9.5434661955876049e-02;
Q(1,5) =  7.5605826437527921e-03;
Q(1,6) =  1.7023091162777847e-02;
Q(2,1) = -6.4266291311515422e-01;
Q(2,2) = -4.1211301753077210e-02;
Q(2,3) =  5.4838020236138163e-01;
Q(2,4) =  1.4779568643771213e-01;
Q(2,5) =  4.2867607443325462e-02;
Q(2,6) = -5.5169281374187813e-02;
Q(3,1) =  5.9650728894691037e-02;
Q(3,2) = -3.8964589924065557e-01;
Q(3,3) = -1.8981412408305359e-01;
Q(3,4) =  7.7944473800544301e-01;
Q(3,5) = -3.6499678090218779e-01;
Q(3,6) =  1.2202800399242958e-01;
Q(4,1) =  1.1563856597179052e-01;
Q(4,2) = -3.0041908566685549e-01;
Q(4,3) = -3.0432351458134749e-01;
Q(4,4) = -3.9749782419495211e-01;
Q(4,5) =  1.1826071123498234e+00;
Q(4,6) = -4.1267192054512553e-01;
Q(5,1) = -1.4165705110494057e-02;
Q(5,2) =  3.0388640225454847e-02;
Q(5,3) =  3.8023014894479214e-02;
Q(5,4) = -4.0850540900530019e-01;
Q(5,5) = -5.3864074350366775e-01;
Q(5,6) =  1.2762335358328614e+00;
Q(6,1) = -1.6401432577672785e-02;
Q(6,2) =  4.1129122229590025e-02;
Q(6,3) = -5.4172865984223684e-03;
Q(6,4) = -2.5802529287026674e-02;
Q(6,5) = -3.6273111136437941e-01;
Q(6,6) = -5.8077676240208875e-01;
do i = 1,BP
    do j = 1,BP
        !Q(N+1-i,N+1-j) = -Q(i,j)
        Q(N+1-i,N+1-j) = Q(j,i)
    end do
end do

Qp = Q
Qm = -transpose(Qp)

Qtmp = Qp
Qtmp(1,1) = Qtmp(1,1)-0.5d0
Qtmp(N,N) = Qtmp(N,N)+0.5d0
Dp = matmul(HI,Qtmp)

Qtmp = Qm
Qtmp(1,1) = Qtmp(1,1)-0.5d0
Qtmp(N,N) = Qtmp(N,N)+0.5d0
Dm = matmul(HI,Qtmp)
!do i = 1,BP
!    !print*,sngl(Q(i,1:BP))
!    print*,sngl(Q(N-BP+i,N-BP+1:N))
!end do

!do i = 1,N
!    print*,sngl(Dm(i,:))
!end do

end subroutine


!end program
end module
