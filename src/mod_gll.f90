module mod_gll

!real*8,parameter :: PI = 4d0*atan(1d0)
use mod_para,only : PI, steph

contains

subroutine gll_nodes(N,xi,wts)
implicit none
integer,intent(in) :: N
real*8,intent(out) :: xi(:),wts(:)
integer :: N1, k
real*8,allocatable :: xold(:),x(:),P(:,:)

N1 = N+1

allocate(xold(N1))
allocate(x(N1))
allocate(P(N1,N1))

xold(:) = 2d0
do k = 0,N
    x(k+1) = dcos(PI*k/N)
end do

do while (maxval(abs(x-xold)) > 1d-8)
    xold = x
    P(:,1) = 1d0
    P(:,2) = x

    do k = 2,N
        P(:,k+1) = ((2d0*k-1d0)*x*P(:,k)-(k-1d0)*P(:,k-1))/dble(k)
    end do

    x = xold - (x*P(:,N1)-P(:,N))/(N1*P(:,N1))
end do

wts = 2d0/(N*N1*P(:,N1)**2)
! convert [-1 1] to [0 1]
xi = 0.5d0*(x(N1:1:-1)+1d0)
wts = 0.5d0*wts(N1:1:-1)

! covert to [0,steph]
xi = xi * steph
wts = wts * steph

end subroutine

subroutine gl_nodes(N,xi,wts)
implicit none
integer,intent(in) :: N
real*8,intent(out) :: xi(:),wts(:)

if (N == 1) then
    xi(1) = -1d0/3d0*dsqrt(3d0)
    xi(2) =  1d0/3d0*dsqrt(3d0)

    wts(1) = 1d0
    wts(2) = 1d0
else if (N == 2) then
    xi(1) = -1d0/5d0*dsqrt(15d0)
    xi(2) = 0d0
    xi(3) =  1d0/5d0*dsqrt(15d0)

    wts(1) = 5d0/9d0
    wts(2) = 8d0/9d0
    wts(3) = 5d0/9d0
else if (N == 3) then
    xi(1) = -1d0/(35d0)*dsqrt(525d0+70d0*dsqrt(30d0))
    xi(2) = -1d0/(35d0)*dsqrt(525d0-70d0*dsqrt(30d0))
    xi(3) =  1d0/(35d0)*dsqrt(525d0-70d0*dsqrt(30d0))
    xi(4) =  1d0/(35d0)*dsqrt(525d0+70d0*dsqrt(30d0))

    wts(1) = 1d0/(36d0)*(18d0-dsqrt(30d0))
    wts(2) = 1d0/(36d0)*(18d0+dsqrt(30d0))
    wts(3) = 1d0/(36d0)*(18d0+dsqrt(30d0))
    wts(4) = 1d0/(36d0)*(18d0-dsqrt(30d0))
else if (N == 4) then
    xi(1) = -1.0/(21d0)*dsqrt(245d0+14d0*dsqrt(70d0))
    xi(2) = -1.0/(21d0)*dsqrt(245d0-14d0*dsqrt(70d0))
    xi(3) =  0.0d0
    xi(4) =  1.0/(21d0)*dsqrt(245d0-14d0*dsqrt(70d0))
    xi(5) =  1.0/(21d0)*dsqrt(245d0+14d0*dsqrt(70d0))

    wts(1) = 1d0/(900d0)*(322d0-13d0*dsqrt(70d0))
    wts(2) = 1d0/(900d0)*(322d0+13d0*dsqrt(70d0))
    wts(3) = (128d0)/(225d0)
    wts(4) = 1d0/(900d0)*(322d0+13d0*dsqrt(70d0))
    wts(5) = 1d0/(900d0)*(322d0-13d0*dsqrt(70d0))
else if (N == 5) then
    xi(1) = -0.93246951420315
    xi(2) = -0.66120938646626
    xi(3) = -0.23861918608320
    xi(4) =  0.23861918608320
    xi(5) =  0.66120938646626
    xi(6) =  0.93246951420315

    wts(1) = 0.17132449237917
    wts(2) = 0.36076157304814
    wts(3) = 0.46791393457269
    wts(4) = 0.46791393457269
    wts(5) = 0.36076157304814
    wts(6) = 0.17132449237917
else if (N == 6) then
    xi(1) = -0.94910791234276
    xi(2) = -0.74153118559939    
    xi(3) = -0.40584515137740
    xi(4) =  0.00000000000000
    xi(5) =  0.40584515137740
    xi(6) =  0.74153118559939
    xi(7) =  0.94910791234276
    
    wts(1) = 0.12948496616887
    wts(2) = 0.27970539148928
    wts(3) = 0.38183005050512
    wts(4) = 0.41795918367347
    wts(5) = 0.38183005050512
    wts(6) = 0.27970539148928
    wts(7) = 0.12948496616887
else if (N == 7) then
    xi(1)  = -0.960289856497536231684
    xi(2)  = -0.7966664774136267395916
    xi(3)  = -0.5255324099163289858177
    xi(4)  = -0.1834346424956498049395
    xi(5)  =  0.1834346424956498049395
    xi(6)  =  0.525532409916328985818
    xi(7)  =  0.796666477413626739592
    xi(8)  =  0.9602898564975362316836
       
    wts(1) = 0.1012285362903762591525
    wts(2) = 0.222381034453374470544
    wts(3) = 0.313706645877887287338
    wts(4) = 0.3626837833783619829652
    wts(5) = 0.3626837833783619829652
    wts(6) = 0.31370664587788728734
    wts(7) = 0.222381034453374470544
    wts(8) = 0.1012285362903762591525
else if (N == 8) then
    xi(1)   = -0.968160239507626089836
    xi(2)   = -0.8360311073266357942994
    xi(3)   = -0.6133714327005903973087
    xi(4)   = -0.3242534234038089290385 
    xi(5)   =  0.000000000000000000000
    xi(6)   =  0.324253423403808929039
    xi(7)   =  0.6133714327005903973087
    xi(8)   =  0.8360311073266357942994
    xi(9)   =  0.9681602395076260898356
       
    wts(1)  = 0.0812743883615744119723
    wts(2)  = 0.180648160694857404058
    wts(3)  = 0.2606106964029354623187
    wts(4)  = 0.312347077040002840069
    wts(5)  = 0.330239355001259763165
    wts(6)  = 0.312347077040002840069
    wts(7)  = 0.2606106964029354623187
    wts(8)  = 0.180648160694857404058
    wts(9)  = 0.081274388361574411972
else if (N == 9) then
    xi(1)   = -0.973906528517171720078
    xi(2)   = -0.8650633666889845107321
    xi(3)   = -0.6794095682990244062343
    xi(4)   = -0.4333953941292471907993
    xi(5)   = -0.148874338981631210885
    xi(6)   =  0.1488743389816312108848
    xi(7)   =  0.433395394129247190799
    xi(8)   =  0.6794095682990244062343
    xi(9)   =  0.8650633666889845107321
    xi(10)  =  0.973906528517171720078
       
    wts(1)  = 0.0666713443086881375936
    wts(2)  = 0.149451349150580593146
    wts(3)  = 0.2190863625159820439955
    wts(4)  = 0.2692667193099963550912
    wts(5)  = 0.295524224714752870174
    wts(6)  = 0.295524224714752870174
    wts(7)  = 0.2692667193099963550913
    wts(8)  = 0.219086362515982043995
    wts(9)  = 0.1494513491505805931458
    wts(10) = 0.0666713443086881375936
else if (N == 10) then
    xi(1)   = -0.978228658146056992804
    xi(2)   = -0.8870625997680952990752
    xi(3)   = -0.730152005574049324093
    xi(4)   = -0.519096129206811815926
    xi(5)   = -0.2695431559523449723315
    xi(6)   =  0.0000000000000000000000
    xi(7)   =  0.2695431559523449723315
    xi(8)   =  0.5190961292068118159257
    xi(9)   =  0.7301520055740493240934
    xi(10)  =  0.887062599768095299075
    xi(11)  =  0.9782286581460569928039

    wts(1)  = 0.0556685671161736664828
    wts(2)  = 0.1255803694649046246347
    wts(3)  = 0.1862902109277342514261
    wts(4)  = 0.2331937645919904799185
    wts(5)  = 0.262804544510246662181
    wts(6)  = 0.2729250867779006307145
    wts(7)  = 0.262804544510246662181
    wts(8)  = 0.233193764591990479919
    wts(9)  = 0.186290210927734251425
    wts(10) = 0.125580369464904624635
    wts(11) = 0.055668567116173666483
else if (N == 11) then
    xi(1)   = -0.9815606342467192506905
    xi(2)   = -0.904117256370474856678
    xi(3)   = -0.769902674194304687037
    xi(4)   = -0.587317954286617447297
    xi(5)   = -0.3678314989981801937527
    xi(6)   = -0.1252334085114689154724
    xi(7)   =  0.1252334085114689154724
    xi(8)   =  0.3678314989981801937527
    xi(9)   =  0.587317954286617447297
    xi(10)  =  0.769902674194304687037
    xi(11)  =  0.9041172563704748566785
    xi(12)  =  0.981560634246719250691
            
    wts(1)  = 0.047175336386511827195
    wts(2)  = 0.10693932599531843096
    wts(3)  = 0.160078328543346226335
    wts(4)  = 0.2031674267230659217491
    wts(5)  = 0.2334925365383548087608
    wts(6)  = 0.249147045813402785001
    wts(7)  = 0.249147045813402785001
    wts(8)  = 0.2334925365383548087609
    wts(9)  = 0.203167426723065921749
    wts(10) = 0.160078328543346226335
    wts(11) = 0.1069393259953184309603
    wts(12) = 0.0471753363865118271946
else if (N == 12) then
    xi(1)   =  -0.984183054718588149473
    xi(2)   =  -0.917598399222977965207
    xi(3)   =  -0.8015780907333099127942
    xi(4)   =  -0.642349339440340220644
    xi(5)   =  -0.4484927510364468528779
    xi(6)   =  -0.2304583159551347940655
    xi(7)   =   0.000000000000000000000
    xi(8)   =   0.230458315955134794066
    xi(9)   =   0.448492751036446852878
    xi(10)  =   0.642349339440340220644
    xi(11)  =   0.8015780907333099127942
    xi(12)  =   0.9175983992229779652065
    xi(13)  =   0.9841830547185881494728
             
    wts(1)  = 0.04048400476531587952
    wts(2)  = 0.092121499837728447914
    wts(3)  = 0.138873510219787238464
    wts(4)  = 0.17814598076194573828
    wts(5)  = 0.2078160475368885023125
    wts(6)  = 0.226283180262897238412
    wts(7)  = 0.2325515532308739101946
    wts(8)  = 0.2262831802628972384122
    wts(9)  = 0.207816047536888502312
    wts(10) = 0.17814598076194573828
    wts(11) = 0.1388735102197872384636
    wts(12) = 0.092121499837728447914
    wts(13) = 0.04048400476531587952
end if

    xi = 0.5d0*(xi+1d0)
    wts = 0.5d0*wts

end subroutine

subroutine igllm(x1,x2,G)
implicit none
real*8,intent(in) :: x1(:),x2(:)
real*8,intent(out) :: G(:,:)
real*8,allocatable :: w(:)
integer :: i,j,k,p1,p2
real*8 :: ell

p1 = size(x1)-1
p2 = size(x2)-1

allocate(w(p1+1))

w(:) = 1d0

do j = 1,p1+1
    do k = 1,p1+1
        if (j .ne. k) then
            w(j) = w(j) * (x1(j) - x1(k))
        end if
    end do
end do

w = 1d0/w

G(:,:) = 0d0

do i = 1,p2+1
    ell = 1d0
    do j = 1,p1+1
        ell = ell*(x2(i) - x1(j))
    end do
    do j = 1,p1+1
        !if (x2(i) == x1(j)) then
        if (abs(x2(i)-x1(j))<1e-20) then
            G(i,j) = 1d0;
        else
            G(i,j) = ell*w(j)/(x2(i) - x1(j))
        end if
    end do
end do

end subroutine

subroutine lagint(x1,x2,G,GD)
implicit none
real*8,intent(in) :: x1(:),x2(:)
real*8,intent(out) :: G(:,:),GD(:,:)
real*8,allocatable :: D(:,:),w(:)
integer :: i,j,k,p1,p2
real*8 :: ell

p1 = size(x1)-1
p2 = size(x2)-1

allocate(D(p1+1,p1+1))
!allocate(G(p2+1,p1+1))
allocate(w(p1+1))

w(:) = 1d0

do j = 1,p1+1
    do k = 1,p1+1
        if (j .ne. k) then
            w(j) = w(j) * (x1(j) - x1(k))
        end if
    end do
end do

w = 1d0/w

do i = 1,p1+1
    do j = 1,p1+1
        if (i .ne. j) then
            D(i,j) = w(j)/w(i) / (x1(i) - x1(j))
        end if
    end do
    D(i,i) = 0d0;
    do j = 1,p1+1
        if (i .ne. j) then
            D(i,i) = D(i,i) - D(i,j)
        end if
    end do
end do

G(:,:) = 0d0

do i = 1,p2+1
    ell = 1d0
    do j = 1,p1+1
        ell = ell*(x2(i) - x1(j))
    end do
    do j = 1,p1+1
        !if (x2(i) == x1(j)) then
        if (abs(x2(i)-x1(j))<1e-20) then
            G(i,j) = 1d0;
        else
            G(i,j) = ell*w(j)/(x2(i) - x1(j))
        end if
    end do
end do

GD = matmul(G,D)

end subroutine


function lagrange_basis(Np,i,x,xi)
    real*8,intent(in) :: x,xi(:)
    integer,intent(in) :: Np,i
    real*8 :: lagrange_basis
    real*8 :: h,num,den
    integer :: j
    h = 1d0
    do j = 1,Np
        if (j .ne. i) then
            num = x - xi(j)
            den = xi(i)-xi(j)
            h = h * num/den
        end if
    end do
    lagrange_basis = h
end function

function lagrange_interpol(Np,x,xi,f)
    real*8,intent(in) :: x,xi(:),f(:)
    integer,intent(in) :: Np
    real*8 :: lagrange_interpol
    real*8 :: a_i, g
    integer :: i
    g = 0d0
    do i = 1,Np
        a_i = lagrange_basis(Np,i,x,xi)
        g = g + a_i * f(i)
    end do
    lagrange_interpol = g
end function

function legendreP(i,x) result(y)
    implicit none
    integer,intent(in) :: i
    real*8,intent(in) :: x
    real*8 :: y

    if (i == 0) then
        y = 1d0
    elseif (i == 1) then
        y = x
    elseif (i == 2) then
        y = (3d0*x**2)/2d0 - 1d0/2d0
    elseif (i == 3) then
        y = (5d0*x**3)/2d0 - (3d0*x)/2d0
    elseif (i == 4) then
        y = (35d0*x**4)/8d0 - (15d0*x**2)/4d0 + 3d0/8d0
    elseif (i == 5) then
        y = (63d0*x**5)/8d0 - (35d0*x**3)/4d0 + (15d0*x)/8d0
    elseif (i == 6) then
        y = (231d0*x**6)/16d0 - (315d0*x**4)/16d0 + (105d0*x**2)/16d0 - 5d0/16d0
    else 
        y = 1e38
        print*,'order must less than 6!'
    end if
end function

function linspace(a,b,n)
implicit none
real*8 :: a,b
integer :: n,i
real*8,allocatable :: linspace(:)
allocate( linspace(n))
do i = 1,n
    linspace(i) = a+(b-a)/dble(n-1)*(i-1)
end do
end function

function xnode_sub(xnode,n,nsub)
implicit none
integer,intent(in) :: n,nsub
real*8,intent(in) :: xnode(:)
real*8,allocatable :: xnode_sub(:),xeq(:)
integer :: i,j,Ngrid

Ngrid = n*nsub+1
allocate( xeq(nsub+1) )
allocate( xnode_sub(Ngrid) )

xeq = linspace(0d0,1d0,nsub+1)

do i = 1,nsub
    do j = 1,n
        xnode_sub(j+(i-1)*(n-1))=xeq(i)+xnode(j)*(xeq(i+1)-xeq(i))
    end do
end do


end function

end module
