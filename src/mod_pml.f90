module mod_pml
use mod_para
use mod_gll
use mod_mesh
!use mod_jacobi

implicit none

contains

subroutine init_pml(mesh,rank)
    type(meshvar) :: mesh
    integer :: i,j,ie
    integer :: nelem,rank
    real(kind=rkind) :: x1,x2,y1,y2,L,dpx,dpy
    real(kind=rkind) :: pax,pay,pbx,pby,pdx,pdy
    real(kind=rkind) :: x,y
    real(kind=rkind) :: amax,bmax,dmax
    real(kind=rkind) :: Vp, f0
    real(kind=rkind) :: r,r1
    real(kind=rkind) :: axisLimits(4)
    character(len=80) :: filename

    !axisLimits = (/-10d3,20d3,-10d3*dsqrt(3d0),0d0/)
    axisLimits = (/-40d3,50d3,-50d3,0d0/)
    !axisLimits = (/-20d3,30d3,-30d3,0d0/)
    axisLimits = (/-50d3,50d3,-40d3,40d0/)
    axisLimits = (/-1e30,1e30,-1e30,1e30/) ! no damp
    L = 20d3
    bmax = 3
    Vp = 8000
    f0 = 1
    dmax = -3*Vp/(2*L)*log(1e-4)
    amax = pi * f0
    print*,'amax=',amax
    print*,'dmax=',dmax

#ifdef TPV5
    axisLimits = (/-30d3,30d3,-20d3,20d3/)
    axisLimits = (/-50d3,50d3,-40d3,40d3/)
    axisLimits = (/-100d3,100d3,-80d3,80d3/)
    axisLimits = (/-100d3,100d3,-50d3,40000d0/)
    axisLimits = (/-100d3,100d3,-20d3,40000d0/)
    ! acoustic case
    axisLimits = (/-100d3,100d3,-100d3,1d10/)
    ! fullspace case
    axisLimits = (/-100d3,100d3,-50d3,1d10/)
    axisLimits = (/-50d3,50d3,-100d3,0d0/)
    !axisLimits = (/-1e30,1e30,-1e30,1e30/) ! no damp
#endif
#ifdef TPV14
    axisLimits = (/-30d3,30d3,-20d3,20d3/)
#endif

    x1 = axisLimits(1)
    x2 = axisLimits(2)
    y1 = axisLimits(3)
    y2 = axisLimits(4)

    nelem = mesh%nelem

    allocate(mesh%pax(Ngrid,Ngrid,nelem))
    allocate(mesh%pbx(Ngrid,Ngrid,nelem))
    allocate(mesh%pdx(Ngrid,Ngrid,nelem))
    allocate(mesh%pay(Ngrid,Ngrid,nelem))
    allocate(mesh%pby(Ngrid,Ngrid,nelem))
    allocate(mesh%pdy(Ngrid,Ngrid,nelem))

    do ie = 1,nelem
        do j = 1,Ngrid
            do i = 1,Ngrid
                x = mesh%vx(i,j,ie)
                y = mesh%vy(i,j,ie)
                ! x direction
                if (x < x1+L) then
                    r = x1+L-x
                elseif (x > x2-L) then
                    r = x-(x2-L)
                else
                    r = 0
                endif
                mesh%pax(i,j,ie) = amax*(1-r/L)
                mesh%pbx(i,j,ie) = 1+(bmax-1)*(r/L)**2
                mesh%pdx(i,j,ie) = dmax*(r/L)**2
                ! y direction
                if (y < y1+L) then
                    r = y1+L-y
                elseif (y > y2-L) then
                    r = y-(y2-L)
                else
                    r = 0
                endif
                mesh%pay(i,j,ie) = amax*(1-r/L)
                mesh%pby(i,j,ie) = 1+(bmax-1)*(r/L)**2
                mesh%pdy(i,j,ie) = dmax*(r/L)**2
            end do
        end do
    end do

    write(filename,'(a,i6.6)') 'data/pax',rank
    open(100,file=trim(filename),access='stream',form='unformatted',status='unknown')
    write(100) sngl(mesh%pax)
    close(100)
    write(filename,'(a,i6.6)') 'data/pbx',rank
    open(100,file=trim(filename),access='stream',form='unformatted',status='unknown')
    write(100) sngl(mesh%pbx)
    close(100)
    write(filename,'(a,i6.6)') 'data/pdx',rank
    open(100,file=trim(filename),access='stream',form='unformatted',status='unknown')
    write(100) sngl(mesh%pdx)
    close(100)
    write(filename,'(a,i6.6)') 'data/pay',rank
    open(100,file=trim(filename),access='stream',form='unformatted',status='unknown')
    write(100) sngl(mesh%pay)
    close(100)
    write(filename,'(a,i6.6)') 'data/pby',rank
    open(100,file=trim(filename),access='stream',form='unformatted',status='unknown')
    write(100) sngl(mesh%pby)
    close(100)
    write(filename,'(a,i6.6)') 'data/pdy',rank
    open(100,file=trim(filename),access='stream',form='unformatted',status='unknown')
    write(100) sngl(mesh%pdy)
    close(100)

end subroutine

end module
