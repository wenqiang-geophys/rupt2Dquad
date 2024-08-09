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
    real(kind=rkind) :: Vp, f0, tmp
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
        mesh%ispml(ie) = 0
        do j = 1,Ngrid
            do i = 1,Ngrid
                x = mesh%vx(i,j,ie)
                y = mesh%vy(i,j,ie)
                ! x direction
                if (x < x1+L) then
                    r = x1+L-x
                    mesh%ispml(ie) = 1
                elseif (x > x2-L) then
                    r = x-(x2-L)
                    mesh%ispml(ie) = 1
                else
                    r = 0
                endif
                pax = amax*(1-r/L)
                if (r==0) pax = 0
                pbx = 1+(bmax-1)*(r/L)**2
                pdx = dmax*(r/L)**2
                ! y direction
                if (y < y1+L) then
                    r = y1+L-y
                    mesh%ispml(ie) = 1
                elseif (y > y2-L) then
                    r = y-(y2-L)
                    mesh%ispml(ie) = 1
                else
                    r = 0
                endif
                pay = amax*(1-r/L)
                if (r==0) pay = 0
                pby = 1+(bmax-1)*(r/L)**2
                pdy = dmax*(r/L)**2

                !if (&
                !        (x<x1+L.and.y<y1+L) .or. &
                !        (x<x1+L.and.y>y2-L) .or. &
                !        (x>x2-L.and.y<y1+L) .or. &
                !        (x>x2-L.and.y>y2-L) ) then
                !tmp = max(mesh%pax(i,j,ie),mesh%pay(i,j,ie))
                !tmp = 0
                !mesh%pax(i,j,ie) = tmp
                !mesh%pay(i,j,ie) = tmp
                !tmp = max(mesh%pbx(i,j,ie),mesh%pby(i,j,ie))
                !tmp = 1
                !mesh%pbx(i,j,ie) = tmp
                !mesh%pby(i,j,ie) = tmp
                !tmp = max(mesh%pdx(i,j,ie),mesh%pdy(i,j,ie))
                !tmp = 0
                !mesh%pdx(i,j,ie) = tmp
                !mesh%pdy(i,j,ie) = tmp
                !end if
                ! turn off pml
                mesh%pax(i,j,ie) = pax
                mesh%pay(i,j,ie) = pay
                mesh%pbx(i,j,ie) = pbx
                mesh%pby(i,j,ie) = pby
                mesh%pdx(i,j,ie) = pdx+0.0*pdy
                mesh%pdy(i,j,ie) = pdy+0.0*pdx
                if (.false.) then
                mesh%pax(i,j,ie) = 0
                mesh%pay(i,j,ie) = 0
                mesh%pbx(i,j,ie) = 1
                mesh%pby(i,j,ie) = 1
                mesh%pdx(i,j,ie) = 0
                mesh%pdy(i,j,ie) = 0
                end if

            end do
        end do
    end do

    write(filename,'(a,i6.6)') 'data/ispml',rank
    open(100,file=trim(filename),access='stream',form='unformatted',status='unknown')
    write(100) real(mesh%ispml)
    close(100)

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
