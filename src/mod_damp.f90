module mod_damp
use mod_para
use mod_gll
use mod_mesh
!use mod_jacobi

contains

subroutine init_damp(mesh,rank)
    type(meshvar) :: mesh
    integer :: i,j,ie
    integer :: nelem,rank
    real(kind=rkind) :: x1,x2,y1,y2,L,dpx,dpy
    real(kind=rkind) :: r,r1
    real(kind=rkind) :: axisLimits(4)
    character(len=80) :: filename

    !axisLimits = (/-10d3,20d3,-10d3*dsqrt(3d0),0d0/)
    axisLimits = (/-40d3,50d3,-50d3,0d0/)
    !axisLimits = (/-20d3,30d3,-30d3,0d0/)
    axisLimits = (/-50d3,50d3,-40d3,40d0/)
    axisLimits = (/-1e30,1e30,-1e30,1e30/) ! no damp
    L = 20d3
#ifdef TPV5
    axisLimits = (/-30d3,30d3,-20d3,20d3/)
    axisLimits = (/-50d3,50d3,-40d3,40d3/)
    axisLimits = (/-100d3,100d3,-80d3,80d3/)
    axisLimits = (/-100d3,100d3,-50d3,40000d0/)
    axisLimits = (/-100d3,100d3,-20d3,40000d0/)
    ! acoustic case
    axisLimits = (/-100d3,100d3,-100d3,1d10/)
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

    allocate(mesh%damp(Ngrid,Ngrid,nelem))

    do ie = 1,nelem
        do j = 1,Ngrid
            do i = 1,Ngrid
                if (mesh%vx(i,j,ie) < x1+L) then
                    dpx = dexp(-0.1*(mesh%vx(i,j,ie)-(x1+L))**2/L**2)
                else if (mesh%vx(i,j,ie) >= x1+L .and. mesh%vx(i,j,ie) <= x2-L) then
                    dpx = 1d0
                else
                    dpx = dexp(-0.1*(mesh%vx(i,j,ie)-(x2-L))**2/L**2)
                end if

                if (mesh%vy(i,j,ie) < y1+L) then
                    dpy = dexp(-0.1*(mesh%vy(i,j,ie)-(y1+L))**2/L**2)
                else
                    dpy = 1d0
                end if
#if defined(TPV5) || defined (TPV14)
                if (mesh%vy(i,j,ie) < y1+L) then
                    dpy = dexp(-0.1*(mesh%vy(i,j,ie)-(y1+L))**2/L**2)
                else if (mesh%vy(i,j,ie) >= y1+L .and. mesh%vy(i,j,ie) <= y2-L) then
                    dpy = 1d0
                else
                    dpy = dexp(-0.1*(mesh%vy(i,j,ie)-(y2-L))**2/L**2)
                end if
#endif

                ! circle-type damping distribution
#ifdef TPV14
                r = dsqrt(mesh%vx(i,j,ie)**2+mesh%vy(i,j,ie)**2)
                r1 = 40e3
                L = 10e3
                if (r > r1) then
                    dpx = dexp(-0.1*(r-r1)**2/L**2)
                else
                    dpx = 1d0
                end if
                dpy = 1e30
#endif

                mesh%damp(i,j,ie) = min(dpx,dpy)

            end do
        end do
    end do

    write(filename,'(a,i6.6)') 'data/damp',rank
    open(100,file=trim(filename),access='stream',form='unformatted',status='unknown')
    write(100) sngl(mesh%damp)
    close(100)

end subroutine

end module
