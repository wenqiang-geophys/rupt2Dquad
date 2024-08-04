module mod_fault

use mod_para
use mod_mesh
use mod_wave
use mod_types

implicit none

contains

subroutine fault_init(mesh)
    type(meshvar) :: mesh
    real(rkind) :: x,y, asp_size
    integer :: element,surface,i,j,k,is
    integer :: nelem_fault
    real(rkind) :: sxx,syy,sxy,vec_n(2),vec_m(2),Tx,Ty,Tn,Tm,tmp
    !real(rkind) :: xtmp,ytmp
    !real(rkind) :: vx1,vx2,vy1,vy2
    real(rkind) :: damp,r1,r2
    real(rkind) :: xface(Nfp),yface(Nfp)
    integer :: ie,ief
    logical :: isfault

    !mesh%nelem_fault = 0
    !do ie = 1,mesh%nelem
    !do is = 1,Nfaces
    !if (mesh%bctype(is,ie) >= BC_FAULT) then
    !    mesh%nelem_fault = mesh%nelem_fault + 1
    !end if
    !end do
    !end do

    mesh%nelem_fault = 0
    do ie = 1,mesh%nelem
      isfault=.false.
      do is = 1,Nfaces
        if (mesh%bctype(is,ie) >= BC_FAULT) then
          isfault=.true.
        end if
      end do
      if (isfault) then
          mesh%nelem_fault = mesh%nelem_fault + 1
      end if
    end do

    mesh%nface_fault = 0
    do ie = 1,mesh%nelem
      do is = 1,Nfaces
        if (mesh%bctype(is,ie) >= BC_FAULT) then
          mesh%nface_fault = mesh%nface_fault + 1
        end if
      end do
    end do

    ! count free faces for ground surface output
    mesh%nface_free = 0
    do ie = 1,mesh%nelem
      do is = 1,Nfaces
        if (mesh%bctype(is,ie) == BC_FREE) then
          mesh%nface_free = mesh%nface_free + 1
        end if
      end do
    end do

    !print*,'rank=',myrank,'nfault=',mesh%nelem_fault,&
    !'nfree=',mesh%nface_free

    allocate(mesh%fault2wave(mesh%nelem_fault))
    allocate(mesh%wave2fault(mesh%nelem))

    mesh%wave2fault(:) = 0

    k = 0
    do ie = 1,mesh%nelem
      isfault=.false.
      do is = 1,Nfaces
        if (mesh%bctype(is,ie) >= BC_FAULT) then
          isfault=.true.
        end if
      end do
      if (isfault) then
        k = k + 1
        mesh%fault2wave(k) = ie
        mesh%wave2fault(ie) = k
      end if
    end do

    allocate(mesh%tau_n   (Nfp,Nfaces,mesh%nelem_fault))
    allocate(mesh%tau_0   (Nfp,Nfaces,mesh%nelem_fault))
    allocate(mesh%stress  (Nfp,Nfaces,mesh%nelem_fault))
    allocate(mesh%Slip    (Nfp,Nfaces,mesh%nelem_fault))
    allocate(mesh%mSlip   (Nfp,Nfaces,mesh%nelem_fault))
    allocate(mesh%tSlip   (Nfp,Nfaces,mesh%nelem_fault))
    allocate(mesh%sliprate(Nfp,Nfaces,mesh%nelem_fault))
    allocate(mesh%mu_s    (Nfp,Nfaces,mesh%nelem_fault))
    allocate(mesh%mu_d    (Nfp,Nfaces,mesh%nelem_fault))
    allocate(mesh%Dc      (Nfp,Nfaces,mesh%nelem_fault))
    allocate(mesh%C0      (Nfp,Nfaces,mesh%nelem_fault))
    allocate(mesh%isslip  (Nfp,Nfaces,mesh%nelem_fault))

    mesh%slip    (:,:,:) = 0
    mesh%mslip   (:,:,:) = 0
    mesh%sliprate(:,:,:) = 0
    nelem_fault = 0
    !do element = 1,mesh%nelem
    do ief = 1,mesh%nelem_fault
        element = mesh%fault2wave(ief)
        do surface = 1,Nfaces
            !if(mesh%face(surface,element)==FACE_FAULT) then
            if(mesh%bctype(surface,element)>=BC_FAULT) then
                !nelem_fault = nelem_fault + 1
                !!!!print*,element,surface
                !!! constant in one element
                !%%%if (surface == 1) then
                !%%%    x = sum(mesh%vx(:,1,element))/NGLL
                !%%%    y = sum(mesh%vy(:,1,element))/NGLL
                !%%%else if (surface == 2) then
                !%%%    x = sum(mesh%vx(NGLL,:,element))/NGLL
                !%%%    y = sum(mesh%vy(NGLL,:,element))/NGLL
                !%%%else if (surface == 3) then
                !%%%    x = sum(mesh%vx(:,NGLL,element))/NGLL
                !%%%    y = sum(mesh%vy(:,NGLL,element))/NGLL
                !%%%else if (surface == 4) then
                !%%%    x = sum(mesh%vx(1,:,element))/NGLL
                !%%%    y = sum(mesh%vy(1,:,element))/NGLL
                !%%%end if 
                do i = 1,Nfp
                    if (surface == 1) then
                        x = mesh%vx(i,1,element)
                        y = mesh%vy(i,1,element)
                    else if (surface == 2) then
                        x = mesh%vx(Ngrid,i,element)
                        y = mesh%vy(Ngrid,i,element)
                    else if (surface == 3) then
                        x = mesh%vx(i,Ngrid,element)
                        y = mesh%vy(i,Ngrid,element)
                    else if (surface == 4) then
                        x = mesh%vx(1,i,element)
                        y = mesh%vy(1,i,element)
                    end if 

#if defined(TPV10) || defined(TPV11)
                    mesh%mu_s (i,surface,ief) = 0.76
                    mesh%mu_d (i,surface,ief) = 0.448
                    mesh%Dc   (i,surface,ief) = 0.5
                    mesh%C0   (i,surface,ief) = 0.2e6

                    asp_size = 1.5e3*(1.0+1e-15)*dsqrt(3d0)/2d0
                    if ( abs(y+10.392304845413264e3) < asp_size ) then
                        mesh%tau_0(i,surface,ief) = +81.24e6 + 1*0.36e6
                    endif

                    mesh%tau_n(i,surface,ief) = -7378.0*abs(y)/sqrt(3.0)*2.0
                    mesh%tau_0(i,surface,ief) = -0.55*(-mesh%tau_n(i,surface,ief))

                    if ( abs(y+10.392304845413264e3) < asp_size ) then
                        mesh%tau_0(i,surface,ief) = -mesh%C0(i,surface,ief)+ &
                            -(mesh%mu_s(i,surface,ief)+0.0057)*(-mesh%tau_n(i,surface,ief))
                    endif
                    if ( abs(y) > 12.990381056766578e3*(1.0+1e-5) ) then
                        mesh%mu_s(i,surface,ief) = 1e4
                        mesh%C0  (i,surface,ief) = 1e9
                    endif
#endif
#ifdef TPV11
                    mesh%mu_s (i,surface,ief) = 0.57
                    if ( abs(y+10.392304845413264e3) < asp_size ) then
                        mesh%tau_0(i,surface,ief) = -mesh%C0(i,surface,ief)+ &
                            -(mesh%mu_s(i,surface,ief)+0.0057)*(-mesh%tau_n(i,surface,ief))
                    endif
#endif

#ifdef TPV14
                    mesh%mu_s (i,surface,ief) = 0.677
                    mesh%mu_d (i,surface,ief) = 0.525
                    mesh%Dc   (i,surface,ief) = 0.4
                    mesh%C0   (i,surface,ief) = 0.0e6
                    mesh%tau_n(i,surface,ief) = -120e6
                    mesh%tau_0(i,surface,ief) = -70e6

                    asp_size = 1.5e3*(1.0+1e-15)
                    if ( abs(x+8e3) < asp_size ) then
                        mesh%tau_0(i,surface,ief) = -81.24e6 - 1*0.36e6
                    endif
                    ! 100 gap on the branch, main fault normal vector = (0,1)
                    !if ( &
                    !        abs(mesh%nx(i,surface,ief)) >1e-3 .and. &
                    !        abs(mesh%ny(i,surface,ief)) >1e-3 .and. &
                    !        y>-50 .and. &
                    !        x<86.602540378443862+1e-3 ) then
                    !    mesh%mu_s(i,surface,ief) = 1e4
                    !    mesh%C0  (i,surface,ief) = 1e9
                    !endif

                    vec_n = (/mesh%nx(i,surface,element),mesh%ny(i,surface,element)/)
                    vec_m = (/-vec_n(2),vec_n(1)/)

                    sxx = 0
                    syy = -120e6
                    sxy = 70e6
                    Tx = sxx*vec_n(1)+sxy*vec_n(2)
                    Ty = sxy*vec_n(1)+syy*vec_n(2)

                    call rotate_xy2nm(vec_n,vec_m,Tx,Ty,Tn,Tm)

                    mesh%tau_n(i,surface,ief) = Tn
                    mesh%tau_0(i,surface,ief) = Tm
                    !mesh%tau_0(i,surface,element) = Tm

                    !if (abs((abs(vec_n(2))-1d0))< 1e-6 ) then
                    !!print*,sngl(x),sngl(y),sngl(mesh%nx(i,surface,element)),sngl(mesh%ny(i,surface,element)),&
                    !print*,sngl(x),sngl(y),sngl(vec_n(1)),sngl(vec_n(2)),&
                    !            sngl(Tn),sngl(Tm)
                    !end if
#endif

#ifdef TPV5
                    mesh%mu_s (i,surface,ief) = 0.677
                    mesh%mu_d (i,surface,ief) = 0.525
                    mesh%Dc   (i,surface,ief) = 0.4
                    mesh%C0   (i,surface,ief) = 0.0e6

                    vec_n = (/mesh%nx(i,surface,element),mesh%ny(i,surface,element)/)
                    vec_m = (/-vec_n(2),vec_n(1)/)

                    sxx = -60e6
                    !syy = -120e6
                    !sxy = 70e6
                    syy = -60e6
                    sxy = 29.38e6

                    sxx = 0
                    syy = -120e6
                    sxy = 70e6
                    r1 = 1e3
                    r2 = 3e3
                    asp_size = r2*(1.0+1e-15)
                    !if ( abs(x+0e3) < asp_size ) then
                    !    !mesh%tau_0(i,surface,element) = -81.6e0
                    !    !sxy = 81.6e6
                    !    !sxy = -(mesh%mu_s(i,surface,element)+0.01)*syy
                    !    if (abs(x)<r1) then
                    !        damp = 1d0
                    !    else if (x>r1) then
                    !        damp = 0.5*(1+cos(pi*(x-r1)/(r2-r1)))
                    !    else if (x<-r1) then
                    !        damp = 0.5*(1+cos(pi*(x+r1)/(r2-r1)))
                    !    end if
                    !    sxy = 70e6 + damp * 11.6e6
                    !endif
                    asp_size = 1.5e3
                    if ( abs(x+0e3) < asp_size ) then
                        sxy = 81.6e6
                    endif

                    Tx = sxx*vec_n(1)+sxy*vec_n(2)
                    Ty = sxy*vec_n(1)+syy*vec_n(2)
                    call rotate_xy2nm(vec_n,vec_m,Tx,Ty,Tn,Tm)
                    mesh%tau_n(i,surface,ief) = Tn
                    mesh%tau_0(i,surface,ief) = Tm

                    !if ( abs(x+7.5e0) < asp_size ) then
                    !    mesh%tau_0(i,surface,element) = -78.0e0
                    !endif
                    !if ( abs(x-7.5e0) < asp_size ) then
                    !    mesh%tau_0(i,surface,element) = -62.0e0
                    !endif
                    !if ( x<-15e3 .or. x>15e3) then
                    !    mesh%mu_s(i,surface,element) = 1e4
                    !    mesh%C0  (i,surface,element) = 1e9
                    !endif
                    !print*,mesh%tau_0(i,surface,element)
#endif
#ifdef TPV28
                    mesh%mu_s (i,surface,ief) = 0.677
                    mesh%mu_d (i,surface,ief) = 0.373
                    mesh%Dc   (i,surface,ief) = 0.4
                    mesh%C0   (i,surface,ief) = 0.0e6

                    vec_n = (/mesh%nx(i,surface,element),mesh%ny(i,surface,element)/)
                    vec_m = (/-vec_n(2),vec_n(1)/)

                    sxx = -60e6
                    syy = -60e6
                    sxy = 29.38e6

                    asp_size = 1.5e3*(1.0+1e-15)
                    if ( abs(x+0e3) < asp_size ) then
                        sxy = -(mesh%mu_s(i,surface,ief)+0.01)*syy
                    endif

                    Tx = sxx*vec_n(1)+sxy*vec_n(2)
                    Ty = sxy*vec_n(1)+syy*vec_n(2)
                    call rotate_xy2nm(vec_n,vec_m,Tx,Ty,Tn,Tm)
                    mesh%tau_n(i,surface,ief) = Tn
                    mesh%tau_0(i,surface,ief) = Tm
#endif

                    !print*,sngl(x),sngl(y),sngl(mesh%nx(i,surface,element)),sngl(mesh%ny(i,surface,element)),&
                    tmp = abs(Tm)-mesh%mu_s(i,surface,ief)*abs(Tn)
                    if (tmp>0) then
                    !print*,'fault:',sngl(x),sngl(Tn*1e-6),sngl(Tm*1e-6),sngl(tmp*1e-6)
                    end if

                    !!!if (abs(y)>1e-3) then
                    !!!    mesh%mu_s (i,surface,element) = 1e4
                    !!!    mesh%C0   (i,surface,element) = 1e9
                    !!!end if
                enddo

                !!! constant in one element
                !%%%mesh%mu_s (:,surface,element) = 0.677
                !%%%mesh%tau_0(:,surface,element) = +70e6
                !%%%mesh%C0   (:,surface,element) = 0

                !%%%mesh%mu_s (:,surface,element) = 0.76
                !%%%mesh%C0   (:,surface,element) = 0.2e6
                !%%%!y = sum(mesh%vy(mesh%vmapM(1:NGLL,surface,element)))/NGLL
                !%%%if (surface == 1) then
                !%%%    x = sum(mesh%vx(:,1,element))/NGLL
                !%%%    y = sum(mesh%vy(:,1,element))/NGLL
                !%%%else if (surface == 2) then
                !%%%    x = sum(mesh%vx(NGLL,:,element))/NGLL
                !%%%    y = sum(mesh%vy(NGLL,:,element))/NGLL
                !%%%else if (surface == 3) then
                !%%%    x = sum(mesh%vx(:,NGLL,element))/NGLL
                !%%%    y = sum(mesh%vy(:,NGLL,element))/NGLL
                !%%%else if (surface == 4) then
                !%%%    x = sum(mesh%vx(1,:,element))/NGLL
                !%%%    y = sum(mesh%vy(1,:,element))/NGLL
                !%%%end if 
                !%%%asp_size = 1.5e3*(1.0+1e-15)*dsqrt(3d0)/2d0
                !%%%!if ( abs(y+10.3923e3) < asp_size .and. abs(x) < 1e-3 ) then
                !%%%if ( abs(y+10.3923e3) < asp_size ) then
                !%%%    mesh%tau_0(:,surface,element) = +81.24e6 + 1*0.36e6
                !%%%endif

                !%%%mesh%tau_n(:,surface,element) = -7378.0*abs(y)/sqrt(3.0)*2.0
                !%%%mesh%tau_0(:,surface,element) = -0.55*(-mesh%tau_n(:,surface,element))

                !%%%if ( abs(y+10.3923e3) < asp_size ) then
                !%%%    mesh%tau_0(:,surface,element) = -mesh%C0(:,surface,element)+ &
                !%%%        -(mesh%mu_s(:,surface,element)+0.0057)*(-mesh%tau_n(:,surface,element))
                !%%%endif
                !%%%!print*,'y=',sngl(y),'elem=',element,'face=',surface,'tau=',sngl(mesh%tau_0(1,surface,element))
                !%%%!print*,'y=',sngl(y),'tau=',sngl(mesh%tau_0(1,surface,element))
                !%%%!write(*,'(f10.2,1e20.10)') -y/sqrt(3.0)*2.0d-3,mesh%tau_0(1,surface,element)-0.76*mesh%tau_n(1,surface,element)
                !%%%!!!if ( y >= asp_size .and. y < asp_size + 2) then
                !%%%!!!    mesh%tau_0(((surface-1)*NGLL+1):(surface*NGLL),element) = (81.6-70.0)*0.5*(1+cos(pi*(y-asp_size)/2.0))+70.0
                !%%%!!!endif
                !%%%!!!if ( y < -asp_size .and. y >= -asp_size - 2) then
                !%%%!!!    mesh%tau_0(((surface-1)*NGLL+1):(surface*NGLL),element) = (81.6-70.0)*0.5*(1+cos(pi*(y+asp_size)/2.0))+70.0
                !%%%!!!endif
                !%%%if ( abs(y) > 13e3*(1.0+1e-5) ) then
                !%%%    mesh%mu_s(:,surface,element) = 1e4
                !%%%    mesh%C0  (:,surface,element) = 1e9
                !%%%endif
            else
                mesh%mu_s (:,surface,ief) = 1e4
                mesh%C0   (:,surface,ief) = 1e9
            endif

        enddo
    enddo

    mesh%stress = mesh%tau_0

    !stop 2
    !!!mesh%nelem_fault = nelem_fault
    !!!allocate(mesh%fault_elem(mesh%nelem_fault))
    !!!allocate(mesh%fault_face(mesh%nelem_fault))
    !!!i = 0
    !!!do element = 1,mesh%nelem
    !!!    do surface = 1,3
    !!!        !if(mesh%face(surface,element)==FACE_FAULT) then
    !!!        if(mesh%face_flag(surface,element)==FACE_FAULT) then
    !!!            i = i + 1
    !!!            mesh%fault_elem(i) = element
    !!!            mesh%fault_face(i) = surface
    !!!        endif
    !!!    enddo
    !!!enddo
    !!!!print*,mesh%fault_elem
    !!!!print*,mesh%fault_face
    !!!
    !!!open(100,file='fault_elem_face')
    !!!do i = 1,mesh%nelem_fault
    !!!    write(100,*) mesh%fault_elem(i),mesh%fault_face(i)
    !!!enddo
    !!!close(100)

    !!!!!! debug
    !!! if (2>3) then
    !!!    do element = 1,mesh%nelem
    !!!        do surface = 1,3
    !!!            !if(mesh%face(surface,element)==FACE_FAULT) then
    !!!            if(mesh%face_flag(surface,element)==FACE_FAULT) then
    !!!                write(*,'(2i6,5f10.3)')element,mesh%neigh(surface,element),&
    !!!                mesh%tau_0((surface-1)*NGLL+1:surface*NGLL,element)

    !!!            endif
    !!!        enddo
    !!!    enddo

    !!!    stop 8
    !!! endif

    !!!!!! debug
    !!! if (2>3) then
    !!! print*,'nelem=',mesh%nelem
    !!! do element = 1, mesh%nelem
    !!!     do surface = 1, 3
    !!!         !if ( mesh%face(surface, element) == FACE_FAULT) then
    !!!         if ( mesh%face_flag(surface, element) == FACE_FAULT) then
    !!!         do i = 1, NGLL
    !!!             vx1 = mesh%vx(mesh%vmapM(i,surface,element))
    !!!             vx2 = mesh%vx(mesh%vmapP(i,surface,element))
    !!!             vy1 = mesh%vy(mesh%vmapM(i,surface,element))
    !!!             vy2 = mesh%vy(mesh%vmapP(i,surface,element))
    !!!             if (abs(vx1-vx2) > 1e-5 .or. abs(vy1-vy2) > 1e-5) then
    !!!                 write(*,'(2i2,i6,2f15.6)')i,surface,element,vx1-vx2,vy1-vy2
    !!!                 !write(*,'(2i2,i6,4f15.6)')i,surface,element,vx1,vx2,vy1,vy2
    !!!             endif
    !!!         enddo
    !!!         !if (element == 498) then
    !!!         !    write(*,'(i2,5f15.6)')surface,mesh%vy(mesh%vmapM(:,surface,element))
    !!!         !    write(*,'(i2,5f15.6)')surface,mesh%vy(mesh%vmapP(:,surface,element))
    !!!         !endif

    !!!         endif
    !!!     enddo
    !!! enddo
    !!!
    !!! stop 9
    !!! endif
    k = 0
    do i = 1,mesh%nelem
        do is = 1,4
            if (mesh%bctype(is,i) >= BC_FAULT) then
                k = k + 1
                if (is == 1) then
                    xface = mesh%vx(:,1,i)
                    yface = mesh%vy(:,1,i)
                else if (is == 2) then
                    xface = mesh%vx(Ngrid,:,i)
                    yface = mesh%vy(Ngrid,:,i)
                else if (is == 3) then
                    xface = mesh%vx(:,Ngrid,i)
                    yface = mesh%vy(:,Ngrid,i)
                else if (is == 4) then
                    xface = mesh%vx(1,:,i)
                    yface = mesh%vy(1,:,i)
                end if
                mesh%fault_buffer(:,k,1) = xface
                mesh%fault_buffer(:,k,2) = yface
                mesh%fault_buffer(:,k,3) = mesh%nx(:,is,i)
                mesh%fault_buffer(:,k,4) = mesh%ny(:,is,i)
            end if
        end do
    end do

    k = 0
    do i = 1,mesh%nelem
        do is = 1,4
            if (mesh%bctype(is,i) == BC_FREE .or. mesh%bctype(is,i) == BC_FREE_G) then
                k = k + 1
                if (is == 1) then
                    xface = mesh%vx(:,1,i)
                    yface = mesh%vy(:,1,i)
                else if (is == 2) then
                    xface = mesh%vx(Ngrid,:,i)
                    yface = mesh%vy(Ngrid,:,i)
                else if (is == 3) then
                    xface = mesh%vx(:,Ngrid,i)
                    yface = mesh%vy(:,Ngrid,i)
                else if (is == 4) then
                    xface = mesh%vx(1,:,i)
                    yface = mesh%vy(1,:,i)
                end if
                mesh%surface_buffer(:,k,1) = xface
                mesh%surface_buffer(:,k,2) = yface
                mesh%surface_buffer(:,k,3) = mesh%nx(:,is,i)
                mesh%surface_buffer(:,k,4) = mesh%ny(:,is,i)
            end if
        end do
    end do

end subroutine fault_init

end module
