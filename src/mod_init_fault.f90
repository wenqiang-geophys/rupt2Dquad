module mod_init_fault

  use mod_para
  use mod_types
  use mod_funcs
  use netcdf
  use mod_check

  implicit none

contains


subroutine fault_init_external(mesh)
  implicit none
  type(meshvar) :: mesh
  integer :: is,ie,ief,i
  integer :: myrank
  integer :: ncid,varid,ierr
  character(len=128) :: filename
  real(kind=rkind) :: vec_n(2),vec_m(2)
  real(kind=rkind) :: Tx,Ty,Tn,Tm,Sxx,Syy,Sxy,sigma,tau
  real(kind=rkind),allocatable,dimension(:,:,:) :: Sxx0,Syy0,Sxy0,sigma0,tau0
  real(kind=rkind),allocatable,dimension(:,:,:) :: mu_s,mu_d,Dc,C0
  real(kind=rkind) :: c(2),cp,p(2),v1(2),v2(2)

  integer :: FToV(2,4)

  if (mesh%nelem_fault == 0) return;

  ! point outward
  ! Face 1: 1,2
  ! Face 2: 2,3
  ! Face 3: 3,4
  ! Face 4: 4,1
  FToV = reshape((/1,2,2,3,3,4,4,1/),shape(FToV))

  myrank = mesh%rank

  allocate(Sxx0(2,Nfaces,mesh%nelem_fault))
  allocate(Syy0(2,Nfaces,mesh%nelem_fault))
  allocate(Sxy0(2,Nfaces,mesh%nelem_fault))
  allocate(sigma0(2,Nfaces,mesh%nelem_fault))
  allocate(tau0(2,Nfaces,mesh%nelem_fault))
  allocate(mu_s(2,Nfaces,mesh%nelem_fault))
  allocate(mu_d(2,Nfaces,mesh%nelem_fault))
  allocate(Dc  (2,Nfaces,mesh%nelem_fault))
  allocate(C0  (2,Nfaces,mesh%nelem_fault))

  write(filename,'(a,a,i6.6,a)') trim(mesh_dir),'/meshVar',myrank,'.nc'
  ierr = nf90_open(trim(filename),NF90_NOWRITE,ncid)

  ierr = nf90_inq_varid(ncid,'Sxx0',varid)
  call check2(ierr,'inq_varid Sxx0')
  ierr = nf90_get_var(ncid,varid,Sxx0)
  call check2(ierr,'get_var Sxx0')

  ierr = nf90_inq_varid(ncid,'Syy0',varid)
  call check2(ierr,'inq_varid Syy0')
  ierr = nf90_get_var(ncid,varid,Syy0)
  call check2(ierr,'get_var Syy0')

  ierr = nf90_inq_varid(ncid,'Sxy0',varid)
  call check2(ierr,'inq_varid Sxy0')
  ierr = nf90_get_var(ncid,varid,Sxy0)
  call check2(ierr,'get_var Sxy0')

  ierr = nf90_inq_varid(ncid,'sigma0',varid)
  call check2(ierr,'inq_varid sigma0')
  ierr = nf90_get_var(ncid,varid,sigma0)
  call check2(ierr,'get_var sigma0')

  ierr = nf90_inq_varid(ncid,'tau0',varid)
  call check2(ierr,'inq_varid tau0')
  ierr = nf90_get_var(ncid,varid,tau0)
  call check2(ierr,'get_var tau0')

  ierr = nf90_inq_varid(ncid,'mu_s',varid)
  call check2(ierr,'inq_varid mu_s')
  ierr = nf90_get_var(ncid,varid,mu_s)
  call check2(ierr,'get_var mu_s')

  ierr = nf90_inq_varid(ncid,'mu_d',varid)
  call check2(ierr,'inq_varid mu_d')
  ierr = nf90_get_var(ncid,varid,mu_d)
  call check2(ierr,'get_var mu_d')

  ierr = nf90_inq_varid(ncid,'Dc',varid)
  call check2(ierr,'inq_varid Dc')
  ierr = nf90_get_var(ncid,varid,Dc)
  call check2(ierr,'get_var Dc')

  ierr = nf90_inq_varid(ncid,'C0',varid)
  call check2(ierr,'inq_varid C0')
  ierr = nf90_get_var(ncid,varid,C0)
  call check2(ierr,'get_var C0')

  do ief = 1,mesh%nelem_fault
    ie = mesh%fault2wave(ief)
    do is = 1,Nfaces
      ! linear interpolation,
      if (mesh%bctype(is,ie) >= BC_FAULT) then

        v1=mesh%coord(1:2,mesh%elem(FtoV(1,is),ie))
        v2=mesh%coord(1:2,mesh%elem(FtoV(2,is),ie))

        do i = 1,Nfp
          !p(1)=mesh%vx(mesh%vmapM(i,is,ie))
          !p(2)=mesh%vy(mesh%vmapM(i,is,ie))
          if (is == 1) then
              p(1) = mesh%vx(i,1,ie)
              p(2) = mesh%vy(i,1,ie)
          else if (is == 2) then
              p(1) = mesh%vx(Ngrid,i,ie)
              p(2) = mesh%vy(Ngrid,i,ie)
          else if (is == 3) then
              p(1) = mesh%vx(i,Ngrid,ie)
              p(2) = mesh%vy(i,Ngrid,ie)
          else if (is == 4) then
              p(1) = mesh%vx(1,i,ie)
              p(2) = mesh%vy(1,i,ie)
          end if

          c=mu_s(1:2,is,ief)
          cp=interp1d(v1,v2,c(1),c(2),p)
          mesh%mu_s(i,is,ief)=cp

          c=mu_d(1:2,is,ief)
          cp=interp1d(v1,v2,c(1),c(2),p)
          mesh%mu_d(i,is,ief)=cp

          c=Dc(1:2,is,ief)
          cp=interp1d(v1,v2,c(1),c(2),p)
          mesh%Dc(i,is,ief)=cp

          c=C0(1:2,is,ief)
          cp=interp1d(v1,v2,c(1),c(2),p)
          mesh%C0(i,is,ief)=cp

          if (input_stress_type==0) then
            c=Sxx0(1:2,is,ief)
            Sxx=interp1d(v1,v2,c(1),c(2),p)

            c=Syy0(1:2,is,ief)
            Syy=interp1d(v1,v2,c(1),c(2),p)

            c=Sxy0(1:2,is,ief)
            Sxy=interp1d(v1,v2,c(1),c(2),p)
          elseif (input_stress_type==1) then
            c=sigma0(1:2,is,ief)
            sigma=interp1d(v1,v2,c(1),c(2),p)

            c=tau0(1:2,is,ief)
            tau=interp1d(v1,v2,c(1),c(2),p)
          endif

          vec_n = (/mesh%nx(i,is,ie),mesh%ny(i,is,ie)/)
          vec_m = (/-vec_n(2),vec_n(1)/)

          if (input_stress_type==0) then
          Tx = sxx*vec_n(1)+sxy*vec_n(2)
          Ty = sxy*vec_n(1)+syy*vec_n(2)

          call rotate_xy2nm(vec_n,vec_m,Tx,Ty,Tn,Tm)

          mesh%tau_n(i,is,ief) = Tn
          mesh%tau_0(i,is,ief) = Tm
          elseif (input_stress_type==1) then
          mesh%tau_n(i,is,ief) = sigma
          mesh%tau_0(i,is,ief) = tau
          endif

          ! write initial stress
          !mesh%stress1(i,is,ief) = Tm+0*dTm
          !mesh%stress2(i,is,ief) = Tl+0*dTl
          !mesh%stress(i,is,ief) = sqrt((Tm+0*dTm)**2+(Tl+0*dTl)**2)
          !mesh%sigma (i,is,ief) = Tn+0*dTn
        end do

      else! not fault
        mesh%mu_s(:,is,ief) = 1e4
        mesh%C0  (:,is,ief) = 1e9
      end if
    end do
  end do

  ierr = nf90_close(ncid)

  deallocate(Sxx0,Syy0,Sxy0)
  deallocate(mu_s,mu_d,Dc,C0)
end subroutine

end module
