module mod_mesh

use mod_para
use mod_types

implicit none


contains

!old subroutine read_mesh(mesh,basename)
!old     type(meshvar) :: mesh
!old     integer :: i,j,k
!old     character(len=*),intent(in) :: basename
!old 
!old     !print*, 'read mesh: nodes, elements, neighbours, faces ...'
!old     ! read in nodes
!old     open(100,file=trim(adjustl(basename))//'.node')
!old     read(100,*) mesh%ncoord,i
!old     allocate(mesh%coord(2,mesh%ncoord))
!old     do j = 1,mesh%ncoord
!old     read(100,*) (mesh%coord(i,j),i=1,2)
!old     enddo
!old     close(100)
!old 
!old     ! read in elements
!old     open(100,file=trim(adjustl(basename))//'.elem')
!old     read(100,*) mesh%nelem,i
!old     allocate(mesh%elem(nsurface,mesh%nelem))
!old     do j = 1,mesh%nelem
!old     read(100,*) (mesh%elem(i,j),i=1,nsurface)
!old     enddo
!old     close(100)
!old 
!old     ! read in neighbours
!old     open(100,file=trim(adjustl(basename))//'.neigh')
!old     read(100,*) i,j
!old     allocate(mesh%neigh(nsurface,mesh%nelem))
!old     do j = 1,mesh%nelem
!old     read(100,*) (mesh%neigh(i,j),i=1,nsurface)
!old     enddo
!old     close(100)
!old 
!old     ! read in face number of neighbour elements
!old     open(100,file=trim(adjustl(basename))//'.face')
!old     read(100,*) i,j
!old     allocate(mesh%face(nsurface,mesh%nelem))
!old     !allocate(mesh%face_flag(nsurface,mesh%nelem))
!old     do j = 1,mesh%nelem
!old     read(100,*) (mesh%face(i,j),i=1,nsurface)
!old     enddo
!old     !do j = 1,mesh%nelem
!old     !do i = 1,nsurface
!old     !    mesh%face_flag(i,j) = mesh%face(i,j)
!old     !enddo
!old     !enddo
!old     close(100)
!old 
!old     ! read in direction of face neighbours
!old     open(100,file=trim(adjustl(basename))//'.direction')
!old     read(100,*) i,j
!old     allocate(mesh%direction(nsurface,mesh%nelem))
!old     do j = 1,mesh%nelem
!old     read(100,*) (mesh%direction(i,j),i=1,nsurface)
!old     enddo
!old     close(100)
!old 
!old     ! read in boundary condition types
!old     open(100,file=trim(adjustl(basename))//'.bctype')
!old     read(100,*) i,j
!old     allocate(mesh%bctype(nsurface,mesh%nelem))
!old     do j = 1,mesh%nelem
!old     read(100,*) (mesh%bctype(i,j),i=1,nsurface)
!old     enddo
!old     close(100)
!old 
!old     !! read free nodes
!old     !open(100,file='free')
!old     !read(100,*) mesh%nfree
!old     !allocate(mesh%free_node(mesh%nfree))
!old     !read(100,*) (mesh%free_node(i),i=1,mesh%nfree)
!old     !close(100)
!old 
!old     !! read fault nodes
!old     !open(100,file='fault')
!old     !read(100,*) mesh%nfault
!old     !allocate(mesh%fault_node(mesh%nfault))
!old     !read(100,*) (mesh%fault_node(i),i=1,mesh%nfault)
!old     !close(100)
!old 
!old     if (.false.) then
!old 
!old         print*,'ncoord=',mesh%ncoord
!old         do j = 1,mesh%ncoord
!old         print*,(mesh%coord(i,j),i=1,2)
!old         enddo
!old 
!old         print*,'nelem=',mesh%nelem
!old         do j = 1,mesh%nelem
!old         print*,(mesh%elem(i,j),i=1,nsurface)
!old         enddo
!old 
!old         print*,'EtoE=',mesh%nelem
!old         do j = 1,mesh%nelem
!old         print*,(mesh%neigh(i,j),i=1,nsurface)
!old         enddo
!old 
!old         print*,'EtoF=',mesh%nelem
!old         do j = 1,mesh%nelem
!old         print*,(mesh%face(i,j),i=1,nsurface)
!old         enddo
!old     endif
!old 
!old end subroutine read_mesh

subroutine readMeshVar(this,myrank)
    implicit none

    type(meshVar),intent(inout) :: this
    integer,intent(in) :: myrank
    integer :: ios
    !character(len=*) :: filename
    character(len=80) :: filename
    integer :: i,j,k,k1,k2,k3,k4,k5!,is
    !real(kind=rkind),dimension(NGLL) :: xface,yface

    write(filename,'(a,i6.6)') 'data/meshVar',myrank
    !print*,trim(filename)
    !open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown', iostat = ios)
    open(unit=27,file=trim(filename),status='old',iostat=ios)
    if (ios /= 0) then
        print*,'could not open: '//trim(filename)
    end if

    read(27,*) this%ncoord

    allocate(this%coord(2,this%ncoord))

    read(27,*) this%coord

    read(27,*) this%nelem

    allocate(this%elem     (4,this%nelem))
    allocate(this%neigh    (4,this%nelem))
    allocate(this%face     (4,this%nelem))
    allocate(this%direction(4,this%nelem))
    allocate(this%bctype   (4,this%nelem))
    allocate(this%fluxtype (4,this%nelem))
    allocate(this%elemtype (  this%nelem))
    allocate(this%rho      (  this%nelem))
    allocate(this%vp       (  this%nelem))
    allocate(this%vs       (  this%nelem))

    !this%elemtype = ELEM_FLUID
    !this%elemtype = ELEM_SOLID

    read(27,*) this%elem
    read(27,*) this%neigh
    read(27,*) this%face
    read(27,*) this%direction
    read(27,*) this%bctype
    read(27,*) this%fluxtype
    read(27,*) this%elemtype
    read(27,*) this%rho
    read(27,*) this%vp
    read(27,*) this%vs

    read(27,*) this%mpi_nn
    read(27,*) this%mpi_ne
    read(27,*) this%mpi_nemax
    read(27,*) this%pinterfaces

    allocate(this%mpi_neighbor(this%mpi_nn))
    allocate(this%mpi_connection(this%mpi_nn,this%mpi_ne,2))
    allocate(this%mpi_ibool(nsurface,this%mpi_nemax))
    allocate(this%mpi_interface(4,nsurface,this%nelem))
    allocate(this%mpi_rho(this%pinterfaces,3))
    allocate(this%mpi_vp(this%pinterfaces,3))
    allocate(this%mpi_vs(this%pinterfaces,3))

    read(27,*) this%mpi_neighbor
    read(27,*) this%mpi_connection
    read(27,*) this%mpi_ibool
    read(27,*) this%mpi_interface
    if (this%nproc > 1) then
        read(27,*) this%mpi_rho
        read(27,*) this%mpi_vp
        read(27,*) this%mpi_vs
    end if


    read(27,*) this%nrecv
    if (this%nrecv > 0) then
        allocate(this%recv_fid   (this%nrecv))
        allocate(this%recv_i     (this%nrecv))
        allocate(this%recv_ie    (this%nrecv))
        allocate(this%recv_refx  (this%nrecv))
        allocate(this%recv_buffer(this%nrecv,Ngrid,10))

        read(27,*) this%recv_fid
        read(27,*) this%recv_i
        read(27,*) this%recv_ie
        read(27,*) this%recv_refx
    end if

    read(27,*) this%body_nrecv
   ! print*,'rank=',myrank,'body_nrecv=',this%body_nrecv
    if (this%body_nrecv > 0) then
        allocate(this%body_recv_fid (this%body_nrecv))
        allocate(this%body_recv_i   (this%body_nrecv))
        allocate(this%body_recv_j   (this%body_nrecv))
        allocate(this%body_recv_ie  (this%body_nrecv))
        allocate(this%body_recv_refx(this%body_nrecv))
        allocate(this%body_recv_refy(this%body_nrecv))

        read(27,*) this%body_recv_fid
        read(27,*) this%body_recv_i
        read(27,*) this%body_recv_j
        read(27,*) this%body_recv_ie
        read(27,*) this%body_recv_refx
        read(27,*) this%body_recv_refy
    end if

    close(27)

    !print*,sngl(this%coord)
    !print*, 'elem      =',        this%elem
    !print*, 'neigh     =',        this%neigh
    !print*, 'face      =',        this%face
    !print*, 'direction =',        this%direction
    !print*,'mpi_nn=',this%mpi_nn
    !print*,'mpi_ne=',this%mpi_ne
    !print*,'mpi_nemax=',this%mpi_nemax
    !print*,'mpi_neighbor=',this%mpi_neighbor
    !print*,'mpi_connection=',this%mpi_connection
    !print*,'mpi_ibool=',this%mpi_ibool
    !print*,'mpi_interface=',this%mpi_interface

    allocate(this%eta (Ngrid,Nfaces,this%nelem))
    allocate(this%deta(Ngrid,Nfaces,this%nelem))

    this%eta(:,:,:) = 0
    this%deta(:,:,:) = 0

    k1=0; k2=0; k3=0; k4=0;
    do i = 1,this%nelem
        do j = 1,4
            if (this%bctype(j,i) >= BC_FAULT) then
                k1=k1+1
            end if
            if (this%bctype(j,i) == BC_FREE .or. this%bctype(j,i) == BC_FREE_G) then
                k2=k2+1
            end if
            if (this%elemtype(i)==ELEM_SOLID .and. this%bctype(j,i) == BC_SOLID_FLUID) then
                k3=k3+1
            end if
            if (this%elemtype(i)==ELEM_FLUID .and. this%bctype(j,i) == BC_SOLID_FLUID) then
                k4=k4+1
            end if
        end do
    end do
    this%nface_fault   = k1
    this%nface_free    = k2
    this%nface_inter_s = k3
    this%nface_inter_f = k4

    allocate(this%fault_buffer(Ngrid,this%nface_fault,5))
    allocate(this%surface_buffer(Ngrid,this%nface_free,5))

    print*,'rank=',myrank,'ncoord=',this%ncoord,'nelem=',this%nelem,'fault face=',this%nface_fault
    print*,'rank=',myrank,'ncoord=',this%ncoord,'nelem=',this%nelem,'free face=',this%nface_free
    print*,'rank=',myrank,'ncoord=',this%ncoord,'nelem=',this%nelem,'interface(solid) face=',this%nface_inter_s
    print*,'rank=',myrank,'ncoord=',this%ncoord,'nelem=',this%nelem,'interface(fluid) face=',this%nface_inter_f
    !stop 2
    !stop 222
end subroutine
end module
