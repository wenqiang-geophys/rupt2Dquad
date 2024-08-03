module mod_smooth

use mod_para
use mod_gll
use mod_fd
use mod_mesh

contains

subroutine smooth_fault(mesh,u)
    implicit none
    !type(meshvar),intent(inout) :: mesh
    type(meshvar),intent(in) :: mesh
    real(kind=rkind),intent(inout) :: u(:,:,:)
    real(kind=rkind),dimension(Ngrid,Ngrid,8) :: u1,u2,u3
    integer :: i,j,ii,jj,is,ie,iv,N
    real(kind=rkind) :: beta,Dc,slip
    real(kind=rkind) :: wts(NGLL,NGLL)

    !do i = 1,NGLL
    !    wts(i) = exp(-(i-1.0)**2/0.5**2)
    !end do
    !wts(1) = 2.0
    !wts(2:NGLL) = 0.0
    !wts = wts/sum(wts)

    return

    !print*,wts

    N = NGLL
    do ie = 1,mesh%nelem
        do is = 1,Nfaces
            if (mesh%bctype(is,ie) == BC_FAULT) then
                Dc = sum(mesh%Dc(:,is,ie))/dble(NGLL)
                slip = sum(mesh%slip(:,is,ie))/dble(NGLL)
                beta = min(slip/Dc,1.0)
                !if(beta>0) print*,beta

                !u1 = reshape(u(:,ie,1:8),(/Ngrid,Ngrid,8/))
                u1 = reshape(u(:,ie,1:8),(/N,N,8/))
                u2 = u1
                !u2(:,:,:) = 0.0
                
                !       3
                !    4------3
                !    |      | 2
                !  4 |      | 
                !    1------2
                !       1
                !mesh%vmapM(:,1) = iv2d(:,1)
                !mesh%vmapM(:,2) = iv2d(N,:)
                !mesh%vmapM(:,3) = iv2d(:,N)
                !mesh%vmapM(:,4) = iv2d(1,:)
                if(is==1) then
                    u2(:,1,:) = 0.
                    do i = 1,NGLL
                        do ii = 1,NGLL
                        do jj = 1,NGLL
                            wts(ii,jj) = exp( -( (ii-i)**2+(jj-1)**2 )/0.5**2  )
                        enddo
                        enddo
                        wts = wts/sum(wts)
                        do ii = 1,NGLL
                            do jj = 1,NGLL
                                u2(i,1,:) = u2(i,1,:) + wts(ii,jj) * u1(ii,jj,:)
                            end do
                        end do
                    end do
                elseif(is==2) then
                    u2(N,:,:) = 0.
                    do j = 1,NGLL
                        do ii = 1,NGLL
                        do jj = 1,NGLL
                            wts(ii,jj) = exp( -( (ii-N)**2+(jj-j)**2 )/0.5**2  )
                        enddo
                        enddo
                        wts = wts/sum(wts)
                        do ii = 1,NGLL
                            do jj = 1,NGLL
                                u2(N,j,:) = u2(N,j,:) + wts(ii,jj) * u1(ii,jj,:)
                            end do
                        end do
                    end do
                elseif(is==3) then
                    u2(:,N,:) = 0.
                    do i = 1,NGLL
                        do ii = 1,NGLL
                        do jj = 1,NGLL
                            wts(ii,jj) = exp( -( (ii-i)**2+(jj-N)**2 )/0.5**2  )
                        enddo
                        enddo
                        wts = wts/sum(wts)
                        do ii = 1,NGLL
                            do jj = 1,NGLL
                                u2(i,N,:) = u2(i,N,:) + wts(ii,jj) * u1(ii,jj,:)
                            end do
                        end do
                    end do
                else!if (is==4) then
                    u2(1,:,:) = 0.
                    do j = 1,NGLL
                        do ii = 1,NGLL
                        do jj = 1,NGLL
                            wts(ii,jj) = exp( -( (ii-1)**2+(jj-j)**2 )/0.5**2  )
                        enddo
                        enddo
                        wts = wts/sum(wts)
                        do ii = 1,NGLL
                            do jj = 1,NGLL
                                u2(1,j,:) = u2(1,j,:) + wts(ii,jj) * u1(ii,jj,:)
                            end do
                        end do
                    end do
                endif

                !u2 = u1

                u3 = beta*u2 + (1.0-beta)*u1
                ! only smooth stress
                u3(:,:,1:2) = u1(:,:,1:2)

                !u3 = u1

                !u3(:,:,1:2) = u1(:,:,1:2)
                !u3(:,:,6:8) = u1(:,:,6:8)

                u(:,ie,1:8) = reshape(u3,(/Np,8/))
            end if
        end do ! face
    end do ! elem

end subroutine

end module
