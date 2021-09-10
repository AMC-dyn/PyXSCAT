module integrals
    implicit none
    contains


    subroutine tot_integration(ncap,px,py,pz,l,m,n,p0matrix,dx,dy,dz,z1,z2,group_start,group_count,group, &
            cutoffz,cutoffmd, cutoffcentre,q,e12,tsi)


        use omp_lib
        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(kind=ikind), INTENT(IN) :: ncap
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:) :: l,m,n,group,group_start,group_count

        REAL(kind=dp), intent(in), dimension(:,:,:,:,:),target :: dx,dy,dz
        real(kind=dp),dimension(:,:,:), pointer :: dx1,dy1,dz1,dx2,dy2,dz2
        REAL(kind=dp), intent(in), dimension(:,:,:,:) :: p0matrix
        REAL(kind=dp), intent(in), dimension(:,:,:) ::e12
        real(kind=dp), intent(in), dimension(:,:,:)::z1,z2
        real(kind=dp), dimension(:,:,:,:), allocatable::zcontrred,zcontrred2
        REAL(kind=dp), intent(in), dimension(:,:) :: px,py,pz
        REAL(kind=dp), intent(in), dimension(:) :: q
        REAL(kind=dp), intent(in) :: cutoffz, cutoffmd,cutoffcentre


        REAL(kind=dp), dimension(size(q)) :: f
        integer(kind=ikind), dimension(:,:), allocatable :: posits
        real(kind=dp),dimension(:,:), allocatable :: za,zb,cmat
        integer(kind=ikind),dimension(:), allocatable ::posi,posj,posk,posr
        REAL(kind=dp), intent(out), dimension(size(q)) :: tsi
        real(kind=dp) :: hx,hy,hz,h
        integer(kind=ikind),dimension(size(l)) :: ll
        integer(kind=ikind) :: nq,i,j,k,r,count,ng,ii,jj
        integer(kind=ikind) :: spi, spj, spk, spr, szo,nt

        ng=maxval(group)
        ll=l+m+n
        allocate(posits(ng,maxval(group_count)))
        posits=1
        count=1
        do i=1,ncap
            count=1
            do j=group_start(i),group_start(i)+group_count(i)-1
               posits(i,count)=j
                count=count+1
            end do
        end do
        print*,'posits created'
        print*,size(group_start),sum(group_count)


        szo = size(z1(:,1,1))

        nq= size(q)

        print*,OMP_get_num_threads()
        call omp_set_num_threads(32)
        print*,OMP_get_num_threads()
        !First big loop

        tsi=0.0_dp
       !$OMP PARALLEL do private(posI,posK,posJ,posR,spi,spj,spk,spr,zcontrred,zcontrred2,za,zb,cmat), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,j,k,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz, posits,cutoffmd,group_count,group_start) REDUCTION(+:tsi)

        do i=1,ncap
            do j=i+1,ncap
                do k=i+1,ncap
                    do r=k+1,ncap
                        hx = px(k, r) - px(i, j)
                        hy = py(k, r) - py(i, j)
                        hz = pz(k, r) - pz(i, j)
                        h = (hx * hx + hy * hy + hz * hz)**0.5

                        allocate(posI(size(posits(i,:group_count(i)))),posJ(size(posits(j,:group_count(j)))) &
                                ,posK(size(posits(k,:group_count(k)))),posR(size(posits(r,:group_count(r)))))


                        posI = posits(i,:group_count(i))
                        posJ = posits(j,:group_count(j))
                        posK = posits(k,:group_count(k))
                        posR = posits(r,:group_count(r))

                        spi = size(posI)
                        spj = size(posJ)
                        spk = size(posK)
                        spr = size(posR)

                        allocate(zcontrred(spj, spk, spr, spi), za(szo, spj), zb(szo, spk), &
                       &         cmat(spj, spk), zcontrred2(spr, spi, spj, spk))

                          do ii = 1, spi
                            do jj = 1, spr
                                za = z1(:,posj,posi(ii))
                              !  za = transpose(z1(:,posj,posi(ii)))
                                zb = z2(:,posk,posr(jj))
                               ! cmat = matmul(za,zb)
                                call dgemm('t','n', spj, spk, szo, 1.0_dp/8.0_dp, za, &
                              &           szo, zb, szo, 0.0_dp, cmat, spj)
                                zcontrred(:,:,jj,ii) = cmat

                                za = z2(:,posj,posi(ii))
                              !  za = transpose(z1(:,posj,posi(ii)))
                                zb = z1(:,posk,posr(jj))
                               ! cmat = matmul(za,zb)
                                call dgemm('t','n', spj, spk, szo, 1.0_dp/8.0_dp, za, &
                              &           szo, zb, szo, 0.0_dp, cmat, spj)

                                zcontrred2(jj,ii,:,:) = cmat
                            enddo
                          enddo



                         deallocate(posI,posJ,posK,posR,za,zb,cmat)
                         dx1=>dx(:,:,:,j,i)
                         dy1=>dy(:,:,:,j,i)
                         dz1=>dz(:,:,:,j,i)
                         dx2=>dx(:,:,:,r,k)
                         dy2=>dy(:,:,:,r,k)
                         dz2=>dz(:,:,:,r,k)

                        if (h < cutoffcentre) then

                            call tot_integral_ijkr_pzero(nq, l,m,n,group_start, group_count, p0matrix, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, k, r, &
                                    zcontrred, zcontrred2,  cutoffz, cutoffmd,f)



                        else

                            call tot_integral_k_ijkr(q,l,m,n,group_start, group_count, hx, hy, hz, h,dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, k, r, &
                                    zcontrred, zcontrred2,  cutoffz, cutoffmd,f)





                        end if
                        tsi = tsi + 8.000 * f * e12(:, i, j) * e12(:, k, r)


                        deallocate(zcontrred, zcontrred2)
               !        deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)


                    end do
                end do
            end do
        end do

       !$OMP END parallel DO


        print*,'la rata ', tsi(1)

        !$OMP PARALLEL do private(posI,posJ,posR,spi,spj,spk,spr,zcontrred,zcontrred2,za,zb,cmat), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,j,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz,posits, cutoffmd,group_count,group_start) REDUCTION(+:tsi)
        do i=1,ncap
            do j=i+1,ncap
                do r=i+1,ncap
                    hx = px(i, r) - px(i, j)
                    hy = py(i, r) - py(i, j)
                    hz = pz(i, r) - pz(i, j)
                    h = sqrt((hx * hx + hy * hy + hz * hz))
                    allocate(posI(size(posits(i,:group_count(i)))),posJ(size(posits(j,:group_count(j)))) &
                               ,posR(size(posits(r,:group_count(r)))))


                        posI = posits(i,:group_count(i))
                        posJ = posits(j,:group_count(j))

                        posR = posits(r,:group_count(r))

                        spi = size(posI)
                        spj = size(posJ)
                        spr = size(posR)

                        allocate(zcontrred(spj, spi, spr, spi), za(szo, spj), zb(szo, spi), &
                       &         cmat(spj, spi), zcontrred2(spr, spi, spj, spi))

                          do ii = 1, spi
                            do jj = 1, spr
                              ! za = transpose(z1(:,posj,posi(ii)))
                                za = z1(:,posj,posi(ii))
                                zb = z2(:,posi,posr(jj))
                               !cmat = matmul(za,zb)
                                call dgemm('t','n', spj, spi, szo, 1.0_dp/8.0_dp, za, &
                               &           szo, zb, szo, 0.0_dp, cmat, spj)
                                zcontrred(:,:,jj,ii) = cmat
                                  za = z2(:,posj,posi(ii))
                              !  za = transpose(z1(:,posj,posi(ii)))
                                zb = z1(:,posi,posr(jj))
                               ! cmat = matmul(za,zb)
                                call dgemm('t','n', spj, spi, szo, 1.0_dp/8.0_dp, za, &
                              &           szo, zb, szo, 0.0_dp, cmat, spj)

                                zcontrred2(jj,ii,:,:) = cmat
                            enddo
                          enddo

!
                    dx1=>dx(:,:,:,j,i)
                    dy1=>dy(:,:,:,j,i)
                    dz1=>dz(:,:,:,j,i)
                    dx2=>dx(:,:,:,r,i)
                    dy2=>dy(:,:,:,r,i)
                    dz2=>dz(:,:,:,r,i)


                    if (h < cutoffcentre) then
                        call tot_integral_ijkr_pzero(nq, l,m,n,group_start, group_count, p0matrix, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, i, r, &
                               zcontrred,  zcontrred2,  cutoffz, cutoffmd,f)
                    else

                        call tot_integral_k_ijkr(q, l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, i, r, &
                                zcontrred,  zcontrred2,  cutoffz, cutoffmd,f)



                    end if
                    tsi = tsi + 4.000 * f * e12(:, i, j) * e12(:, i, r)
                    count=count+1
                    deallocate(posI,posJ,posR,za,zb,cmat)
                    deallocate(zcontrred, zcontrred2)
                  !   deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

                end do
            end do

        end do
    !$OMP END parallel DO
        print*,'here'
      !$OMP PARALLEL do private(posI,posK,posR,spi,spj,spk,spr,zcontrred,zcontrred2,za,zb,cmat), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,k,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz,posits, cutoffmd,group_count,group_start) REDUCTION(+:tsi)
        do i=1,ncap
            do k=1,ncap
                do r=k+1,ncap
                    hx = px(k, r) - px(i, i)
                    hy = py(k, r) - py(i, i)
                    hz = pz(k, r) - pz(i, i)
                    h = sqrt((hx * hx + hy * hy + hz * hz))

                    allocate(posI(size(posits(i,:group_count(i)))) &
                                ,posK(size(posits(k,:group_count(k)))),posR(size(posits(r,:group_count(r)))))


                        posI = posits(i,:group_count(i))

                        posK = posits(k,:group_count(k))
                        posR = posits(r,:group_count(r))

                    spi = size(posI)
                    spk = size(posK)
                    spr = size(posR)

                    allocate(zcontrred(spi, spk, spr, spi), za(szo, spi), zb(szo, spk), &
                   &         cmat(spi, spk), zcontrred2(spr, spi, spi, spk))

                    do ii = 1, spi
                        do jj = 1, spr
                          !  za = transpose(z1(:,posi,posi(ii)))
                            za = z1(:,posi,posi(ii))
                            zb = z2(:,posk,posr(jj))
                           ! cmat = matmul(za,zb)
                            call  dgemm('t','n', spi, spk, szo, 1.0_dp/8.0_dp, za, &
                           &           szo, zb, szo, 0.0_dp, cmat, spi)
                            zcontrred(:,:,jj,ii) = cmat
                             za = z2(:,posi,posi(ii))
                              !  za = transpose(z1(:,posj,posi(ii)))
                                zb = z1(:,posk,posr(jj))
                               ! cmat = matmul(za,zb)
                                call dgemm('t','n', spi, spk, szo, 1.0_dp/8.0_dp, za, &
                              &           szo, zb, szo, 0.0_dp, cmat, spi)

                                zcontrred2(jj,ii,:,:) = cmat
                        enddo
                    enddo


                    dx1=>dx(:,:,:,i,i)
                    dy1=>dy(:,:,:,i,i)
                    dz1=>dz(:,:,:,i,i)
                    dx2=>dx(:,:,:,r,k)
                    dy2=>dy(:,:,:,r,k)
                    dz2=>dz(:,:,:,r,k)
!                    zcontrred=zcontrred
!                    zcontrred2=zcontrred2

                    if (h < cutoffcentre) then
                        call tot_integral_ijkr_pzero(nq,l,m,n,group_start, group_count, p0matrix, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, k, r, &
                                    zcontrred, zcontrred2, cutoffz, cutoffmd, f)
                    else

                        call tot_integral_k_ijkr(q, l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, k, r, &
                                    zcontrred, zcontrred2,  cutoffz, cutoffmd, f)



                    end if
                    tsi = tsi+ 4.000 * f * e12(:, i, i) * e12(:, k, r)
                    count=count+1
                    deallocate(posI,posK,posR,za,zb,cmat)
                    deallocate(zcontrred, zcontrred2)
               !     deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

                end do
            end do
        end do

        !$OMP END parallel DO
        print*,'here 2'
          !$OMP PARALLEL do private(posI,posK,spi,spk,zcontrred,zcontrred2,za,zb,cmat), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,k,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
          !$OMP& shared( cutoffz, cutoffmd,posits,group_count,group_start) REDUCTION(+:tsi)
        do i=1,ncap
            do k=i+1,ncap

                hx = px(k, k) - px(i, i)
                hy = py(k, k) - py(i, i)
                hz = pz(k, k) - pz(i, i)
                h = sqrt((hx * hx + hy * hy + hz * hz))

                allocate(posI(size(posits(i,:group_count(i)))), &
                        posK(size(posits(k,:group_count(k)))))


                posI = posits(i,:group_count(i))

                posK = posits(k,:group_count(k))


                spi = size(posI)
                spk = size(posK)

                allocate(zcontrred(spi, spk, spk, spi), za(szo, spi), zb(szo, spk), &
               &         cmat(spi, spk), zcontrred2(spk, spi, spi, spk))

                do ii = 1, spi
                    do jj = 1, spk
                       ! za = transpose(z1(:,posi,posi(ii)))
                        za = z1(:,posi,posi(ii))
                        zb = z2(:,posk,posk(jj))
!                        cmat = matmul(za,zb)
                        call dgemm('t','n', spi, spk, szo, 1.0_dp/8.0_dp, za, &
                       &           szo, zb, szo, 0.0_dp, cmat, spi)
                        zcontrred(:,:,jj,ii) = cmat
                        za = z2(:,posi,posi(ii))
                              !  za = transpose(z1(:,posj,posi(ii)))
                        zb = z1(:,posk,posk(jj))
                               ! cmat = matmul(za,zb)
                        call dgemm('t','n', spi, spk, szo, 1.0_dp/8.0_dp, za, &
                              &           szo, zb, szo, 0.0_dp, cmat, spi)

                        zcontrred2(jj,ii,:,:) = cmat
                    enddo
                enddo


                dx1=>dx(:,:,:,i,i)
                dy1=>dy(:,:,:,i,i)
                dz1=>dz(:,:,:,i,i)
                dx2=>dx(:,:,:,k,k)
                dy2=>dy(:,:,:,k,k)
                dz2=>dz(:,:,:,k,k)


                if (h < cutoffcentre) then
                    call tot_integral_ijkr_pzero(nq, l,m,n,group_start, group_count, p0matrix, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, k, k, &
                                    zcontrred, zcontrred2,  cutoffz, cutoffmd, f)
                else

                    call tot_integral_k_ijkr(q,l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, k, k, &
                                    zcontrred, zcontrred2,  cutoffz, cutoffmd, f)
                end if
                tsi = tsi+ 2.000 * f * e12(:, i, i) * e12(:, k, k)
                deallocate(posI,posK,za,zb,cmat)
                deallocate(zcontrred, zcontrred2)
               ! deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

            end do
        end do


        !$OMP END parallel DO
        print*,'here 3'
        !$OMP PARALLEL do private(posI,spi,zcontrred,zcontrred2,za,zb,cmat), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,ll, p0matrix), &
          !$OMP& shared( cutoffz, cutoffmd,posits,group_count,group_start) REDUCTION(+:tsi)
        do i=1,ncap
                              allocate(posI(size(posits(i,:group_count(i)))))


                        posI = posits(i,:group_count(i))



            spi = size(posI)

            allocate(zcontrred(spi, spi, spi, spi), za(szo, spi), zb(szo, spi), &
           &         cmat(spi, spi), zcontrred2(spi, spi, spi, spi))

            do ii = 1, spi
                do jj = 1, spi
                  !  za = transpose(z1(:,posi,posi(ii)))
                    za = z1(:,posi,posi(ii))
                    zb = z2(:,posi,posi(jj))
                   ! cmat = matmul(za,zb)
                    call dgemm('t','n', spi, spi, szo, 1.0_dp/8.0_dp, za, &
                   &           szo, zb, szo, 0.0_dp, cmat, spi)
                    zcontrred(:,:,jj,ii) = cmat
                     za = z2(:,posi,posi(ii))
                              !  za = transpose(z1(:,posj,posi(ii)))
                                zb = z1(:,posi,posi(jj))
                               ! cmat = matmul(za,zb)
                                call dgemm('t','n', spi, spi, szo, 1.0_dp/8.0_dp, za, &
                              &           szo, zb, szo, 0.0_dp, cmat, spi)

                                zcontrred2(jj,ii,:,:) = cmat
                enddo
            enddo

            dx1=>dx(:,:,:,i,i)
            dy1=>dy(:,:,:,i,i)
            dz1=>dz(:,:,:,i,i)
            dx2=>dx(:,:,:,i,i)
            dy2=>dy(:,:,:,i,i)
            dz2=>dz(:,:,:,i,i)
!            zcontrred=zcontrred/8.0
!            zcontrred2=zcontrred2/8.0

            call tot_integral_ijkr_pzero(nq, l,m,n,group_start, group_count, p0matrix, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, i, i, &
                                    zcontrred, zcontrred2,  cutoffz, cutoffmd,f)

            tsi = tsi + f * e12(:, i, i) * e12(:, i, i)
            count=count+1
            deallocate(posI,za,zb,cmat)
            deallocate(zcontrred, zcontrred2)


        end do
        !$OMP END parallel DO
    print*, tsi(1)

    end subroutine tot_integration

    subroutine elastic_integration(ncap,px,py,pz,l,m,n,p0matrix,dx,dy,dz,z,group_start,group_count,group, &
            cutoffz,cutoffmd, cutoffcentre,q,e12,tsi)


        use omp_lib
        implicit none

              INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(kind=ikind), INTENT(IN) :: ncap
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:) :: l,m,n,group,group_start,group_count

        REAL(kind=dp), intent(in), dimension(:,:,:,:,:),target :: dx,dy,dz
        real(kind=dp),dimension(:,:,:), pointer :: dx1,dy1,dz1,dx2,dy2,dz2
        REAL(kind=dp), intent(in), dimension(:,:,:,:) :: p0matrix
        REAL(kind=dp), intent(in), dimension(:,:,:) ::e12
        real(kind=dp), intent(in), dimension(:,:)::z

        REAL(kind=dp), intent(in), dimension(:,:) :: px,py,pz
        REAL(kind=dp), intent(in), dimension(:) :: q
        REAL(kind=dp), intent(in) :: cutoffz, cutoffmd,cutoffcentre


        REAL(kind=dp), dimension(size(q)) :: f
        integer(kind=ikind), dimension(:,:), allocatable :: posits
        real(kind=dp),dimension(:,:), allocatable :: za,zb
        integer(kind=ikind),dimension(:), allocatable ::posi,posj,posk,posr
        REAL(kind=dp), intent(out), dimension(size(q)) :: tsi
        real(kind=dp) :: hx,hy,hz,h
        integer(kind=ikind),dimension(size(l)) :: ll
        integer(kind=ikind) :: nq,i,j,k,r,count,ng,ii,jj
        integer(kind=ikind) :: spi, spj, spk, spr,nt

        ng=maxval(group)
        ll=l+m+n
        allocate(posits(ng,maxval(group_count)))
        posits=1
        count=1
        do i=1,ncap
            count=1
            do j=group_start(i),group_start(i)+group_count(i)-1
               posits(i,count)=j
                count=count+1
            end do
        end do
        print*,'posits created'
        print*,size(group_start),sum(group_count)




        nq= size(q)

        print*,OMP_get_num_threads()
        call omp_set_num_threads(16)
        print*,OMP_get_num_threads()
        !First big loop

        tsi=0.0_dp
       !$OMP PARALLEL do private(posI,posK,posJ,posR,spi,spj,spk,spr,za,zb), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,j,k,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz, posits,cutoffmd,group_count,group_start,z) REDUCTION(+:tsi)

        do i=1,ncap
            do j=i+1,ncap
                do k=i+1,ncap
                    do r=k+1,ncap
                        hx = px(k, r) - px(i, j)
                        hy = py(k, r) - py(i, j)
                        hz = pz(k, r) - pz(i, j)
                        h = (hx * hx + hy * hy + hz * hz)**0.5

                        allocate(posI(size(posits(i,:group_count(i)))),posJ(size(posits(j,:group_count(j)))) &
                                ,posK(size(posits(k,:group_count(k)))),posR(size(posits(r,:group_count(r)))))


                        posI = posits(i,:group_count(i))
                        posJ = posits(j,:group_count(j))
                        posK = posits(k,:group_count(k))
                        posR = posits(r,:group_count(r))

                        spi = size(posI)
                        spj = size(posJ)
                        spk = size(posK)
                        spr = size(posR)

                        allocate(zb(spi, spj), &
                       &         za(spk,spr))

                         za=z(posI,posJ)
                         zb=z(posK,posR)

                         deallocate(posI,posJ,posK,posR)
                         dx1=>dx(:,:,:,j,i)
                         dy1=>dy(:,:,:,j,i)
                         dz1=>dz(:,:,:,j,i)
                         dx2=>dx(:,:,:,r,k)
                         dy2=>dy(:,:,:,r,k)
                         dz2=>dz(:,:,:,r,k)

                        if (h < cutoffcentre) then

                            call integral_ijkr_pzero(nq, l,m,n,group_start, group_count, p0matrix, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, k, r, &
                                    za,zb,  cutoffz, cutoffmd,f)



                        else

                            call integral_k_ijkr(q,l,m,n,group_start, group_count, hx, hy, hz, h,dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, k, r, &
                                    za,zb,  cutoffz, cutoffmd,f)





                        end if
                        tsi = tsi + 8.000 * f * e12(:, i, j) * e12(:, k, r)

                        count=count+1
                        deallocate(za,zb)
               !        deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)


                    end do
                end do
            end do
        end do

       !$OMP END parallel DO


        print*,'la rata ', tsi(1)

        !$OMP PARALLEL do private(posI,posJ,posR,spi,spj,spk,spr,za,zb), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,j,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz,posits, cutoffmd,group_count,group_start,z) REDUCTION(+:tsi)
        do i=1,ncap
            do j=i+1,ncap
                do r=i+1,ncap

                    hx = px(i, r) - px(i, j)
                    hy = py(i, r) - py(i, j)
                    hz = pz(i, r) - pz(i, j)
                    h = sqrt((hx * hx + hy * hy + hz * hz))
                    allocate(posI(size(posits(i,:group_count(i)))),posJ(size(posits(j,:group_count(j)))) &
                               ,posR(size(posits(r,:group_count(r)))))


                        posI = posits(i,:group_count(i))
                        posJ = posits(j,:group_count(j))

                        posR = posits(r,:group_count(r))

                        spi = size(posI)
                        spj = size(posJ)
                        spr = size(posR)

                        allocate(za(spi, spj), zb(spi, spr))

                        za=z(posI,posJ)
                        zb=z(posI,posR)
                        deallocate(posI,posJ,posR)
                    dx1=>dx(:,:,:,j,i)
                    dy1=>dy(:,:,:,j,i)
                    dz1=>dz(:,:,:,j,i)
                    dx2=>dx(:,:,:,r,i)
                    dy2=>dy(:,:,:,r,i)
                    dz2=>dz(:,:,:,r,i)


                    if (h < cutoffcentre) then
                        call integral_ijkr_pzero(nq, l,m,n,group_start, group_count, p0matrix, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, i, r, &
                               za,zb,  cutoffz, cutoffmd,f)
                    else

                        call  integral_k_ijkr(q, l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, i, r, &
                                za,zb,  cutoffz, cutoffmd,f)



                    end if
                    tsi = tsi + 4.000 * f * e12(:, i, j) * e12(:, i, r)
                    count=count+1

                    deallocate(za,zb)
                  !   deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

                end do
            end do

        end do
    !$OMP END parallel DO
        print*,'up to second'
      !$OMP PARALLEL do private(posI,posK,posR,spi,spj,spk,spr,za,zb), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,k,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz,posits, cutoffmd,group_count,group_start,z) REDUCTION(+:tsi)
        do i=1,ncap
            do k=1,ncap
                do r=k+1,ncap
                    hx = px(k, r) - px(i, i)
                    hy = py(k, r) - py(i, i)
                    hz = pz(k, r) - pz(i, i)
                    h = sqrt((hx * hx + hy * hy + hz * hz))

                    allocate(posI(size(posits(i,:group_count(i)))) &
                                ,posK(size(posits(k,:group_count(k)))),posR(size(posits(r,:group_count(r)))))


                        posI = posits(i,:group_count(i))

                        posK = posits(k,:group_count(k))
                        posR = posits(r,:group_count(r))

                    spi = size(posI)
                    spk = size(posK)
                    spr = size(posR)

                    allocate(za(spi,spi), zb(spk,spr))

                    za=z(posi,posi)
                    zb=z(posk,posr)

                    dx1=>dx(:,:,:,i,i)
                    dy1=>dy(:,:,:,i,i)
                    dz1=>dz(:,:,:,i,i)
                    dx2=>dx(:,:,:,r,k)
                    dy2=>dy(:,:,:,r,k)
                    dz2=>dz(:,:,:,r,k)
!                    zcontrred=zcontrred
!                    zcontrred2=zcontrred2

                    if (h < cutoffcentre) then
                        call integral_ijkr_pzero(nq,l,m,n,group_start, group_count, p0matrix, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, k, r, &
                                    za,zb, cutoffz, cutoffmd, f)
                    else

                        call  integral_k_ijkr(q, l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, k, r, &
                                    za,zb,  cutoffz, cutoffmd, f)



                    end if
                    tsi = tsi+ 4.000 * f * e12(:, i, i) * e12(:, k, r)
                    count=count+1
                    deallocate(posI,posK,posR,za,zb)

               !     deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

                end do
            end do
        end do

        !$OMP END parallel DO

          !$OMP PARALLEL do private(posI,posK,spi,spk,za,zb), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,k,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
          !$OMP& shared( cutoffz, cutoffmd,posits,group_count,group_start) REDUCTION(+:tsi)
        do i=1,ncap
            do k=i+1,ncap

                hx = px(k, k) - px(i, i)
                hy = py(k, k) - py(i, i)
                hz = pz(k, k) - pz(i, i)
                h = sqrt((hx * hx + hy * hy + hz * hz))

                allocate(posI(size(posits(i,:group_count(i)))), &
                        posK(size(posits(k,:group_count(k)))))


                posI = posits(i,:group_count(i))

                posK = posits(k,:group_count(k))


                spi = size(posI)
                spk = size(posK)

                allocate(za(spk,spk), zb(spk,spk))

               za=z(posi,posi)
                zb=z(posk,posk)

                dx1=>dx(:,:,:,i,i)
                dy1=>dy(:,:,:,i,i)
                dz1=>dz(:,:,:,i,i)
                dx2=>dx(:,:,:,k,k)
                dy2=>dy(:,:,:,k,k)
                dz2=>dz(:,:,:,k,k)


                if (h < cutoffcentre) then
                    call integral_ijkr_pzero(nq, l,m,n,group_start, group_count, p0matrix, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, k, k, &
                                    za,zb,  cutoffz, cutoffmd, f)
                else

                    call  integral_k_ijkr(q,l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, k, k, &
                                    za,zb,  cutoffz, cutoffmd, f)
                end if
                tsi = tsi+ 2.000 * f * e12(:, i, i) * e12(:, k, k)
                deallocate(posI,posK,za,zb)
                count=count+1
               ! deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

            end do
        end do


        !$OMP END parallel DO
        !$OMP PARALLEL do private(posI,spi,za,zb), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,ll, p0matrix), &
        !$OMP& shared( cutoffz, cutoffmd,posits,group_count,group_start) REDUCTION(+:tsi)
        do i=1,ncap
                              allocate(posI(size(posits(i,:group_count(i)))))


                        posI = posits(i,:group_count(i))



            spi = size(posI)

            allocate( za(spi, spi), zb(spi, spi))

            za=z(posi,posi)
            zb=z(posi,posi)
            dx1=>dx(:,:,:,i,i)
            dy1=>dy(:,:,:,i,i)
            dz1=>dz(:,:,:,i,i)
            dx2=>dx(:,:,:,i,i)
            dy2=>dy(:,:,:,i,i)
            dz2=>dz(:,:,:,i,i)
!            zcontrred=zcontrred/8.0
!            zcontrred2=zcontrred2/8.0

            call integral_ijkr_pzero(nq, l,m,n,group_start, group_count, p0matrix, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, i, i, &
                                    za,zb,  cutoffz, cutoffmd,f)

            tsi = tsi + f * e12(:, i, i) * e12(:, i, i)
            count=count+1
            deallocate(posI,za,zb)



        end do
        !$OMP END parallel DO

    print*,count
    end subroutine elastic_integration

    SUBROUTINE set_P0(P0, lmax4, q)

        ! P0    (Nq,4lmax+1,4lmax+1,4lmax+1) will store the values of <qx^L*qy^M*qz^N>
        !   for all possible combinations of angular momenta of the 4 GTOs (L,M,N)
              INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        REAL(KIND=dp), INTENT(OUT), DIMENSION(:,:,:,:)  :: P0
        ! the maximum value of the orbital angular momentum of any GTO times 4
        INTEGER(KIND=ikind), INTENT(IN)                 :: lmax4
        ! the radial component of the scattering vector
        REAL(KIND=dp), INTENT(IN), DIMENSION(:)         :: q


        ! loop variables
        INTEGER(KIND=ikind) :: i, j, k


        ! loop for all combinations of angular momenta
        do k = 0, lmax4
            do j = 0, lmax4
                do i = 0, lmax4
                    ! Fill the current LMN
                    CALL P_LMN(P0, i, j, k, q)
                end do
            end do
        end do

    END SUBROUTINE set_P0

!--------------------------------------------------------------------
    ! populate just the current part of the matrix
    SUBROUTINE P_LMN(P0, L, M, N, q)

        ! the three values of the angular momentum
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(KIND=ikind), INTENT(IN) :: L, M, N
        ! The value of <qx^L*qy^M*qz^N>
        REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:) :: P0
        ! The radial part of the scattering vector
        REAL(KIND=dp), INTENT(IN), DIMENSION(:)     :: q


        ! For ordering
        INTEGER(KIND=ikind), DIMENSION(1:3)     :: A
        ! helper variable
        !INTEGER(KIND=ikind) :: i

        ! Order L, M and N
        A(1) = L
        A(2) = M
        A(3) = N
        CALL Bubble_Sort(A)

        ! These are analytical solutions to the integral
        ! They have been obtained in Matematica by Andres Moreno
        ! They have been extensively debugged in the MATLAB version of the code
        ! however one should be careful with the precision as there is nasty
        ! divisions and exponentiations. Needs to be retested if it works well
        ! with the current data kinds.
        if (mod(L,2)/=0 .or. mod(M,2)/=0 .or. mod(N,2)/=0) then
            P0(:,L+1,M+1,N+1) = 0.0_dp
        elseif (sum(A)==0) then
            P0(:,L+1,M+1,N+1) = 1.0_dp
        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==2) then
            P0(:,L+1,M+1,N+1)=-q**2/3.
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==2) then
            P0(:,L+1,M+1,N+1)=(q**4./15._dp)
        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==2) then
            P0(:,L+1,M+1,N+1)=-(q**6./105._dp)
        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==4) then
            P0(:,L+1,M+1,N+1)=(q**4./5._dp)
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==4) then
            P0(:,L+1,M+1,N+1)=-(q**6./35._dp)
        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==4) then
            P0(:,L+1,M+1,N+1)=(q**8./315._dp)
        elseif (A(1)==2 .and. A(2)==4 .and. A(3)==4) then
            P0(:,L+1,M+1,N+1)=-(q**10./1155._dp)
        elseif (A(1)==4 .and. A(2)==4 .and. A(3)==4) then
            P0(:,L+1,M+1,N+1)=(q**12./5005._dp)
        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=-(q**6./7._dp)
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=(q**8./63._dp)
        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=-(q**10./693._dp)
        elseif (A(1)==2 .and. A(2)==4 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=(q**12./3003._dp)
        elseif (A(1)==0 .and. A(2)==4 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=-(q**10./231._dp)
        elseif (A(1)==4 .and. A(2)==4 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=-(q**14./15015._dp)
        elseif (A(1)==2 .and. A(2)==6 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=-(q**14./9009._dp)
        elseif (A(1)==4 .and. A(2)==6 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=(q**16./51051._dp)
        elseif (A(1)==6 .and. A(2)==6 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=-5*(q**18./969969._dp)
        elseif (A(1)==0 .and. A(2)==6 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=5*(q**12./3003._dp)
        elseif (A(1)==0 .and. A(2)==4 .and. A(3)==4) then
            P0(:,L+1,M+1,N+1)=(q**8./105._dp)
        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(q**8./9._dp)
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=-(q**10./99._dp)
        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(q**12./1287._dp)


        elseif (A(1)==0 .and. A(2)==4 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(q**12./429._dp)
        elseif (A(1)==2 .and. A(2)==4 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=-(q**14./6435._dp)
        elseif (A(1)==4 .and. A(2)==4 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(q**16./36465._dp)


        elseif (A(1)==0 .and. A(2)==6 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=-(q**14./1287._dp)
        elseif (A(1)==2 .and. A(2)==6 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(q**16./21879._dp)
        elseif (A(1)==4 .and. A(2)==6 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=-(q**18./138567._dp)
        elseif (A(1)==6 .and. A(2)==6 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=5*(q**20./2909907._dp)


        elseif (A(1)==0 .and. A(2)==8 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=7*(q**16./21879._dp)
        elseif (A(1)==2 .and. A(2)==8 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=-7*(q**18./415701._dp)
        elseif (A(1)==4 .and. A(2)==8 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(q**20./415701._dp)
        elseif (A(1)==6 .and. A(2)==8 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=-5*(q**22./9561123._dp)
        elseif (A(1)==8 .and. A(2)==8 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=7*(q**24./47805615._dp)

        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=-(q**10/11._dp)
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(q**12./143._dp)
        elseif (A(1)==0 .and. A(2)==4 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=-(q**14./715._dp)
        elseif (A(1)==0 .and. A(2)==6 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(q**16./2431._dp)
        elseif (A(1)==0 .and. A(2)==8 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=-7*(q**18./46189._dp)
        elseif (A(1)==0 .and. A(2)==10 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=3*(q**20./46189._dp)

        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(q**12/13._dp)
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=-(q**14./195._dp)
        elseif (A(1)==0 .and. A(2)==4 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(q**16./1105._dp)
        elseif (A(1)==0 .and. A(2)==6 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=-(q**18./4199._dp)
        elseif (A(1)==0 .and. A(2)==8 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(q**20./12597._dp)
        elseif (A(1)==0 .and. A(2)==10 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=-3*(q**22./96577._dp)
        elseif (A(1)==0 .and. A(2)==12 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=33*(q**24./2414425._dp)

        elseif (A(1)==0 .and. A(2)==10 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=3*(q**20/46189._dp)
        elseif (A(1)==2 .and. A(2)==10.and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=-3*(q**22./1062347._dp)
        elseif (A(1)==4 .and. A(2)==10 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=9*(q**24./26558675._dp)
        elseif (A(1)==6 .and. A(2)==10 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=-(q**26./15935205._dp)
        elseif (A(1)==8 .and. A(2)==10 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=7*(q**28./462120945._dp)
        elseif (A(1)==10 .and. A(2)==10 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=-21*(q**30./4775249765._dp)

        elseif (A(1)==2 .and. A(2)==12 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=-11*(q**26/21729825._dp)
        elseif (A(1)==4 .and. A(2)==12.and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=11*(q**28./210054975._dp)
        elseif (A(1)==6 .and. A(2)==12 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=-11*(q**30./1302340845._dp)
        elseif (A(1)==8 .and. A(2)==12 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=7*(q**32./3907022535._dp)
        elseif (A(1)==10 .and. A(2)==12 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=-(q**34./2170568075._dp)
        elseif (A(1)==12 .and. A(2)==12 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=11*(q**36./80311018775._dp)


        elseif (A(1)==2 .and. A(2)==10 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=3*(q**24/2414425._dp)
        elseif (A(1)==4 .and. A(2)==10.and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=-(q**26./7243275._dp)
        elseif (A(1)==6 .and. A(2)==10 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(q**28./42010995._dp)
        elseif (A(1)==8 .and. A(2)==10 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=-7*(q**30./1302340845._dp)
        elseif (A(1)==10 .and. A(2)==10 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=7*(q**32./4775249765._dp)

        elseif (A(1)==2 .and. A(2)==8 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=-(q**22/289731._dp)
        elseif (A(1)==4 .and. A(2)==8.and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(q**24./2414425._dp)
        elseif (A(1)==6 .and. A(2)==8 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=-(q**26./13037895._dp)
        elseif (A(1)==8 .and. A(2)==8 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=7*(q**28./378098955._dp)

        elseif (A(1)==2 .and. A(2)==6 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(q**20/88179._dp)
        elseif (A(1)==4 .and. A(2)==6 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=-(q**22./676039._dp)
        elseif (A(1)==6 .and. A(2)==6 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(q**24./3380195._dp)

        elseif (A(1)==2 .and. A(2)==4 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=-(q**18/20995._dp)
        elseif (A(1)==4 .and. A(2)==4 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(q**20./146965._dp)

        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(q**16./3315._dp)

        elseif (A(1)==2 .and. A(2)==8 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(q**20./138567)
        elseif (A(1)==4 .and. A(2)==8.and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=-(q**22./1062347._dp)
        elseif (A(1)==6 .and. A(2)==8 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(q**24./5311735)
        elseif (A(1)==8 .and. A(2)==8 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=-7*(q**26./143416845._dp)

        elseif (A(1)==2 .and. A(2)==6 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=-(q**18/46189._dp)
        elseif (A(1)==4 .and. A(2)==6.and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(q**20./323323._dp)
        elseif (A(1)==6 .and. A(2)==6 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=-5*(q**22./7436429._dp)

        elseif (A(1)==2 .and. A(2)==4 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(q**16/12155._dp)
        elseif (A(1)==4 .and. A(2)==4.and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=-3*(q**18./230945._dp)

        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=-(q**14./2145._dp)
        else
            print*, "CASE NOT PROGRAMED YET", L, M, N
            P0(:,L+1,M+1,N+1) = 0.0_dp
            STOP

        end if

    END SUBROUTINE P_LMN

SUBROUTINE BesselDeriv(BD, LL, MM,NN,a,b,c,LLmax)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        ! The three nested loops give the formula
        ! in the form of coefficients multiplying the h functions
        REAL(kind=dp), INTENT(out), DIMENSION(LLmax+1)  :: BD
        INTEGER(kind=ikind), INTENT(in)                 :: LL, MM, NN, LLmax
        REAL(kind=dp), INTENT(in), DIMENSION(LLmax+1,LLmax+1)  :: a, b, c

        ! loop and temp variables
        INTEGER(kind=ikind) :: ii, jj, kk, hOrder, temp, ceil
        REAL(kind=dp)       :: C1, C2, C3, Ct2, Ct3
        ! set this to 0 initially and accumulate
        BD = 0
        !TempBesselPre=zeros(LLmax+1,1);
        do ii = 0, LL
            C1=a(LL+1,ii+1)
            if (abs(C1)<1.0E-30) cycle

            do jj = 0, MM
                Ct2=b(MM+1,jj+1)
                if (abs(Ct2)<1.0E-30) cycle
                C2 = C1 * Ct2

                do kk = 0, NN
                    Ct3=c(NN+1,kk+1)
                    if (abs(Ct3)<1.0E-30) cycle
                    C3 = C2 * Ct3
                    ! hOrder = ceiling((LL+MM+NN-ii-jj-kk)/2.0_dp)+ii+jj+kk
                    temp = LL+MM+NN-ii-jj-kk
                    ceil = temp/2_ikind + mod(temp, 2_ikind) ! integer division with rounding up
                    hOrder = ceil+ii+jj+kk
                    BD(hOrder+1)=BD(hOrder+1) + C3
                end do
            end do
        end do
    END SUBROUTINE




    SUBROUTINE integral_ijkr_pzero(nq,l,m,n,gs,gc,p0mat,dx1,dy1,dz1,dx2,dy2,dz2,gi,gj,gk,gr,za,zb, &
             cutoff1,cutoff2,f)



        implicit none

              INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

        ! definition of input

        INTEGER(kind=ikind), INTENT(IN)                       :: nq, gi, gj, gk, gr
        REAL(kind=dp), INTENT(IN)                             :: cutoff1, cutoff2
        REAL(kind=dp), INTENT(IN), DIMENSION(:,:,:,:)               :: p0mat
        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1,dx2,dy2,dz2
         real(kind=dp), intent(in),  dimension(:,:) :: za,zb

        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:)         :: gs,gc,l,m,n
        ! definition of output
        REAL(kind=dp), INTENT(OUT), DIMENSION(nq)             :: f
        ! definition of loop indices
        INTEGER(kind=ikind)                                   :: h1
        INTEGER(kind=ikind)                                   :: ll, mm, nn, llp, mmp, nnp
        ! definition of internal variables
        INTEGER(kind=ikind)                                   :: i,j,k,r, posi, posj, posk, posr
        REAL(kind=dp)                                         :: z1,ztot
        real(dp), external :: ddot
        REAL(kind=dp)                                         :: mdl, mdm, mdn, mdlp, mdmp,mdnp
        real(kind=dp) :: prod6, prod5,prod4,prod3,prod2,prod1
!        REAL(kind=dp), DIMENSION(size(Z1(:,1,1)))                          :: zij1, zkr1
!        REAL(kind=dp), DIMENSION(size(Z2(:,1,1)))                          :: zij2, zkr2


        posI=1
       ! posI=apos(i)

        f=0.0
        ! loop through all possible ways to get total angular momentum lmax1
        do i = gs(gi), gs(gi) + gc(gi) - 1
            posj=1

            do j = gs(gj), gs(gj) + gc(gj) - 1
                Z1=Za(posI, posJ)
                posk=1
                do k = gs(gk), gs(gk) + gc(gk) - 1
                posr=1
                    do r = gs(gr), gs(gr) + gc(gr) - 1
                        ztot=Z1*Zb(posk,posr)
!

                        ! continue only if larger
                        if (abs(ztot) < cutoff1) then
                            posr=posr+1
                            cycle
                        end if
                     do ll = 0, l(i)+l(j)
                            MDL = Dx1(ll+1,l(i)+1,l(j)+1)
                            if (abs(MDL)<1.0e-30) cycle
                            prod1 = MDL * ztot
                            ! MD coeff 2
                            do mm = 0, m(i)+m(j)
                                MDM = Dy1(mm+1,m(i)+1,m(j)+1)
                                if (abs(MDM)<1.0e-30) cycle
                                prod2 = MDM * prod1
                                ! MD coeff 3
                                do nn =0, n(i)+n(j)
                                    H1=(-1.0)**(ll+mm+nn)
                                    MDN=Dz1(nn+1,n(i)+1,n(j)+1)
                                    if (abs(MDN)<1.0e-30) cycle ! check if MD coeff is 0
                                    prod3 = MDN * H1  * prod2
                                    ! MD coeff 4
                                    do llp = 0, l(k)+l(r)
                                        MDLp=Dx2(llp+1,l(k)+1,l(r)+1)
                                        if (abs(MDLp)<1.0e-30) cycle ! check if MD coeff is 0
                                        prod4 = MDLp * prod3
                                        ! MD coeff 5
                                        do mmp = 0, m(k)+m(r)
                                            MDMp = Dy2(mmp+1,m(k)+1,m(r)+1)
                                            if (abs(MDMp)<1.0e-30) cycle ! check if MD coeff is 0
                                            prod5 = MDMp * prod4
                                            ! MD coeff 6
                                            do nnp = 0, n(k)+n(r)
                                                MDNp = Dz2(nnp+1,n(k)+1,n(r)+1)
                                                prod6 = MDNp * prod5
                                                ! cutoff after MD
                                                if (abs(prod6)>cutoff2) then
                                                    ! add the contribution to the total
                                                    F = F + prod6 * P0mat(:,ll+llp+1,mm+mmp+1,nn+nnp+1)
                                                    !Int=Int+prodD.*f;
                                                end if
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do

                    posr=posr+1
                    end do
                posk=posk+1
                end do
            posj=posj+1
            end do
        posi=posi+1
        end do

    END SUBROUTINE




SUBROUTINE tot_integral_ijkr_pzero(nq,l,m,n,gs,gc,p0mat,dx1,dy1,dz1,dx2,dy2,dz2,gi,gj,gk,gr,zcontr,zcontr2, &
             cutoff1,cutoff2,f)



        implicit none

        ! definition of input
              INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(kind=ikind), INTENT(IN)                       :: nq, gi, gj, gk, gr
        REAL(kind=dp), INTENT(IN)                             :: cutoff1, cutoff2
        REAL(kind=dp), INTENT(IN), DIMENSION(:,:,:,:)               :: p0mat
        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1,dx2,dy2,dz2
        REAL(kind=dp), INTENT(IN), DIMENSION(:,:,:,:)               :: zcontr,zcontr2


        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:)         :: gs,gc,l,m,n
        ! definition of output
        REAL(kind=dp), INTENT(OUT), DIMENSION(nq)             :: f
        ! definition of loop indices
        INTEGER(kind=ikind)                                   :: h1
        INTEGER(kind=ikind)                                   :: ll, mm, nn, llp, mmp, nnp
        ! definition of internal variables
        INTEGER(kind=ikind)                                   :: i,j,k,r, posi, posj, posk, posr
        REAL(kind=dp)                                         :: ztot
        real(dp), external :: ddot
        REAL(kind=dp)                                         :: mdl, mdm, mdn, mdlp, mdmp,mdnp
        real(kind=dp) :: prod6, prod5,prod4,prod3,prod2,prod1
!        REAL(kind=dp), DIMENSION(size(Z1(:,1,1)))                          :: zij1, zkr1
!        REAL(kind=dp), DIMENSION(size(Z2(:,1,1)))                          :: zij2, zkr2


        posI=1
       ! posI=apos(i)

        f=0.0
        ! loop through all possible ways to get total angular momentum lmax1
        do i = gs(gi), gs(gi) + gc(gi) - 1
            posj=1
            do j = gs(gj), gs(gj) + gc(gj) - 1
                posk=1
                do k = gs(gk), gs(gk) + gc(gk) - 1
                posr=1
                    do r = gs(gr), gs(gr) + gc(gr) - 1
                        ztot=zcontr(posJ,posK,posR,posI)+zcontr2(posR,posI,posJ,posK)
!

                        ! continue only if larger
                        if (abs(ztot) < cutoff1) then
                            posr=posr+1
                            cycle
                        end if
                     do ll = 0, l(i)+l(j)
                            MDL = Dx1(ll+1,l(i)+1,l(j)+1)
                            if (abs(MDL)<1.0e-30) cycle
                            prod1 = MDL * ztot
                            ! MD coeff 2
                            do mm = 0, m(i)+m(j)
                                MDM = Dy1(mm+1,m(i)+1,m(j)+1)
                                if (abs(MDM)<1.0e-30) cycle
                                prod2 = MDM * prod1
                                ! MD coeff 3
                                do nn =0, n(i)+n(j)
                                    H1=(-1.0)**(ll+mm+nn)
                                    MDN=Dz1(nn+1,n(i)+1,n(j)+1)
                                    if (abs(MDN)<1.0e-30) cycle ! check if MD coeff is 0
                                    prod3 = MDN * H1  * prod2
                                    ! MD coeff 4
                                    do llp = 0, l(k)+l(r)
                                        MDLp=Dx2(llp+1,l(k)+1,l(r)+1)
                                        if (abs(MDLp)<1.0e-30) cycle ! check if MD coeff is 0
                                        prod4 = MDLp * prod3
                                        ! MD coeff 5
                                        do mmp = 0, m(k)+m(r)
                                            MDMp = Dy2(mmp+1,m(k)+1,m(r)+1)
                                            if (abs(MDMp)<1.0e-30) cycle ! check if MD coeff is 0
                                            prod5 = MDMp * prod4
                                            ! MD coeff 6
                                            do nnp = 0, n(k)+n(r)
                                                MDNp = Dz2(nnp+1,n(k)+1,n(r)+1)
                                                prod6 = MDNp * prod5
                                                ! cutoff after MD
                                                if (abs(prod6)>cutoff2) then
                                                    ! add the contribution to the total
                                                    F = F + prod6 * P0mat(:,ll+llp+1,mm+mmp+1,nn+nnp+1)
                                                    !Int=Int+prodD.*f;
                                                end if
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do

                    posr=posr+1
                    end do
                posk=posk+1
                end do
            posj=posj+1
            end do
        posi=posi+1
        end do

    END SUBROUTINE

    subroutine tot_integral_k_ijkr(mu,l,m,n,group_start,group_count,hx,hy,hz,h,dx1,dy1,dz1,dx2,dy2,dz2, gi,gj,gk, gr,&
            zcontr, zcontr2, &
            cutoff1, cutoff2,f)



        implicit none


        !INTEGER, PARAMETER :: dp = selected_real_kind(2*precision(1.0_dp))
              INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=selected_int_kind(8)), intent(in)  :: gi,gj,gk,gr
        integer(kind=selected_int_kind(8)), dimension(:), intent(in) :: l,m,n,group_start,group_count
        real(kind=dp), intent(in)              :: cutoff1,cutoff2, hx, hy, hz, h
        real(kind=dp), intent(in), dimension(:,:,:,:) :: zcontr,zcontr2
        real(kind=dp), intent(in), dimension(:)   ::  mu

        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1,dx2,dy2,dz2
        !real(kind=dp), pointer,dimension(:,:,:) :: dxx,dyy,dzz, dxx2,dyy2,dzz2



        integer(kind=selected_int_kind(8)) :: ka, posi, posj, posk, posr, ra,i,j,k,r
        integer(kind=selected_int_kind(8)) ::  h1
        integer(kind=selected_int_kind(8)) ::  ll, mm, nn, llp, mmp, nnp,llmax

        real(kind=dp)  ::coeff,prodd,ztot,mdn, mdl, mdm, mdlp, mdmp,mdnp,z11,z22
        INTEGER(kind=ikind), parameter :: dim = 13
        real(kind=dp) :: prod1,prod2,prod3,prod4,prod5,prod6
        REAL(kind=dp), DIMENSION(dim,dim)           :: a, b, c
        REAL(kind=dp), DIMENSION(dim)               :: h_saved
        REAL(kind=dp), DIMENSION(dim, dim, dim, dim)  :: h_pre2
        REAL(kind=dp), DIMENSION(dim)               :: BD



      !  real(kind=dp), dimension(:), allocatable :: pmu, h_0, h_1, h_r, muoh,zij, zij2, zkr, zkr2

        real(kind=dp), intent(out), dimension(size(mu)) :: f
        real(kind=dp), allocatable,dimension(:) :: l1vec,l2vec,suml1l2


        LLmax = l(group_start(gi)) + m(group_start(gi)) + n(group_start(gi)) + &
                l(group_start(gj)) + m(group_start(gj)) + n(group_start(gj)) + &
                l(group_start(gk)) + m(group_start(gk)) + n(group_start(gk)) + &
                l(group_start(gr)) + m(group_start(gr)) + n(group_start(gr))


        if (LLmax + 1 > dim) then
            print*, "only s,p,d and f type GTOs are supported"
            stop
        end if



        call Hermite_like_coeffs(a, LLmax, Hx)
        call Hermite_like_coeffs(b, LLmax, Hy)
        call Hermite_like_coeffs(c, LLmax, Hz)




        bd=0.0_dp

        do k = 0, LLmax
            do j = 0, LLmax - k
                do i = 0, LLmax - k - j
                    call BesselDeriv(BD, i, j, k, a, b, c, LLmax)
                    h_pre2(:,i+1,j+1,k+1) = BD
                end do
            end do
        end do

       ! posi=apos(i)



        h_saved=0.0_dp

! loop through all possible ways to get total angular momentum lmax1


        posI=1


        ! loop through all possible ways to get total angular momentum lmax1
        do i = group_start(gi), group_start(gi) + group_count(gi) - 1
            posj=1
            do j = group_start(gj), group_start(gj) + group_count(gj) - 1
                posk=1
                do k = group_start(gk), group_start(gk) + group_count(gk) - 1
                    posr=1
                    do r = group_start(gr), group_start(gr) + group_count(gr) - 1

                        ztot=zcontr(posJ,posK,posR,posI)+zcontr2(posR,posI,posJ,posK)
!
                        ! continue only if larger
                        if (abs(ztot) < cutoff1) then
                                posr=posr+1
                                cycle
                            end if
                           do ll = 0, l(i)+l(j)
                            MDL = Dx1(ll+1,l(i)+1,l(j)+1)
                            if (abs(MDL)<1.0e-30) cycle
                            prod1 = MDL * ztot
                            ! MD coeff 2
                            do mm = 0, m(i)+m(j)
                                MDM = Dy1(mm+1,m(i)+1,m(j)+1)
                                if (abs(MDM)<1.0e-30) cycle
                                prod2 = MDM * prod1
                                ! MD coeff 3
                                do nn =0, n(i)+n(j)
                                    H1=(-1.0)**(ll+mm+nn)
                                    MDN=Dz1(nn+1,n(i)+1,n(j)+1)
                                    if (abs(MDN)<1.0e-30) cycle ! check if MD coeff is 0
                                    prod3 = MDN * H1  * prod2
                                    ! MD coeff 4
                                    do llp = 0, l(k)+l(r)
                                        MDLp=Dx2(llp+1,l(k)+1,l(r)+1)
                                        if (abs(MDLp)<1.0e-30) cycle ! check if MD coeff is 0
                                        prod4 = MDLp * prod3
                                        ! MD coeff 5
                                        do mmp = 0, m(k)+m(r)
                                            MDMp = Dy2(mmp+1,m(k)+1,m(r)+1)
                                            if (abs(MDMp)<1.0e-30) cycle ! check if MD coeff is 0
                                            prod5 = MDMp * prod4
                                            ! MD coeff 6
                                            do nnp = 0, n(k)+n(r)
                                                MDNp = Dz2(nnp+1,n(k)+1,n(r)+1)
                                                prod6 = MDNp * prod5
                                                ! cutoff after MD
                                                if (abs(prod6)>cutoff2) then
                                                    ! add the contribution to the total
                                                    h_saved = h_saved + h_pre2(:,ll+llp+1,mm+mmp+1,nn+nnp+1)*prod6

                                                end if
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do

                    posr=posr+1
                    end do
                    posk=posk+1
                end do
                    posj=posj+1
            end do
        posi=posi+1
        end do

        CALL BesselSum(F, mu, H, LLmax, h_saved)

    end subroutine


    subroutine integral_k_ijkr(mu,l,m,n,group_start,group_count,hx,hy,hz,h,dx1,dy1,dz1,dx2,dy2,dz2, gi,gj,gk, gr,&
            za,zb,&
            cutoff1, cutoff2,f)



        implicit none

              INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        !INTEGER, PARAMETER :: dp = selected_real_kind(2*precision(1.0_dp))
        integer(kind=selected_int_kind(8)), intent(in)  :: gi,gj,gk,gr
        integer(kind=selected_int_kind(8)), dimension(:), intent(in) :: l,m,n,group_start,group_count
        real(kind=dp), intent(in)              :: cutoff1,cutoff2, hx, hy, hz, h
        real(kind=dp), intent(in), dimension(:,:) :: za,zb
        real(kind=dp), intent(in), dimension(:)   ::  mu

        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1,dx2,dy2,dz2
        !real(kind=dp), pointer,dimension(:,:,:) :: dxx,dyy,dzz, dxx2,dyy2,dzz2



        integer(kind=selected_int_kind(8)) :: ka, posi, posj, posk, posr, ra,i,j,k,r
        integer(kind=selected_int_kind(8)) ::  h1
        integer(kind=selected_int_kind(8)) ::  ll, mm, nn, llp, mmp, nnp,llmax

        real(kind=dp)  ::coeff,prodd,ztot,mdn, mdl, mdm, mdlp, mdmp,mdnp,z1
        INTEGER(kind=ikind), parameter :: dim = 13
        real(kind=dp) :: prod1,prod2,prod3,prod4,prod5,prod6
        REAL(kind=dp), DIMENSION(dim,dim)           :: a, b, c
        REAL(kind=dp), DIMENSION(dim)               :: h_saved
        REAL(kind=dp), DIMENSION(dim, dim, dim, dim)  :: h_pre2
        REAL(kind=dp), DIMENSION(dim)               :: BD



      !  real(kind=dp), dimension(:), allocatable :: pmu, h_0, h_1, h_r, muoh,zij, zij2, zkr, zkr2

        real(kind=dp), intent(out), dimension(size(mu)) :: f
        real(kind=dp), allocatable,dimension(:) :: l1vec,l2vec,suml1l2


        LLmax = l(group_start(gi)) + m(group_start(gi)) + n(group_start(gi)) + &
                l(group_start(gj)) + m(group_start(gj)) + n(group_start(gj)) + &
                l(group_start(gk)) + m(group_start(gk)) + n(group_start(gk)) + &
                l(group_start(gr)) + m(group_start(gr)) + n(group_start(gr))


        if (LLmax + 1 > dim) then
            print*, "only s,p,d and f type GTOs are supported"
            stop
        end if



        call Hermite_like_coeffs(a, LLmax, Hx)
        call Hermite_like_coeffs(b, LLmax, Hy)
        call Hermite_like_coeffs(c, LLmax, Hz)




        bd=0.0_dp

        do k = 0, LLmax
            do j = 0, LLmax - k
                do i = 0, LLmax - k - j
                    call BesselDeriv(BD, i, j, k, a, b, c, LLmax)
                    h_pre2(:,i+1,j+1,k+1) = BD
                end do
            end do
        end do

       ! posi=apos(i)



        h_saved=0.0_dp

! loop through all possible ways to get total angular momentum lmax1


        posI=1


        ! loop through all possible ways to get total angular momentum lmax1
        do i = group_start(gi), group_start(gi) + group_count(gi) - 1
            posj=1
            do j = group_start(gj), group_start(gj) + group_count(gj) - 1
                posk=1
                do k = group_start(gk), group_start(gk) + group_count(gk) - 1
                    z1=Za(posi,posj)
                    posr=1
                    do r = group_start(gr), group_start(gr) + group_count(gr) - 1

                        ztot=Zb(posK,posR)*Z1;
!
                        ! continue only if larger
                        if (abs(ztot) < cutoff1) then
                                posr=posr+1
                                cycle
                            end if
                           do ll = 0, l(i)+l(j)
                            MDL = Dx1(ll+1,l(i)+1,l(j)+1)
                            if (abs(MDL)<1.0e-30) cycle
                            prod1 = MDL * ztot
                            ! MD coeff 2
                            do mm = 0, m(i)+m(j)
                                MDM = Dy1(mm+1,m(i)+1,m(j)+1)
                                if (abs(MDM)<1.0e-30) cycle
                                prod2 = MDM * prod1
                                ! MD coeff 3
                                do nn =0, n(i)+n(j)
                                    H1=(-1.0)**(ll+mm+nn)
                                    MDN=Dz1(nn+1,n(i)+1,n(j)+1)
                                    if (abs(MDN)<1.0e-30) cycle ! check if MD coeff is 0
                                    prod3 = MDN * H1  * prod2
                                    ! MD coeff 4
                                    do llp = 0, l(k)+l(r)
                                        MDLp=Dx2(llp+1,l(k)+1,l(r)+1)
                                        if (abs(MDLp)<1.0e-30) cycle ! check if MD coeff is 0
                                        prod4 = MDLp * prod3
                                        ! MD coeff 5
                                        do mmp = 0, m(k)+m(r)
                                            MDMp = Dy2(mmp+1,m(k)+1,m(r)+1)
                                            if (abs(MDMp)<1.0e-30) cycle ! check if MD coeff is 0
                                            prod5 = MDMp * prod4
                                            ! MD coeff 6
                                            do nnp = 0, n(k)+n(r)
                                                MDNp = Dz2(nnp+1,n(k)+1,n(r)+1)
                                                prod6 = MDNp * prod5
                                                ! cutoff after MD
                                                if (abs(prod6)>cutoff2) then
                                                    ! add the contribution to the total
                                                    h_saved = h_saved + h_pre2(:,ll+llp+1,mm+mmp+1,nn+nnp+1)*prod6

                                                end if
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do

                    posr=posr+1
                    end do
                    posk=posk+1
                end do
                    posj=posj+1
            end do
        posi=posi+1
        end do

        CALL BesselSum(F, mu, H, LLmax, h_saved)

    end subroutine


        SUBROUTINE BesselSum(h_sum, mu, H, LLmax, h_saved)
                  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
    REAL(kind=dp), INTENT(in),  DIMENSION(:)            :: mu
    REAL(kind=dp), INTENT(in)                           :: H
    REAL(kind=dp), INTENT(in),  DIMENSION(LLmax+1)      :: h_saved
    INTEGER(kind=ikind), INTENT(in)                     :: LLmax
    REAL(kind=dp), INTENT(out), DIMENSION(size(mu))     :: h_sum


    ! helper variables
    REAL(kind=dp), DIMENSION(size(mu))     :: Pmu, h_0, h_1, h_r, muOH
    REAL(kind=dp)       :: coeff
    ! loop varialbes
    INTEGER(kind=ikind) :: ra

    Pmu=H * mu
    h_0=sin(Pmu)/Pmu
    coeff=h_saved(1)
    h_sum=h_0*coeff

    if (LLmax==1) then
        h_1=(sin(Pmu)/Pmu**2-cos(Pmu)/Pmu)*mu/H ! could be written in more efficient way by deleting mu
        coeff=h_saved(2)
        h_sum=h_sum+coeff*h_1
    elseif (LLmax>1) then
        muOH=mu/H
        h_1=(sin(Pmu)/Pmu**2-cos(Pmu)/Pmu)*muOH
        coeff=h_saved(2);
        h_sum=h_sum+coeff*h_1;
        do ra = 2, LLmax
            coeff=h_saved(ra+1);
            h_r= ((2*ra-1)/(Pmu)*h_1-h_0*muOH)*muOH;
            h_sum=h_sum+h_r*coeff;
            h_0=h_1;
            h_1=h_r;

        end do
    end if

    END SUBROUTINE BesselSum


        SUBROUTINE Hermite_like_coeffs(a, LLmax, Hx)
                  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        REAL(kind=dp), INTENT(out), DIMENSION(LLmax+1,LLmax+1) :: a
        INTEGER(kind=ikind), INTENT(in)     :: LLmax
        REAL(kind=dp), INTENT(in)           :: Hx

        ! loop vars
        INTEGER(kind=ikind)     :: L, ka
        a = 0
        a(1,1)=1;               ! L=0

        if (LLmax>0) a(2,2)=-Hx ! L=1

        do L = 2, LLmax
            a(L+1,1)=-(L-1)*a(L-1,1) ! L>1, ka = 0

            do ka=1, LLmax           ! L>1, ka > 0
                a(L+1,ka+1)=-Hx*a(L,ka)-(L-1)*a(L-1,ka+1);
            end do
        end do

    END SUBROUTINE Hermite_like_coeffs








    SUBROUTINE Bubble_Sort(a)

        implicit none
              INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
    INTEGER(kind=ikind), INTENT(inout), DIMENSION(:) :: a
    INTEGER(kind=ikind) :: temp
    INTEGER :: i, j
    LOGICAL :: swapped

     DO j = SIZE(a)-1, 1, -1
     swapped = .FALSE.
      DO i = 1, j
      IF (a(i) > a(i+1)) THEN
        temp = a(i)
        a(i) = a(i+1)
        a(i+1) = temp
        swapped = .TRUE.
      END IF
    END DO
    IF (.NOT. swapped) EXIT
  END DO
    END SUBROUTINE Bubble_Sort


          function kron(A,B) result(C)


       IMPLICIT NONE
             INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
       real(kind=dp), dimension (:), intent(in)  :: A, B
       real(kind=dp), dimension (:), allocatable :: C
       integer(kind=ikind) :: i = 0, j = 0, k = 0, l = 0, n = 0, m = 0



       allocate(C(size(A)*size(B)))
       C = 0.0_dp

       do i = 1,size(A)

         n=(i-1)*size(B) + 1
         m=n+size(B) - 1
         C(n:m) = A(i)*B
       enddo

      end function kron
end module integrals