module integrals
    use calculate_form_factors
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
        call omp_set_num_threads(9)
        print*,OMP_get_num_threads()
        !First big loop

        tsi=0.0_dp
        if (any(isnan(p0matrix))) print*,'ouch'

        !$OMP PARALLEL do private(posI,posK,posJ,posR,spi,spj,spk,spr,zcontrred,zcontrred2,za,zb,cmat), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,j,k,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz, posits,cutoffmd,group_count,group_start) REDUCTION(+:tsi), &
        !$OMP & schedule(dynamic)



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

subroutine tot_integration_j2(ncap,px,py,pz,l,m,n,p0matrix,dx,dy,dz,z1,z2,group_start,group_count,group, &
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
        integer(kind=ikind) :: nq,i,j,k,r,count,ng,ii,jj,start1,stop1,start2,stop2
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
        if (any(isnan(p0matrix))) print*,'ouch'
        print*,maxval(abs(p0matrix(:,1,1,1))), maxval(abs(p0matrix(:,2,1,1)))

        start1=1
        stop1=ncap/2
        start2=ncap/2+1
        stop2=ncap
        !Try with three chops
        !$OMP PARALLEL do private(posI,posK,posJ,posR,spi,spj,spk,spr,zcontrred,zcontrred2,za,zb,cmat), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,j,k,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz, posits,cutoffmd,group_count,group_start) REDUCTION(+:tsi), &
        !$OMP & schedule(dynamic)

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
                            if (maxval(abs(f))>1E3) then
                           !   print*, 'pzero maxval-> ', maxval(abs(f))
                            endif
                        else

                            call tot_integral_k_ijkr_j2(q,l,m,n,group_start, group_count, hx, hy, hz, h,dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, k, r, &
                                    zcontrred, zcontrred2,  cutoffz, cutoffmd,f)

                            if (maxval(abs(f))>1E3) then
                              print*, 'normal int maxval-> ', maxval(abs(f))
                            endif



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

                        call tot_integral_k_ijkr_j2(q, l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
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

                        call tot_integral_k_ijkr_j2(q, l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
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
!                    call tot_integral_ijkr_pzero(nq, l,m,n,group_start, group_count, p0matrix, dx1,dy1,dz1,dx2,&
!                                dy2,dz2,i, i, k, k, &
!                                    zcontrred, zcontrred2,  cutoffz, cutoffmd, f)


                else

                    call tot_integral_k_ijkr_j2(q,l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
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

    end subroutine tot_integration_j2


    subroutine tot_integration_aligned(ncap,px,py,pz,l,m,n,dx,dy,dz,z1,z2,group_start,group_count,group, &
            cutoffz,cutoffmd, cutoffcentre,q_al,e12,exponent1,exponent2,tsi)



        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(kind=ikind), INTENT(IN) :: ncap
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:) :: l,m,n,group,group_start,group_count

        REAL(kind=dp), intent(in), dimension(:,:,:,:,:),target :: dx,dy,dz
        complex(kind=dp), intent(in), dimension(:,:,:,:):: exponent1,exponent2
        real(kind=dp),dimension(:,:,:), pointer :: dx1,dy1,dz1,dx2,dy2,dz2

        REAL(kind=dp), intent(in), dimension(:,:,:) ::e12
        real(kind=dp), intent(in), dimension(:,:,:)::z1,z2
        real(kind=dp), dimension(:,:,:,:), allocatable::zcontrred,zcontrred2
        REAL(kind=dp), intent(in), dimension(:,:) :: px,py,pz
        REAL(kind=dp), intent(in), dimension(:,:) :: q_al
        REAL(kind=dp), intent(in) :: cutoffz, cutoffmd,cutoffcentre


        complex(kind=dp), dimension(size(q_al(1,:))) :: f
        integer(kind=ikind), dimension(:,:), allocatable :: posits
        real(kind=dp),dimension(:,:), allocatable :: za,zb,cmat
        integer(kind=ikind),dimension(:), allocatable ::posi,posj,posk,posr
        Complex(kind=dp), intent(out), dimension(size(q_al(1,:))) :: tsi
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

        nq= size(q_al(1,:))


        !First big loop

        tsi=0.0_dp


        do i=1,ncap
            do j=1,ncap
                do k=1,ncap
                    do r=1,ncap
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

                        call tot_integral_k_ijkr_alig(q_al, l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, k, r, &
                                zcontrred,  zcontrred2,  cutoffz, cutoffmd,f,exponent1,exponent2)





                        tsi = tsi + f * e12(:, i, j) * e12(:, k, r)


                        deallocate(zcontrred, zcontrred2)
               !        deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)


                    end do
                end do
            end do
        end do



        print*,'la rata ', tsi(1)


!        do i=1,ncap
!            do j=i+1,ncap
!                do r=i+1,ncap
!                    hx = px(i, r) - px(i, j)
!                    hy = py(i, r) - py(i, j)
!                    hz = pz(i, r) - pz(i, j)
!                    h = sqrt((hx * hx + hy * hy + hz * hz))
!                    allocate(posI(size(posits(i,:group_count(i)))),posJ(size(posits(j,:group_count(j)))) &
!                               ,posR(size(posits(r,:group_count(r)))))
!
!
!                        posI = posits(i,:group_count(i))
!                        posJ = posits(j,:group_count(j))
!
!                        posR = posits(r,:group_count(r))
!
!                        spi = size(posI)
!                        spj = size(posJ)
!                        spr = size(posR)
!
!                        allocate(zcontrred(spj, spi, spr, spi), za(szo, spj), zb(szo, spi), &
!                       &         cmat(spj, spi), zcontrred2(spr, spi, spj, spi))
!
!                          do ii = 1, spi
!                            do jj = 1, spr
!                              ! za = transpose(z1(:,posj,posi(ii)))
!                                za = z1(:,posj,posi(ii))
!                                zb = z2(:,posi,posr(jj))
!                               !cmat = matmul(za,zb)
!                                call dgemm('t','n', spj, spi, szo, 1.0_dp/8.0_dp, za, &
!                               &           szo, zb, szo, 0.0_dp, cmat, spj)
!                                zcontrred(:,:,jj,ii) = cmat
!                                  za = z2(:,posj,posi(ii))
!                              !  za = transpose(z1(:,posj,posi(ii)))
!                                zb = z1(:,posi,posr(jj))
!                               ! cmat = matmul(za,zb)
!                                call dgemm('t','n', spj, spi, szo, 1.0_dp/8.0_dp, za, &
!                              &           szo, zb, szo, 0.0_dp, cmat, spj)
!
!                                zcontrred2(jj,ii,:,:) = cmat
!                            enddo
!                          enddo
!
!!
!                    dx1=>dx(:,:,:,j,i)
!                    dy1=>dy(:,:,:,j,i)
!                    dz1=>dz(:,:,:,j,i)
!                    dx2=>dx(:,:,:,r,i)
!                    dy2=>dy(:,:,:,r,i)
!                    dz2=>dz(:,:,:,r,i)
!
!
!
!
!
!                        call tot_integral_k_ijkr_alig(q_al, l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
!                                dy2,dz2,i, j, i, r, &
!                                zcontrred,  zcontrred2,  cutoffz, cutoffmd,f,exponent1,exponent2)
!
!
!
!
!                    tsi = tsi + 4.000 * f * e12(:, i, j) * e12(:, i, r)
!                    count=count+1
!                    deallocate(posI,posJ,posR,za,zb,cmat)
!                    deallocate(zcontrred, zcontrred2)
!                  !   deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)
!
!                end do
!            end do
!
!        end do
!
!        print*,'here'
!
!        do i=1,ncap
!            do k=1,ncap
!                do r=k+1,ncap
!                    hx = px(k, r) - px(i, i)
!                    hy = py(k, r) - py(i, i)
!                    hz = pz(k, r) - pz(i, i)
!                    h = sqrt((hx * hx + hy * hy + hz * hz))
!
!                    allocate(posI(size(posits(i,:group_count(i)))) &
!                                ,posK(size(posits(k,:group_count(k)))),posR(size(posits(r,:group_count(r)))))
!
!
!                        posI = posits(i,:group_count(i))
!
!                        posK = posits(k,:group_count(k))
!                        posR = posits(r,:group_count(r))
!
!                    spi = size(posI)
!                    spk = size(posK)
!                    spr = size(posR)
!
!                    allocate(zcontrred(spi, spk, spr, spi), za(szo, spi), zb(szo, spk), &
!                   &         cmat(spi, spk), zcontrred2(spr, spi, spi, spk))
!
!                    do ii = 1, spi
!                        do jj = 1, spr
!                          !  za = transpose(z1(:,posi,posi(ii)))
!                            za = z1(:,posi,posi(ii))
!                            zb = z2(:,posk,posr(jj))
!                           ! cmat = matmul(za,zb)
!                            call  dgemm('t','n', spi, spk, szo, 1.0_dp/8.0_dp, za, &
!                           &           szo, zb, szo, 0.0_dp, cmat, spi)
!                            zcontrred(:,:,jj,ii) = cmat
!                             za = z2(:,posi,posi(ii))
!                              !  za = transpose(z1(:,posj,posi(ii)))
!                                zb = z1(:,posk,posr(jj))
!                               ! cmat = matmul(za,zb)
!                                call dgemm('t','n', spi, spk, szo, 1.0_dp/8.0_dp, za, &
!                              &           szo, zb, szo, 0.0_dp, cmat, spi)
!
!                                zcontrred2(jj,ii,:,:) = cmat
!                        enddo
!                    enddo
!
!
!                    dx1=>dx(:,:,:,i,i)
!                    dy1=>dy(:,:,:,i,i)
!                    dz1=>dz(:,:,:,i,i)
!                    dx2=>dx(:,:,:,r,k)
!                    dy2=>dy(:,:,:,r,k)
!                    dz2=>dz(:,:,:,r,k)
!!                    zcontrred=zcontrred
!!                    zcontrred2=zcontrred2
!
!
!
!                        call tot_integral_k_ijkr_alig(q_al, l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
!                                dy2,dz2,i, i, k, r, &
!                                    zcontrred, zcontrred2,  cutoffz, cutoffmd, f,exponent1,exponent2)
!
!
!
!
!                    tsi = tsi+ 4.000 * f * e12(:, i, i) * e12(:, k, r)
!                    count=count+1
!                    deallocate(posI,posK,posR,za,zb,cmat)
!                    deallocate(zcontrred, zcontrred2)
!               !     deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)
!
!                end do
!            end do
!        end do
!
!
!        print*,'here 2'
!
!        do i=1,ncap
!            do k=i+1,ncap
!
!                hx = px(k, k) - px(i, i)
!                hy = py(k, k) - py(i, i)
!                hz = pz(k, k) - pz(i, i)
!                h = sqrt((hx * hx + hy * hy + hz * hz))
!
!                allocate(posI(size(posits(i,:group_count(i)))), &
!                        posK(size(posits(k,:group_count(k)))))
!
!
!                posI = posits(i,:group_count(i))
!
!                posK = posits(k,:group_count(k))
!
!
!                spi = size(posI)
!                spk = size(posK)
!
!                allocate(zcontrred(spi, spk, spk, spi), za(szo, spi), zb(szo, spk), &
!               &         cmat(spi, spk), zcontrred2(spk, spi, spi, spk))
!
!                do ii = 1, spi
!                    do jj = 1, spk
!                       ! za = transpose(z1(:,posi,posi(ii)))
!                        za = z1(:,posi,posi(ii))
!                        zb = z2(:,posk,posk(jj))
!!                        cmat = matmul(za,zb)
!                        call dgemm('t','n', spi, spk, szo, 1.0_dp/8.0_dp, za, &
!                       &           szo, zb, szo, 0.0_dp, cmat, spi)
!                        zcontrred(:,:,jj,ii) = cmat
!                        za = z2(:,posi,posi(ii))
!                              !  za = transpose(z1(:,posj,posi(ii)))
!                        zb = z1(:,posk,posk(jj))
!                               ! cmat = matmul(za,zb)
!                        call dgemm('t','n', spi, spk, szo, 1.0_dp/8.0_dp, za, &
!                              &           szo, zb, szo, 0.0_dp, cmat, spi)
!
!                        zcontrred2(jj,ii,:,:) = cmat
!                    enddo
!                enddo
!
!
!                dx1=>dx(:,:,:,i,i)
!                dy1=>dy(:,:,:,i,i)
!                dz1=>dz(:,:,:,i,i)
!                dx2=>dx(:,:,:,k,k)
!                dy2=>dy(:,:,:,k,k)
!                dz2=>dz(:,:,:,k,k)
!
!
!
!
!                    call tot_integral_k_ijkr_alig(q_al,l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
!                                dy2,dz2,i, i, k, k, &
!                                    zcontrred, zcontrred2,  cutoffz, cutoffmd, f,exponent1,exponent2)
!
!                tsi = tsi+ 2.000 * f * e12(:, i, i) * e12(:, k, k)
!                deallocate(posI,posK,za,zb,cmat)
!                deallocate(zcontrred, zcontrred2)
!               ! deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)
!
!            end do
!        end do
!
!
!
!        print*,'here 3'
!
!        do i=1,ncap
!                              allocate(posI(size(posits(i,:group_count(i)))))
!
!
!                        posI = posits(i,:group_count(i))
!
!
!
!            spi = size(posI)
!
!            allocate(zcontrred(spi, spi, spi, spi), za(szo, spi), zb(szo, spi), &
!           &         cmat(spi, spi), zcontrred2(spi, spi, spi, spi))
!
!            do ii = 1, spi
!                do jj = 1, spi
!                  !  za = transpose(z1(:,posi,posi(ii)))
!                    za = z1(:,posi,posi(ii))
!                    zb = z2(:,posi,posi(jj))
!                   ! cmat = matmul(za,zb)
!                    call dgemm('t','n', spi, spi, szo, 1.0_dp/8.0_dp, za, &
!                   &           szo, zb, szo, 0.0_dp, cmat, spi)
!                    zcontrred(:,:,jj,ii) = cmat
!                     za = z2(:,posi,posi(ii))
!                              !  za = transpose(z1(:,posj,posi(ii)))
!                                zb = z1(:,posi,posi(jj))
!                               ! cmat = matmul(za,zb)
!                                call dgemm('t','n', spi, spi, szo, 1.0_dp/8.0_dp, za, &
!                              &           szo, zb, szo, 0.0_dp, cmat, spi)
!
!                                zcontrred2(jj,ii,:,:) = cmat
!                enddo
!            enddo
!
!            dx1=>dx(:,:,:,i,i)
!            dy1=>dy(:,:,:,i,i)
!            dz1=>dz(:,:,:,i,i)
!            dx2=>dx(:,:,:,i,i)
!            dy2=>dy(:,:,:,i,i)
!            dz2=>dz(:,:,:,i,i)
!!            zcontrred=zcontrred/8.0
!!            zcontrred2=zcontrred2/8.0
!
!               call tot_integral_k_ijkr_alig(q_al, l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
!                                dy2,dz2,i, i, i,i , &
!                                zcontrred,  zcontrred2,  cutoffz, cutoffmd,f,exponent1,exponent2)
!
!
!            tsi = tsi + f * e12(:, i, i) * e12(:, i, i)
!            count=count+1
!            deallocate(posI,za,zb,cmat)
!            deallocate(zcontrred, zcontrred2)
!
!
!        end do

    print*, tsi(1)

    end subroutine tot_integration_aligned

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
        real(kind=dp),dimension(:,:,:,:), intent(in) :: p0matrix
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
                       ! if(abs(f(1))>=1E-30)print*,f(1)
                        tsi = tsi + 8.000 * f * e12(:, i, j) * e12(:, k, r)

                        count=count+1
                        deallocate(za,zb)
               !        deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)


                    end do
                end do
            end do
        end do

      !$OMP END parallel DO


        print*,'la rata ', tsi

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

    subroutine elastic_integration_j2(ncap,px,py,pz,l,m,n,p0matrix,dx,dy,dz,z,group_start,group_count,group, &
            cutoffz,cutoffmd, cutoffcentre,q,e12,tsi)


        use omp_lib
        implicit none

              INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(kind=ikind), INTENT(IN) :: ncap
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:) :: l,m,n,group,group_start,group_count

        REAL(kind=dp), intent(in), dimension(:,:,:,:,:),target :: dx,dy,dz
        real(kind=dp),dimension(:,:,:), pointer :: dx1,dy1,dz1,dx2,dy2,dz2
        real(kind=dp),dimension(:,:,:,:), intent(in) :: p0matrix
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
     ! !$OMP PARALLEL do private(posI,posK,posJ,posR,spi,spj,spk,spr,za,zb), &
     !   !$OMP& private(f,ii,jj,h,hx,hy,hz,i,j,k,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
      !  !$OMP& shared( cutoffz, posits,cutoffmd,group_count,group_start,z) REDUCTION(+:tsi), &
      !  !$OMP & schedule(dynamic)

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

                            call integral_k_ijkr_j2(q,l,m,n,group_start, group_count, hx, hy, hz, h,dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, k, r, &
                                    za,zb,  cutoffz, cutoffmd,f)





                        end if
                       ! if(abs(f(1))>=1E-30)print*,f(1)
                        tsi = tsi + 8.000 * f * e12(:, i, j) * e12(:, k, r)

                        count=count+1
                        deallocate(za,zb)
               !        deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)


                    end do
                end do
            end do
        end do

      !!$OMP END parallel DO


        print*,'la rata ', tsi

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

                        call  integral_k_ijkr_j2(q, l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
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

                        call  integral_k_ijkr_j2(q, l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
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

                    call  integral_k_ijkr_j2(q,l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
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
    end subroutine elastic_integration_j2

    subroutine elastic_integration_alig(ncap,px,py,pz,l,m,n,dx,dy,dz,z,group_start,group_count,group, &
            cutoffz,cutoffmd, cutoffcentre,q_al,e12,exponent1,tsi)


        use omp_lib
        implicit none

              INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(kind=ikind), INTENT(IN) :: ncap
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:) :: l,m,n,group,group_start,group_count

        REAL(kind=dp), intent(in), dimension(:,:,:,:,:),target :: dx,dy,dz
        complex(kind=dp), intent(in), dimension(:,:,:,:) :: exponent1
        real(kind=dp),dimension(:,:,:), pointer :: dx1,dy1,dz1,dx2,dy2,dz2

        REAL(kind=dp), intent(in), dimension(:,:,:) ::e12
        real(kind=dp), intent(in), dimension(:,:)::z

        REAL(kind=dp), intent(in), dimension(:,:) :: px,py,pz
        REAL(kind=dp), intent(in), dimension(:,:) :: q_al
        REAL(kind=dp), intent(in) :: cutoffz, cutoffmd,cutoffcentre


        complex(kind=dp), dimension(size(q_al(1,:))) :: f
        integer(kind=ikind), dimension(:,:), allocatable :: posits
        real(kind=dp),dimension(:,:), allocatable :: za,zb
        integer(kind=ikind),dimension(:), allocatable ::posi,posj,posk,posr
        complex(kind=dp), intent(out), dimension(size(q_al(1,:))) :: tsi
        real(kind=dp) :: hx,hy,hz,h
        integer(kind=ikind),dimension(size(l)) :: ll
        integer(kind=ikind) :: nq,i,j,k,r,count,ng,ii,jj
        integer(kind=ikind) :: spi, spj, spk, spr,nt,counter1,counter2

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




        nq= size(q_al(1,:))

      !  print*,OMP_get_num_threads()
      !  call omp_set_num_threads(16)
      !  print*,OMP_get_num_threads()
        !First big loop

        tsi=(0.0_dp,0.0_dp)
     !  !$OMP PARALLEL do private(posI,posJ,spi,spj,za), &
      !  !$OMP& private(f,ii,jj,h,hx,hy,hz,i,j,dx1,dy1,dz1) shared(q_al,l,m,n), &
     !   !$OMP& shared( cutoffz, posits,cutoffmd,group_count,group_start,z) REDUCTION(+:tsi)
        print*,Z(2,9), E12(1,16,16)
        counter1=0

        do i=1,ncap
             counter1=counter1+group_count(i)
             counter2=0
            do j=1,ncap

                        counter2=counter2+group_count(j)
                        hx =  px(i, j)
                        hy =  py(i, j)
                        hz =  pz(i, j)
                        h = (hx * hx + hy * hy + hz * hz)**0.5

                        allocate(posI(size(posits(i,:group_count(i)))),posJ(size(posits(j,:group_count(j)))))


                        posI = posits(i,:group_count(i))
                        posJ = posits(j,:group_count(j))


                        spi = size(posI)
                        spj = size(posJ)


                        allocate(za(spi, spj))

                         za=z(posI,posJ)


                         deallocate(posI,posJ)
                         dx1=>dx(:,:,:,j,i)
                         dy1=>dy(:,:,:,j,i)
                         dz1=>dz(:,:,:,j,i)



                        call integral_k_ijkr_alig(q_al,exponent1,l,m,n,group_start, group_count, hx, hy, hz, h,dx1,dy1,dz1,i, j &
                                    ,za, cutoffz, cutoffmd,f)






                        tsi = tsi + f * e12(:, i, j)

                        count=count+1
                        deallocate(za)
               !        deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)


                    end do
                end do
        print*,tsi(1),q_al(1,1), q_al(2,1), q_al(3,1)
      ! !$OMP END parallel DO



    end subroutine elastic_integration_alig


    subroutine nuclei_electron_integration(Zn,geom,ncap,px,py,pz,l,m,n,p0matrix,dx,dy,dz,z,group_start,group_count,group, &
            cutoffz,cutoffmd, cutoffcentre,q,e12,tsi)
             use omp_lib
        implicit none

              INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(kind=ikind), INTENT(IN) :: ncap
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:) :: l,m,n,group,group_start,group_count

        REAL(kind=dp), intent(in), dimension(:,:,:,:,:),target :: dx,dy,dz
        real(kind=dp),dimension(:,:,:), pointer :: dx1,dy1,dz1
        real(kind=dp),dimension(:,:,:,:), intent(in) :: p0matrix
        REAL(kind=dp), intent(in), dimension(:,:,:) ::e12
        real(kind=dp), intent(in), dimension(:,:)::z,geom

        REAL(kind=dp), intent(in), dimension(:,:) :: px,py,pz
        REAL(kind=dp), intent(in), dimension(:) :: q,Zn
        REAL(kind=dp), intent(in) :: cutoffz, cutoffmd,cutoffcentre


        REAL(kind=dp), dimension(size(q)) :: f
        integer(kind=ikind), dimension(:,:), allocatable :: posits
        real(kind=dp),dimension(:,:), allocatable :: za,zb
        integer(kind=ikind),dimension(:), allocatable ::posi,posk
        REAL(kind=dp), intent(out), dimension(size(q)) :: tsi
        real(kind=dp) :: hx,hy,hz,h
        integer(kind=ikind),dimension(size(l)) :: ll
        integer(kind=ikind) :: nq,i,j,k,r,count,ng,ii,jj
        integer(kind=ikind) :: spi, spk, nt,na

        nq=size(q)
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


             tsi=0.0_dp
        do na=1,size(Zn)
            do i=1,ncap
                 do k=1,ncap

                    hx = geom(na,1) - px(i, k)
                    hy = geom(na,2) - py(i, k)
                    hz = geom(na,3) - pz(i, k)
                    h = sqrt((hx * hx + hy * hy + hz * hz))

                allocate(posI(size(posits(i,:group_count(i)))), &
                        posK(size(posits(k,:group_count(k)))))


                posI = posits(i,:group_count(i))

                posK = posits(k,:group_count(k))


                spi = size(posI)
                spk = size(posK)

                allocate(za(spi,spk))

               za=z(posi,posk)



                dx1=>dx(:,:,:,k,i)
                dy1=>dy(:,:,:,k,i)
                dz1=>dz(:,:,:,k,i)


                if (h < cutoffcentre) then
                    call integral_ij_pzero(nq, l,m,n,group_start, group_count, p0matrix, dx1,dy1,dz1,i, k, &
                                    za,  cutoffz, cutoffmd, f)

                else

                    call  integral_k_ij(q,l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,&
                                 i, k, za,  cutoffz, cutoffmd, f)
                end if

                tsi = tsi +  Zn(na) * f * e12(:, i, k)

                deallocate(posI,posK,za)

               ! deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

            end do
        end do

     enddo
     print*,tsi(90)


    end subroutine nuclei_electron_integration



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


     SUBROUTINE set_P0_j2(P0, lmax4, q)

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

        p0=0.0_dp
        ! loop for all combinations of angular momenta
        do k = 0, lmax4
            do j = 0, lmax4
                do i = 0, lmax4

                    if ((i+j+j)<16) then
                    ! Fill the current LMN
                    CALL P_LMN_j2(P0, k, j, i, q)

                    end if

                end do
            end do
        end do

         P0=P0*sqrt(5.0_dp/dacos(-1.0_dp))
    END SUBROUTINE set_P0_j2

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
            P0(:,L+1,M+1,N+1)=-(q**2/3.)
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



    SUBROUTINE P_LMN_j2(P0, L, M, N, q)

        ! the three values of the angular momentumS
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(KIND=ikind), INTENT(IN) :: L, M, N
        ! The value of <qx^L*qy^M*qz^N>
        REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:) :: P0
        ! The radial part of the scattering vector
        REAL(KIND=dp), INTENT(IN), DIMENSION(:)     :: q
        real(kind=dp):: pi

        ! For ordering
        INTEGER(KIND=ikind), DIMENSION(1:3)     :: A
        ! helper variable
        !INTEGER(KIND=ikind) :: i
        pi=dacos(-1.0_dp)
        ! Order L, M and N
        A(1) = L
        A(2) = M
        A(3) = N
        CALL Bubble_Sort(A)

        ! These are analytical solutions to the integralif (mod(L,2)/=0 .or. mod(M,2)/=0 .or. mod(N,2)/=0) then
        ! They have been obtained in Matematica by Andres Moreno
        ! They have been extensively debugged in the MATLAB version of the code
        ! however one should be careful with the precision as there is nasty
        ! divisions and exponentiations. Needs to be retested if it works well
        ! with the current data kinds.
     if (mod(L,2)/=0 .or. mod(M,2)/=0 .or. mod(N,2)/=0) then
    P0(:,L+1,M+1,N+1)=0.0_dp
elseif ((L+M)==2*N) then
    P0(:,L+1,M+1,N+1)=0.0_dp
elseif (L==2 .and. M==0 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/30.0_dp)*q**2.0_dp
elseif (L==0 .and. M==2 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/30.0_dp)*q**2.0_dp
elseif (L==0 .and. M==0 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/15.0_dp)*q**2.0_dp
elseif (L==2 .and. M==2 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/105.0_dp)*q**4.0_dp
elseif (L==2 .and. M==0 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/210.0_dp)*q**4.0_dp
elseif (L==0 .and. M==2 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/210.0_dp)*q**4.0_dp
elseif (L==4 .and. M==0 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/35.0_dp)*q**4.0_dp
elseif (L==0 .and. M==4 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/35.0_dp)*q**4.0_dp
elseif (L==0 .and. M==0 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(-2.0_dp/35.0_dp)*q**4.0_dp
elseif (L==4 .and. M==2 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/210.0_dp)*q**6.0_dp
elseif (L==2 .and. M==4 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/210.0_dp)*q**6.0_dp
elseif (L==2 .and. M==0 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/210.0_dp)*q**6.0_dp
elseif (L==0 .and. M==2 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/210.0_dp)*q**6.0_dp
elseif (L==6 .and. M==0 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/42.0_dp)*q**6.0_dp
elseif (L==0 .and. M==6 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/42.0_dp)*q**6.0_dp
elseif (L==0 .and. M==0 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/21.0_dp)*q**6.0_dp
elseif (L==4 .and. M==2 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/6930.0_dp)*q**8.0_dp
elseif (L==2 .and. M==4 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/6930.0_dp)*q**8.0_dp
elseif (L==2 .and. M==2 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/3465.0_dp)*q**8.0_dp
elseif (L==4 .and. M==4 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(2.0_dp/1155.0_dp)*q**8.0_dp
elseif (L==4 .and. M==0 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/1155.0_dp)*q**8.0_dp
elseif (L==0 .and. M==4 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/1155.0_dp)*q**8.0_dp
elseif (L==6 .and. M==2 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(2.0_dp/693.0_dp)*q**8.0_dp
elseif (L==6 .and. M==0 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/1386.0_dp)*q**8.0_dp
elseif (L==2 .and. M==6 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(2.0_dp/693.0_dp)*q**8.0_dp
elseif (L==0 .and. M==6 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/1386.0_dp)*q**8.0_dp
elseif (L==2 .and. M==0 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(-5.0_dp/1386.0_dp)*q**8.0_dp
elseif (L==0 .and. M==2 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(-5.0_dp/1386.0_dp)*q**8.0_dp
elseif (L==8 .and. M==0 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(2.0_dp/99.0_dp)*q**8.0_dp
elseif (L==0 .and. M==8 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(2.0_dp/99.0_dp)*q**8.0_dp
elseif (L==0 .and. M==0 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(-4.0_dp/99.0_dp)*q**8.0_dp
elseif (L==4 .and. M==4 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/15015.0_dp)*q**10.0_dp
elseif (L==4 .and. M==2 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/30030.0_dp)*q**10.0_dp
elseif (L==2 .and. M==4 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/30030.0_dp)*q**10.0_dp
elseif (L==6 .and. M==2 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/9009.0_dp)*q**10.0_dp
elseif (L==2 .and. M==6 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/9009.0_dp)*q**10.0_dp
elseif (L==2 .and. M==2 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(2.0_dp/9009.0_dp)*q**10.0_dp
elseif (L==6 .and. M==4 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-5.0_dp/6006.0_dp)*q**10.0_dp
elseif (L==4 .and. M==6 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-5.0_dp/6006.0_dp)*q**10.0_dp
elseif (L==6 .and. M==0 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/6006.0_dp)*q**10.0_dp
elseif (L==0 .and. M==6 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/6006.0_dp)*q**10.0_dp
elseif (L==4 .and. M==0 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(2.0_dp/3003.0_dp)*q**10.0_dp
elseif (L==0 .and. M==4 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(2.0_dp/3003.0_dp)*q**10.0_dp
elseif (L==8 .and. M==2 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-5.0_dp/2574.0_dp)*q**10.0_dp
elseif (L==8 .and. M==0 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/1287.0_dp)*q**10.0_dp
elseif (L==2 .and. M==8 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-5.0_dp/2574.0_dp)*q**10.0_dp
elseif (L==0 .and. M==8 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/1287.0_dp)*q**10.0_dp
elseif (L==2 .and. M==0 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(7.0_dp/2574.0_dp)*q**10.0_dp
elseif (L==0 .and. M==2 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(7.0_dp/2574.0_dp)*q**10.0_dp
elseif (L==10 .and. M==0 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-5.0_dp/286.0_dp)*q**10.0_dp
elseif (L==0 .and. M==10 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-5.0_dp/286.0_dp)*q**10.0_dp
elseif (L==0 .and. M==0 .and. N==10) then
    P0(:,L+1,M+1,N+1)=(5.0_dp/143.0_dp)*q**10.0_dp
elseif (L==6 .and. M==4 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/30030.0_dp)*q**12.0_dp
elseif (L==4 .and. M==6 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/30030.0_dp)*q**12.0_dp
elseif (L==4 .and. M==2 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/30030.0_dp)*q**12.0_dp
elseif (L==2 .and. M==4 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/30030.0_dp)*q**12.0_dp
elseif (L==8 .and. M==2 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/12870.0_dp)*q**12.0_dp
elseif (L==6 .and. M==6 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/3003.0_dp)*q**12.0_dp
elseif (L==2 .and. M==8 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/12870.0_dp)*q**12.0_dp
elseif (L==6 .and. M==0 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/6006.0_dp)*q**12.0_dp
elseif (L==0 .and. M==6 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/6006.0_dp)*q**12.0_dp
elseif (L==2 .and. M==2 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/6435.0_dp)*q**12.0_dp
elseif (L==8 .and. M==4 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/2145.0_dp)*q**12.0_dp
elseif (L==4 .and. M==8 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/2145.0_dp)*q**12.0_dp
elseif (L==4 .and. M==0 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/2145.0_dp)*q**12.0_dp
elseif (L==0 .and. M==4 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/2145.0_dp)*q**12.0_dp
elseif (L==10 .and. M==2 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/715.0_dp)*q**12.0_dp
elseif (L==10 .and. M==0 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/1430.0_dp)*q**12.0_dp
elseif (L==2 .and. M==10 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/715.0_dp)*q**12.0_dp
elseif (L==0 .and. M==10 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/1430.0_dp)*q**12.0_dp
elseif (L==2 .and. M==0 .and. N==10) then
    P0(:,L+1,M+1,N+1)=(-3.0_dp/1430.0_dp)*q**12.0_dp
elseif (L==0 .and. M==2 .and. N==10) then
    P0(:,L+1,M+1,N+1)=(-3.0_dp/1430.0_dp)*q**12.0_dp
elseif (L==12 .and. M==0 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/65.0_dp)*q**12.0_dp
elseif (L==0 .and. M==12 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/65.0_dp)*q**12.0_dp
elseif (L==0 .and. M==0 .and. N==12) then
    P0(:,L+1,M+1,N+1)=(-2.0_dp/65.0_dp)*q**12.0_dp
elseif (L==6 .and. M==4 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/510510.0_dp)*q**14.0_dp
elseif (L==4 .and. M==6 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/510510.0_dp)*q**14.0_dp
elseif (L==4 .and. M==4 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/255255.0_dp)*q**14.0_dp
elseif (L==6 .and. M==6 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-2.0_dp/153153.0_dp)*q**14.0_dp
elseif (L==6 .and. M==2 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/153153.0_dp)*q**14.0_dp
elseif (L==2 .and. M==6 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/153153.0_dp)*q**14.0_dp
elseif (L==8 .and. M==4 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-2.0_dp/109395.0_dp)*q**14.0_dp
elseif (L==8 .and. M==2 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/218790.0_dp)*q**14.0_dp
elseif (L==4 .and. M==8 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-2.0_dp/109395.0_dp)*q**14.0_dp
elseif (L==2 .and. M==8 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/218790.0_dp)*q**14.0_dp
elseif (L==4 .and. M==2 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/43758.0_dp)*q**14.0_dp
elseif (L==2 .and. M==4 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/43758.0_dp)*q**14.0_dp
elseif (L==8 .and. M==6 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/43758.0_dp)*q**14.0_dp
elseif (L==6 .and. M==8 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/43758.0_dp)*q**14.0_dp
elseif (L==8 .and. M==0 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/21879.0_dp)*q**14.0_dp
elseif (L==0 .and. M==8 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/21879.0_dp)*q**14.0_dp
elseif (L==6 .and. M==0 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(5.0_dp/43758.0_dp)*q**14.0_dp
elseif (L==0 .and. M==6 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(5.0_dp/43758.0_dp)*q**14.0_dp
elseif (L==10 .and. M==2 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-2.0_dp/36465.0_dp)*q**14.0_dp
elseif (L==2 .and. M==10 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-2.0_dp/36465.0_dp)*q**14.0_dp
elseif (L==2 .and. M==2 .and. N==10) then
    P0(:,L+1,M+1,N+1)=(4.0_dp/36465.0_dp)*q**14.0_dp
elseif (L==10 .and. M==4 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/24310.0_dp)*q**14.0_dp
elseif (L==10 .and. M==0 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/24310.0_dp)*q**14.0_dp
elseif (L==4 .and. M==10 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/24310.0_dp)*q**14.0_dp
elseif (L==0 .and. M==10 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/24310.0_dp)*q**14.0_dp
elseif (L==4 .and. M==0 .and. N==10) then
    P0(:,L+1,M+1,N+1)=(4.0_dp/12155.0_dp)*q**14.0_dp
elseif (L==0 .and. M==4 .and. N==10) then
    P0(:,L+1,M+1,N+1)=(4.0_dp/12155.0_dp)*q**14.0_dp
elseif (L==12 .and. M==2 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/6630.0_dp)*q**14.0_dp
elseif (L==12 .and. M==0 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-2.0_dp/3315.0_dp)*q**14.0_dp
elseif (L==2 .and. M==12 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/6630.0_dp)*q**14.0_dp
elseif (L==0 .and. M==12 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(-2.0_dp/3315.0_dp)*q**14.0_dp
elseif (L==2 .and. M==0 .and. N==12) then
    P0(:,L+1,M+1,N+1)=(11.0_dp/6630.0_dp)*q**14.0_dp
elseif (L==0 .and. M==2 .and. N==12) then
    P0(:,L+1,M+1,N+1)=(11.0_dp/6630.0_dp)*q**14.0_dp
elseif (L==14 .and. M==0 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/510.0_dp)*q**14.0_dp
elseif (L==0 .and. M==14 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/510.0_dp)*q**14.0_dp
elseif (L==0 .and. M==0 .and. N==14) then
    P0(:,L+1,M+1,N+1)=(7.0_dp/255.0_dp)*q**14.0_dp
elseif (L==6 .and. M==6 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/969969.0_dp)*q**16.0_dp
elseif (L==6 .and. M==4 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/1939938.0_dp)*q**16.0_dp
elseif (L==4 .and. M==6 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/1939938.0_dp)*q**16.0_dp
elseif (L==8 .and. M==4 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/692835.0_dp)*q**16.0_dp
elseif (L==4 .and. M==8 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/692835.0_dp)*q**16.0_dp
elseif (L==4 .and. M==4 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(-2.0_dp/692835.0_dp)*q**16.0_dp
elseif (L==8 .and. M==6 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(5.0_dp/831402.0_dp)*q**16.0_dp
elseif (L==6 .and. M==8 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(5.0_dp/831402.0_dp)*q**16.0_dp
elseif (L==8 .and. M==2 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/831402.0_dp)*q**16.0_dp
elseif (L==2 .and. M==8 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/831402.0_dp)*q**16.0_dp
elseif (L==6 .and. M==2 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(-2.0_dp/415701.0_dp)*q**16.0_dp
elseif (L==2 .and. M==6 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(-2.0_dp/415701.0_dp)*q**16.0_dp
elseif (L==10 .and. M==4 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/92378.0_dp)*q**16.0_dp
elseif (L==10 .and. M==2 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/230945.0_dp)*q**16.0_dp
elseif (L==4 .and. M==10 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/92378.0_dp)*q**16.0_dp
elseif (L==2 .and. M==10 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/230945.0_dp)*q**16.0_dp
elseif (L==4 .and. M==2 .and. N==10) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/461890.0_dp)*q**16.0_dp
elseif (L==2 .and. M==4 .and. N==10) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/461890.0_dp)*q**16.0_dp
elseif (L==8 .and. M==8 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(28.0_dp/415701.0_dp)*q**16.0_dp
elseif (L==8 .and. M==0 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(-14.0_dp/415701.0_dp)*q**16.0_dp
elseif (L==0 .and. M==8 .and. N==8) then
    P0(:,L+1,M+1,N+1)=(-14.0_dp/415701.0_dp)*q**16.0_dp
elseif (L==10 .and. M==6 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(4.0_dp/46189.0_dp)*q**16.0_dp
elseif (L==6 .and. M==10 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(4.0_dp/46189.0_dp)*q**16.0_dp
elseif (L==10 .and. M==0 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/92378.0_dp)*q**16.0_dp
elseif (L==0 .and. M==10 .and. N==6) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/92378.0_dp)*q**16.0_dp
elseif (L==6 .and. M==0 .and. N==10) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/92378.0_dp)*q**16.0_dp
elseif (L==0 .and. M==6 .and. N==10) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/92378.0_dp)*q**16.0_dp
elseif (L==12 .and. M==2 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/25194.0_dp)*q**16.0_dp
elseif (L==2 .and. M==12 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/25194.0_dp)*q**16.0_dp
elseif (L==2 .and. M==2 .and. N==12) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/12597.0_dp)*q**16.0_dp
elseif (L==12 .and. M==4 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(4.0_dp/20995.0_dp)*q**16.0_dp
elseif (L==12 .and. M==0 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/20995.0_dp)*q**16.0_dp
elseif (L==4 .and. M==12 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(4.0_dp/20995.0_dp)*q**16.0_dp
elseif (L==0 .and. M==12 .and. N==4) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/20995.0_dp)*q**16.0_dp
elseif (L==4 .and. M==0 .and. N==12) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/4199.0_dp)*q**16.0_dp
elseif (L==0 .and. M==4 .and. N==12) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/4199.0_dp)*q**16.0_dp
elseif (L==14 .and. M==2 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(4.0_dp/4845.0_dp)*q**16.0_dp
elseif (L==14 .and. M==0 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/1938.0_dp)*q**16.0_dp
elseif (L==2 .and. M==14 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(4.0_dp/4845.0_dp)*q**16.0_dp
elseif (L==0 .and. M==14 .and. N==2) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/1938.0_dp)*q**16.0_dp
elseif (L==2 .and. M==0 .and. N==14) then
    P0(:,L+1,M+1,N+1)=(-13.0_dp/9690.0_dp)*q**16.0_dp
elseif (L==0 .and. M==2 .and. N==14) then
    P0(:,L+1,M+1,N+1)=(-13.0_dp/9690.0_dp)*q**16.0_dp
elseif (L==16 .and. M==0 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(4.0_dp/323.0_dp)*q**16.0_dp
elseif (L==0 .and. M==16 .and. N==0) then
    P0(:,L+1,M+1,N+1)=(4.0_dp/323.0_dp)*q**16.0_dp
elseif (L==0 .and. M==0 .and. N==16) then
    P0(:,L+1,M+1,N+1)=(-8.0_dp/323.0_dp)*q**16.0_dp
    else

            print*, "CASE NOT PROGRAMED YET", L, M, N
            P0(:,L+1,M+1,N+1) = 0.0_dp
            !STOP

        end if
        !P0=P0*sqrt(5.0_dp/pi)
    END SUBROUTINE P_LMN_j2

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
                    ! Barrilero style:
                    ! hOrder = ceiling((LL+MM+NN+ii+jj+kk)/2_ikind)
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

    SUBROUTINE integral_ij_pzero(nq,l,m,n,gs,gc,p0mat,dx1,dy1,dz1,gi,gj,za, &
             cutoff1,cutoff2,f)
        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

        ! definition of input

        INTEGER(kind=ikind), INTENT(IN)                       :: nq, gi, gj
        REAL(kind=dp), INTENT(IN)                             :: cutoff1, cutoff2
        REAL(kind=dp), INTENT(IN), DIMENSION(:,:,:,:)               :: p0mat
        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1
         real(kind=dp), intent(in),  dimension(:,:) :: za

        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:)         :: gs,gc,l,m,n
        ! definition of output
        REAL(kind=dp), INTENT(OUT), DIMENSION(nq)             :: f
        ! definition of loop indices
        INTEGER(kind=ikind)                                   :: h1
        INTEGER(kind=ikind)                                   :: ll, mm, nn, llp, mmp, nnp
        ! definition of internal variables
        INTEGER(kind=ikind)                                   :: i,j,k,r, posi, posj
        REAL(kind=dp)                                         :: z1,ztot

        REAL(kind=dp)                                         :: mdl, mdm, mdn, mdlp, mdmp,mdnp
        real(kind=dp) :: prod6, prod5,prod4,prod3,prod2,prod1
!        REAL(kind=dp), DIMENSION(size(Z1(:,1,1)))                          :: zij1, zkr1
!        REAL(kind=dp), DIMENSION(size(Z2(:,1,1)))                          :: zij2, zkr2


        posI=1

       ! posI=apos(i)

        f=0.0_dp
        ! loop through all possible ways to get total angular momentum lmax1
        do i = gs(gi), gs(gi) + gc(gi) - 1

            posj=1

            do j = gs(gj), gs(gj) + gc(gj) - 1
                Ztot=Za(posI, posJ)

                ! continue only if larger
!                if (abs(ztot) < cutoff1) then
!                    posj=posj+1
!                    cycle
!                end if


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
                                    if (abs(prod3)>cutoff2) then
                                        ! add the contribution to the total
                                        F = F + prod3 * P0mat(:,ll+1,mm+1,nn+1)

                                    end if

                                end do

                            end do
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
                                                      !  print*,prod6, maxval(abs(P0mat(:,ll+llp+1,mm+mmp+1,nn+nnp+1))), ll+llp+1,mm+mmp+1,nn+nnp+1
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
    if (maxval(abs(F))>1E3) then
        stop
        end if
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

        a=0.0_dp
        b=0.0_dp
        c=0.0_dp

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


    subroutine tot_integral_k_ijkr_j2(mu,l,m,n,group_start,group_count,hx,hy,hz,h,dx1,dy1,dz1,dx2,dy2,dz2, gi,gj,gk, gr,&
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
        INTEGER(kind=ikind), parameter :: dim = 20
        real(kind=dp) :: prod1,prod2,prod3,prod4,prod5,prod6
        REAL(kind=dp), DIMENSION(dim,dim)           :: a, b, c,ap,bp,cp
        REAL(kind=dp), DIMENSION(dim)               :: h_saved
        REAL(kind=dp), DIMENSION(dim, dim, dim, dim)  :: h_pre2,h_ord2
        REAL(kind=dp), DIMENSION(dim)               :: BD
        integer(kind=ikind),dimension(dim)::horder



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

        a=0.0_dp
        b=0.0_dp
        c=0.0_dp
        ap=0.0_dp
        bp=0.0_dp
        cp=0.0_dp
      !  call rrdj2a(3,0.25_dp,a)
     !   call rrdj2b(a,ap)


       ! do k=1,10
      !  print*,ap(k,1:10)

     !   end do

     !   print*,LLmax
        call rrdj2a(LLmax,Hx,a)
        call rrdj2a(LLmax,Hy,b)
        call rrdj2a(LLmax,Hz,c)
        call rrdj2b(llmax,a,ap)
        call rrdj2b(llmax,b,bp)
        call rrdj2b(llmax,c,cp)








        bd=0.0_dp
        h_pre2=0.0_dp

        do k = 0, LLmax
            do j = 0, LLmax - k
                do i = 0, LLmax - k - j
                    call BesselDeriv2j(BD,Hz,H, i, j, k, a, b, c, ap,bp,cp,LLmax)
                    h_pre2(:,i+1,j+1,k+1) = BD


                end do
            end do
        end do

       ! posi=apos(i)


       ! print*,'max BD', maxval(h_pre2)
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
         !print*,h_saved(2)
    !    CALL BesselSum(F, mu, H, LLmax, h_saved)
        call bessels2(F,LLMAX,mu,H, h_saved)

    end subroutine

     subroutine tot_integral_k_ijkr_alig(q_al,l,m,n,group_start,group_count,hx,hy,hz,h,dx1,dy1,dz1,dx2,dy2,dz2, gi,gj,gk, gr,&
            zcontr, zcontr2, &
            cutoff1, cutoff2,f,exp1,exp2)



        implicit none


        !INTEGER, PARAMETER :: dp = selected_real_kind(2*precision(1.0_dp))
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=selected_int_kind(8)), intent(in)  :: gi,gj,gk,gr
        integer(kind=selected_int_kind(8)), dimension(:), intent(in) :: l,m,n,group_start,group_count
        real(kind=dp), intent(in)              :: cutoff1,cutoff2, hx, hy, hz, h
        real(kind=dp), intent(in), dimension(:,:,:,:) :: zcontr,zcontr2
        complex(kind=dp),intent(in),dimension(:,:,:,:) :: exp1,exp2
        real(kind=dp), intent(in), dimension(:,:)   ::  q_al

        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1,dx2,dy2,dz2
        !real(kind=dp), pointer,dimension(:,:,:) :: dxx,dyy,dzz, dxx2,dyy2,dzz2



        integer(kind=selected_int_kind(8)) :: ka, posi, posj, posk, posr, ra,i,j,k,r
        integer(kind=selected_int_kind(8)) ::  h1
        integer(kind=selected_int_kind(8)) ::  ll, mm, nn, llp, mmp, nnp,llmax

        real(kind=dp)  ::coeff,prodd,ztot,mdn, mdl, mdm, mdlp, mdmp,mdnp,z11,z22
        INTEGER(kind=ikind), parameter :: dim = 13
        real(kind=dp) :: prod1,prod2,prod3,prod4,prod5,prod6
        REAL(kind=dp), DIMENSION(dim,dim)           :: a, b, c
        complex(kind=dp), DIMENSION(size(q_al(1,:)))               :: h_saved
        REAL(kind=dp), DIMENSION(dim, dim, dim, dim)  :: h_pre2
        REAL(kind=dp), DIMENSION(dim)               :: BD



      !  real(kind=dp), dimension(:), allocatable :: pmu, h_0, h_1, h_r, muoh,zij, zij2, zkr, zkr2

        complex(kind=dp), intent(out), dimension(size(q_al(1,:))) :: f
        real(kind=dp), allocatable,dimension(:) :: l1vec,l2vec,suml1l2


        LLmax = l(group_start(gi)) + m(group_start(gi)) + n(group_start(gi)) + &
                l(group_start(gj)) + m(group_start(gj)) + n(group_start(gj)) + &
                l(group_start(gk)) + m(group_start(gk)) + n(group_start(gk)) + &
                l(group_start(gr)) + m(group_start(gr)) + n(group_start(gr))


        if (LLmax + 1 > dim) then
            print*, "only s,p,d and f type GTOs are supported"
            stop
        end if




       ! posi=apos(i)



        h_saved=(0.0_dp,0.0_dp)

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
                                    H1=1
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
                                                    h_saved = h_saved + exp1(:,llp+1,mmp+1,nnp+1)*exp2(:,ll+1,mm+1,nn+1)*prod6

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

       f=h_saved*exp((0,1.00_dp)*q_al(1,:)*(Hx))*exp((0,1.0_dp)*q_al(2,:)*(Hy))*exp((0,1.0_dp)*q_al(3,:)*(Hz))

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
        real(kind=dp), allocatable, dimension(:,:) :: bess
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


        call bessels(F,LLmax,mu,H,h_saved)
!        do i=1,LLmax+1
!            F=F+h_saved(i)*bess(i,:)
!
!        end do
      !   CALL BesselSum(F, mu, H, LLmax, h_saved)

    end subroutine

     subroutine integral_k_ijkr_j2(mu,l,m,n,group_start,group_count,hx,hy,hz,h,dx1,dy1,dz1,dx2,dy2,dz2, gi,gj,gk, gr,&
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



        integer(kind=selected_int_kind(8)) :: ka, posi, posj, posk, posr, ra,i,j,k,r,ind
        integer(kind=selected_int_kind(8)) ::  h1
        integer(kind=selected_int_kind(8)) ::  ll, mm, nn, llp, mmp, nnp,llmax

        real(kind=dp)  ::coeff,prodd,ztot,mdn, mdl, mdm, mdlp, mdmp,mdnp,z1
        INTEGER(kind=ikind), parameter :: dim = 20
        real(kind=dp) :: prod1,prod2,prod3,prod4,prod5,prod6
       REAL(kind=dp), DIMENSION(dim,dim)           :: a, b, c,ap,bp,cp
        REAL(kind=dp), DIMENSION(dim)               :: h_saved
        REAL(kind=dp), DIMENSION(dim, dim, dim, dim)  :: h_pre2
        REAL(kind=dp), DIMENSION(dim)               :: BD



      !  real(kind=dp), dimension(:), allocatable :: pmu, h_0, h_1, h_r, muoh,zij, zij2, zkr, zkr2

        real(kind=dp), intent(out), dimension(size(mu)) :: f
        real(kind=dp), allocatable, dimension(:,:) :: bess
        real(kind=dp), allocatable,dimension(:) :: l1vec,l2vec,suml1l2


        LLmax = l(group_start(gi)) + m(group_start(gi)) + n(group_start(gi)) + &
                l(group_start(gj)) + m(group_start(gj)) + n(group_start(gj)) + &
                l(group_start(gk)) + m(group_start(gk)) + n(group_start(gk)) + &
                l(group_start(gr)) + m(group_start(gr)) + n(group_start(gr))


        if (LLmax + 1 > dim) then
            print*, "only s,p,d and f type GTOs are supported"
            stop
        end if




        call rrdj2a(llmax,Hx,a)
        call rrdj2a(llmax,Hy,b)
        call rrdj2a(llmax,Hz,c)
        call rrdj2b(llmax,a,ap)
        call rrdj2b(llmax,b,bp)
        call rrdj2b(llmax,c,cp)
        cp=-2.0_dp*cp
!        print*,'AP2',a(3,5), ap(3,3)


        bd=0.0_dp

          do k = 0, LLmax
            do j = 0, LLmax-k
                do i = 0, LLmax-k-j
                    call BesselDeriv2j(BD,Hz,H, i, j, k, a, b, c, ap,bp,cp,LLmax)
                    h_pre2(:,i+1,j+1,k+1) = BD
!
                end do
            end do
        end do


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


       ! call bessels(F,LLmax,mu,H,h_saved)
!        do i=1,LLmax+1
!            F=F+h_saved(i)*bess(i,:)
!
!        end do
      !   CALL BesselSum(F, mu, H, LLmax, h_saved)
       if (sum(h_saved)==0.0_dp) then
           F=0.0_dp
           return
        end if
       call bessels2rr(F,LLMAX,mu,H, h_saved)
    end subroutine


    function bessel_funcs(l,x) result(sbj)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=ikind), intent(in):: l
        real(kind=dp),dimension(:), intent(in):: x
        real(kind=dp), dimension(size(x)):: sbj
        sbj=0.0_dp
        select case (l)
            case (0)
                sbj = dsin(x)/x
            case (1)
                sbj = (-(x*dcos(x)) + dsin(x))/x**2
            case (2)
                sbj = -((3*x*dcos(x) + (-3 + x**2)*dsin(x))/x**3)
            case (3)
                sbj = (x*(-15 + x**2)*dcos(x) + 3*(5 - 2*x**2)*dsin(x))/x**4
            case (4)
                sbj = (5*x*(-21 + 2*x**2)*dcos(x) + (105 - 45*x**2 &
                     + x**4)*dsin(x))/x**5
            case (5)
                sbj = (-(x*(945 - 105*x**2 + x**4)*dcos(x)) + 15*(63 - 28*x**2 &
                     + x**4)*dsin(x))/x**6
            case (6)
                sbj = -((21*x*(495 - 60*x**2 + x**4)*dcos(x) + (-10395 &
                     + 4725*x**2 - 210*x**4 &
                     + x**6)*dsin(x))/x**7)
            case (7)
                sbj = (x*(-135135 + 17325*x**2 - 378*x**4 + x**6)*dcos(x) - 7*(-19305 &
                     + 8910*x**2 - 450*x**4 &
                     + 4*x**6)*dsin(x))/x**8
            case (8)
                sbj = (9*x*(-225225 + 30030*x**2 - 770*x**4 + 4*x**6)*dcos(x) &
                     + (2027025 - 945945*x**2 + 51975*x**4 - 630*x**6 &
                     + x**8)*dsin(x))/x**9
            case (9)
                sbj = (-(x*(34459425 - 4729725*x**2 + 135135*x**4 - 990*x**6 &
                     + x**8)*dcos(x)) + 45*(765765 - 360360*x**2 + 21021*x**4 - 308*x**6 &
                     + x**8)*dsin(x))/x**10
            case (10)
                sbj = -((55*x*(11904165 - 1670760*x**2 + 51597*x**4 - 468*x**6 &
                    + x**8)*dcos(x) + (-654729075 + 310134825*x**2 - 18918900*x**4 &
                    + 315315*x**6 - 1485*x**8 & + x**10)*dsin(x))/x**11)
            case (11)
                sbj = (x*(-13749310575 + 1964187225*x**2 - 64324260*x**4 &
                     + 675675*x**6 - 2145*x**8 + x**10)*dcos(x) - 33*(-416645775 + 198402750*x**2 - 12530700*x**4 &
                     + 229320*x**6 - 1365*x**8 &
                     + 2*x**10)*dsin(x))/x**12
            case (12)
                sbj = (39*x*(-8108567775 + 1175154750*x**2 - 40291020*x**4 &
                     + 471240*x**6 - 1925*x**8 + 2*x**10)*dcos(x) + (316234143225 - 151242416325*x**2 &
                     + 9820936125*x**4 - 192972780*x**6 + 1351350*x**8 - 3003*x**10 &
                     + x**12)*dsin(x))/x**13
            case (13)
                sbj = (-(x*(7905853580625 - 1159525191825*x**2 &
                    + 41247931725*x**4 - 523783260*x**6 + 2552550*x**8 - 4095*x**10 + x**12)*dcos(x)) &
                     + 91*(86877511875 - 41701205700*x**2 + 2770007625*x**4 - 57558600*x**6 &
                     + 454410*x**8 - 1320*x**10 &
                     + x**12)*dsin(x))/x**14
            case (14)
                sbj = -((105*x*(2032933777875 - 301175374500*x**2 &
                     + 11043097065*x**4 - 149652360*x**6 + 831402*x**8 - 1768*x**10 + x**12)*dcos(x) &
                     + (-213458046676875 + 102776096548125*x**2 - 6957151150950*x**4 &
                     + 151242416325*x**6 - 1309458150*x**8 + 4594590*x**10 - 5460*x**12 &
                     + x**14)*dsin(x))/x**15)
            case (15)
                sbj = (x*(-6190283353629375 + 924984868933125*x**2 - 34785755754750*x**4 &
                     + 496939367925*x**6 - 3055402350*x**8 + 7936110*x**10 - 7140*x**12 &
                     + x**14)*dcos(x) - 15*(-412685556908625 + 2*x**2*(99613755115875 - 6851739769875*x**2 &
                    + 154603358910*x**4 - 1440403965*x**6 + 5819814*x**8 - 9282*x**10 &
                    + 4*x**12))*dsin(x))/x**16
            case (16)
                sbj = (17*x*(-11288163762500625 + 2*x**2*(849646734811875 - 32646524785875*x**2 &
                     + 487195458750*x**4 - 3247969725*x**6 + 9803430*x**8 - 11970*x**10 &
                     + 4*x**12))*dcos(x) + (191898783962510625 - 92854250304440625*x**2 &
                     + 6474894082531875*x**4 - 150738274937250*x**6 &
                    + 1490818103775*x**8 - 6721885170*x**10 + 13226850*x**12 - 9180*x**14 &
                    + x**16)*dsin(x))/x**17
            case (17)
                sbj = (-(x*(6332659870762850625 - 959493919812553125*x**2 &
                    + 37554385678684875*x**4 - 581419060472250*x**6 + 4141161399375*x**8 - 14054850810*x**10 &
                     + 21366450*x**12 - 11628*x**14 + x**16)*dcos(x)) &
                     + 153*(41389933795835625 - 20067846688890000*x**2 + 1416077891353125*x**4 - 33855655333500*x**6 &
                     + 351863386875*x**8 - 1732250520*x**10 + 3993990*x**12 - 3800*x**14 &
                     + x**16)*dsin(x))/x**18
!            case (18)
!                sbj = -((171*x*(1296158453080115625 - 197509859516970000*x**2 &
!                     + 7855505776243125*x**4 - 125495023153500*x**6 + 944475406875*x**8 - 3522519000*x**10 &
!                     + 6322470*x**12 - 4760*x**14 + x**16)*dcos(x) &
!                     + (-221643095476699771875 &
!                     + 107655217802968460625*x**2 - 7675951358500425000*x**4 &
!                     + 187771928393424375*x**6 - 2034966711652875*x**8 + 10767019638375*x**10 - 28109701620*x**12 &
!                     + 33575850*x**14 - 14535*x**16 &
!                     + x**18)*dsin(x))/x**19)
!            case (19)
!                sbj = (x*(-8200794532637891559375 &
!                     + 1255977541034632040625*x**2 - 50661278966102805000*x**4 &
!                     + 831561397170879375*x**6 - 6557114959770375*x**8 &
!                     + 26428139112375*x**10 - 54057118500*x**12 + 51482970*x**14 - 17955*x**16 &
!                     + x**18)*dcos(x) - 95*(-86324152975135700625 &
!                     + 41995533879795746250*x**2 - 3021900850609641000*x**4 + 75412855451934000*x**6 - 847091406286125*x**8 &
!                     + 4760156050650*x**10 - 13737824100*x**12 + 19509336*x**14 - 11781*x**16 &
!                     + 2*x**18)*dsin(x))/x**20
!            case (20)
!                sbj = (105*x*(-3046009397836931150625 &
!                     + 468616830436450946250*x**2 - 19138705387194393000*x**4 + 321658914070494000*x**6 - 2639877451336125*x**8 &
!                     + 11354311618650*x**10 - 25815032100*x**12 + 29418840*x**14 - 14421*x**16 &
!                     + 2*x**18)*dcos(x) + (319830986772877770815625 - 155815096120119939628125*x**2 &
!                     + 11303797869311688365625*x**4 - 287080580807915895000*x**6 &
!                    + 3326245588683517500*x**8 - 19671344879311125*x**10 &
!                     + 61665657928875*x**12 - 100391791500*x**14 + 77224455*x**16 - 21945*x**18 &
!                     + x**20)*dsin(x))/x**21
!            case (21)
!                sbj = (-(x*(13113070457687988603440625 - 2025596249561559215165625*x**2 &
!                     + 83648104232906493905625*x**4 - 1435402904039579475000*x**6 &
!                     + 12196233825172897500*x**8 - 55437426478058625*x**10 &
!                    + 137561852302875*x**12 - 180705224700*x**14 + 113565375*x**16 - 26565*x**18 + x**20)*dcos(x)) &
!                     + 231*(56766538777870080534375 - 27690994525790283187500*x**2 &
!                    + 2023572676884674540625*x**4 - 52196469237802890000*x**6 &
!                     + 621386538545272500*x**8 - 3839821747398000*x**10 + 12902626839375*x**12 - 23468211000*x**14 &
!                    + 21729825*x**16 - 8740*x**18 &
!                     + x**20)*dsin(x))/x**22
!            case (22)
!                sbj = -((253*x*(2228703674626812292284375 - 345535453430513533687500*x**2 &
!                     + 14411356716248247380625*x**4 - 251904699365048730000*x**6 &
!                    + 2206372491936112500*x**8 - 10517772612438000*x**10 &
!                    + 28092341379375*x**12 - 41426494200*x**14 + 31511025*x**16 - 10500*x**18 + x**20)*dcos(x) &
!                     + (-563862029680583509947946875 + 275374479611447760672253125*x**2 - 20255962495615592151656250*x**4 &
!                    + 529771326808407794735625*x**6 - 6459313068178107637500*x**8 &
!                     + 41467195005587851500*x**10 - 147833137274823000*x**12 &
!                     + 294775397791875*x**14 - 316234143225*x**16 + 164038875*x**18 - 31878*x**20 &
!        &            + x**22)*dsin(x))/x**23)
!            case (23)
!                sbj = (x*(-25373791335626257947657609375 &
!                    + 3947034207764084569635628125*x**2 - 166098892464047855643581250*x**4 &
!                    + 2951583106503986284955625*x**6 - 26554953724732220287500*x**8 + 131941075017779527500*x**10 - 375268733082243000*x**12 &
!                     + 609202488769875*x**14 - 539458244325*x**16 &
!                     + 233107875*x**18 - 37950*x**20 + x**22)*dcos(x) - 69*(-367736106313424028226921875 &
!                     + 179782096419896191577606250*x**2 - 13303114957074771046968750*x**4 &
!                     + 352277608619401602637500*x**6 - 4387340180607932047500*x**8 &
!                     + 29124116893556685000*x**10 - 109267971029217000*x**12 &
!                     + 235440575370000*x**14 - 284807147625*x**16 + 179729550*x**18 - 50050*x**20 &
!                    + 4*x**22)*dsin(x))/x**24
!            case (24)
!                sbj = (75*x*(-15900909236992454980532101875 &
!                     + 2480992930594567443770966250*x**2 - 105254245540375588523616750*x**4 &
!                     + 1898273056731975493069500*x**6 - 17490862853356955762700*x**8 &
!                     + 90125903550606323400*x**10 - 270648359010829800*x**12 &
!                    + 476531724548880*x**14 - 477805873545*x**16 &
!                     + 252378126*x**18 - 59202*x**20 + 4*x**22)*dcos(x) &
!                     + (1192568192774434123539907640625 - 583597200719403932796125015625*x**2 &
!                     + 43417376285404930265991909375*x**4 - 1162692247248334989505068750*x**6 &
!                     + 14757915532519931424778125*x**8 - 100908824153982437092500*x**10 &
!                     + 395823225053338582500*x**12 - 911366923199733000*x**14 + 1218404977539750*x**16 - 899097073875*x**18 &
!                     + 326351025*x**20 - 44850*x**22 &
!                     + x**24)*dsin(x))/x**25
!            case (25)
!                sbj = (-(x*(58435841445947272053455474390625 - 9143022811270661613805958578125*x**2 &
!                     + 390756386568644372393927184375*x**4 - 7142252375954057792673993750*x**6 &
!                     + 67230504092590798712878125*x**8 - 357767649273210458782500*x**10 &
!                     + 1126573794382579042500*x**12 - 2126522820799377000*x**14 &
!                    + 2365139074047750*x**16 - 1466947857375*x**18 + 450675225*x**20 - 52650*x**22 &
!                    + x**24)*dcos(x)) + 325*(179802589064453144779862998125 - 88066574235650519892177795000*x**2 &
!                    + 6584173546577890523853718125*x**4 - 178122569376020226732274500*x**6 & + 2299830818732970308911125*x**8 - 16145411864637189934800*x**10 &
!                     + 65861237210058467100*x**12 - 160604520731701200*x**14 &
!                     + 233683826461470*x**16 - 196022922480*x**18 + 87361659*x**20 - 17388*x**22 &
!                    + x**24)*dsin(x))/x**26
!            case (26)
!                sbj = -((351*x*(8490677816932509614604641578125 - 1331871030107060331702688875000*x**2 &
!                     + 57306695683177936040949028125*x**4 - 1060253389142977540073062500*x**6 &
!                     + 10174148683695239020903125*x**8 - 55720697512636766610000*x**10 &
!                     + 182947881139051297500*x**12 - 366812794263762000*x**14 &
!                     + 445475704038750*x**16 - 315241542000*x**18 + 119409675*x**20 - 20300*x**22 &
!                    + x**24)*dcos(x) + (-2980227913743310874726229193921875 &
!                     + 1460896036148681801336386859765625*x**2 - 109716273735247939365671502937500*x**4 &
!                     + 2995798963692940188353441746875*x**6 - 39282388067747317859706965625*x**8 &
!                     + 282368117188881354594088125*x**10 - 1192558830910701529275000*x**12 &
!                     + 3057843156181285972500*x**14 - 4784676346798598250*x**16 &
!                     + 4467484917645750*x**18 - 2347116571800*x**20 + 614557125*x**22 - 61425*x**24 &
!                    + x**26)*dsin(x))/x**27)
end select




    end function bessel_funcs


    subroutine integral_k_ijkr_alig(q_al,exponent1,l,m,n,group_start,group_count,hx,hy,hz,h,dx1,dy1,dz1, gi,gj,&
            za,&
            cutoff1, cutoff2,f)



        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        !INTEGER, PARAMETER :: dp = selected_real_kind(2*precision(1.0_dp))
        integer(kind=selected_int_kind(8)), intent(in)  :: gi,gj
        integer(kind=selected_int_kind(8)), dimension(:), intent(in) :: l,m,n,group_start,group_count
        real(kind=dp), intent(in)              :: cutoff1,cutoff2, hx, hy, hz, h
        real(kind=dp), intent(in), dimension(:,:) :: za
        real(kind=dp), intent(in), dimension(:,:)   ::  q_al

        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1
        complex(kind=dp), intent(in),  dimension(:,:,:,:) :: exponent1
        !real(kind=dp), pointer,dimension(:,:,:) :: dxx,dyy,dzz, dxx2,dyy2,dzz2



        integer(kind=selected_int_kind(8)) :: ka, posi, posj, posk, posr, ra,i,j,k,r
        integer(kind=selected_int_kind(8)) ::  h1
        integer(kind=selected_int_kind(8)) ::  ll, mm, nn, llp, mmp, nnp,llmax

        real(kind=dp)  ::coeff,prodd,ztot,mdn, mdl, mdm, mdlp, mdmp,mdnp,z1
        INTEGER(kind=ikind), parameter :: dim = 13
        real(kind=dp) :: prod1,prod2,prod3,prod4,prod5,prod6
        REAL(kind=dp), DIMENSION(dim,dim)           :: a, b, c
        complex(kind=dp), DIMENSION(size(q_al(1,:)))               :: h_saved
        REAL(kind=dp), DIMENSION(dim, dim, dim, dim)  :: h_pre2
        REAL(kind=dp), DIMENSION(dim)               :: BD



      !  real(kind=dp), dimension(:), allocatable :: pmu, h_0, h_1, h_r, muoh,zij, zij2, zkr, zkr2

        complex(kind=dp), intent(out), dimension(size(q_al(1,:))) :: f
        real(kind=dp), allocatable,dimension(:) :: l1vec,l2vec,suml1l2








        h_saved=(0.0_dp,0.0_dp)

! loop through all possible ways to get total angular momentum lmax1


        posI=1


        ! loop through all possible ways to get total angular momentum lmax1
        do i = group_start(gi), group_start(gi) + group_count(gi) - 1
            posj=1
            do j = group_start(gj), group_start(gj) + group_count(gj) - 1

                ztot=Za(posi,posj)


!
                        ! continue only if larger
                     !   if (abs(ztot) < cutoff1) then
                      !          posj=posj+1
                       !         cycle
                      !      end if
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
                                    H1=1.0
                                    MDN=Dz1(nn+1,n(i)+1,n(j)+1)
                                    if (abs(MDN)<1.0e-30) cycle ! check if MD coeff is 0
                                    prod3 = MDN  * prod2

                                    if (abs(prod3)>1E-30) then
                                    ! MD coeff 4
                                         h_saved = h_saved + exponent1(:,ll+1,mm+1,nn+1)*prod3

                                    end if

                                end do
                            end do
                           end do

                    posj=posj+1
            end do
        posi=posi+1
        end do

      F=h_saved*exp((0,1.00_dp)*q_al(1,:)*(Hx))*exp((0,1.0_dp)*q_al(2,:) &
              *(Hy))*exp((0,1.0_dp)*q_al(3,:)*(Hz))

    end subroutine

    subroutine integral_k_ij(mu,l,m,n,group_start,group_count,hx,hy,hz,h,dx1,dy1,dz1, gi,gj,&
            za,cutoff1, cutoff2,f)



        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        !INTEGER, PARAMETER :: dp = selected_real_kind(2*precision(1.0_dp))
        integer(kind=selected_int_kind(8)), intent(in)  :: gi,gj
        integer(kind=selected_int_kind(8)), dimension(:), intent(in) :: l,m,n,group_start,group_count
        real(kind=dp), intent(in)              :: cutoff1,cutoff2, hx, hy, hz, h
        real(kind=dp), intent(in), dimension(:,:) :: za
        real(kind=dp), intent(in), dimension(:)   ::  mu

        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1
        !real(kind=dp), pointer,dimension(:,:,:) :: dxx,dyy,dzz, dxx2,dyy2,dzz2



        integer(kind=selected_int_kind(8)) :: ka, posi, posj, ra,i,j,k,r
        integer(kind=selected_int_kind(8)) ::  h1
        integer(kind=selected_int_kind(8)) ::  ll, mm, nn, llp, mmp, nnp,llmax

        real(kind=dp)  ::coeff,prodd,ztot,mdn, mdl, mdm, mdlp, mdmp,mdnp,z1
        INTEGER(kind=ikind), parameter :: dim = 13
        real(kind=dp) :: prod1,prod2,prod3
        REAL(kind=dp), DIMENSION(dim,dim)           :: a, b, c
        REAL(kind=dp), DIMENSION(dim)               :: h_saved
        REAL(kind=dp), DIMENSION(dim, dim, dim, dim)  :: h_pre2
        REAL(kind=dp), DIMENSION(dim)               :: BD
        real(kind=dp),dimension(:), allocatable :: bess



      !  real(kind=dp), dimension(:), allocatable :: pmu, h_0, h_1, h_r, muoh,zij, zij2, zkr, zkr2

        real(kind=dp), intent(out), dimension(size(mu)) :: f
        real(kind=dp), allocatable,dimension(:) :: l1vec,l2vec,suml1l2


        LLmax = l(group_start(gi)) + m(group_start(gi)) + n(group_start(gi)) + &
                l(group_start(gj)) + m(group_start(gj)) + n(group_start(gj))


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


                    ztot=Za(posi,posj)


!
                        ! continue only if larger
!                        if (abs(ztot) < cutoff1) then
!                                posj=posj+1
!                                cycle
!                            end if
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
                                    if (abs(prod3)>cutoff2) then
                                                    ! add the contribution to the total
                                        h_saved = h_saved + h_pre2(:,ll+1,mm+1,nn+1)*prod3

                                    end if

                                    end do
                                    end do
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
        coeff=h_saved(2)
        h_sum=h_sum+coeff*h_1
        do ra = 2, LLmax
            coeff=h_saved(ra+1)
            h_r= ((2*ra-1)/(Pmu)*h_1-h_0*muOH)*muOH
            h_sum=h_sum+h_r*coeff
            h_0=h_1
            h_1=h_r
        end do
    end if

    END SUBROUTINE BesselSum

    SUBROUTINE BesselDeriv2j(BD, hz,h,LL,MM,NN, a, b, c, ap,bp,cp,LLmax)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        ! The three nested loops give the formula
        ! in the form of coefficients multiplying the h functions
        REAL(kind=dp), INTENT(out), DIMENSION(20)  :: BD
        INTEGER(kind=ikind), INTENT(in)                 :: LL, MM, NN, LLmax
        REAL(kind=dp), INTENT(in), DIMENSION(20,20)  :: a, b, c,ap,bp,cp
        real(kind=dp),intent(in) :: h,hz
        ! loop and temp variables
        INTEGER(kind=ikind) :: ii, jj, kk, hOrder, temp, ceil
        REAL(kind=dp)       :: C1, C2, C3, Ct2, Ct3,eta,C1p,Ct2p,Ct3p,pi



        pi=dacos(-1.0_dp)
        eta=3.0_dp*hz**2-h**2
        ! set this to 0 initially and accumulate
        BD = 0.0_dp
        !TempBesselPre=zeros(LLmax+1,1);
        do ii = 0, LL+2
            C1=a(LL+1,ii+1)
            C1p=ap(LL+1,ii+1)

           ! if (abs(C1)<1.0E-30) cycle
           ! if (abs(C1p)<1.0E-30) cycle

            do jj = 0, MM+2
                Ct2=b(MM+1,jj+1)
                Ct2p=bp(MM+1,jj+1)

               ! if (abs(Ct2)<1.0E-30) cycle
                !if (abs(Ct2p)<1.0E-30) cycle


                do kk = 0, NN+2
                    Ct3=c(NN+1,kk+1)
                    Ct3p=cp(NN+1,kk+1)
                  !  if (abs(Ct3)<1.0E-30) cycle
                   ! if (abs(Ct3p)<1.0E-30) cycle

                    C3 = c1*ct2*ct3*eta+c1p*ct2*ct3+c1*ct2p*ct3+c1*ct2*ct3p
                    ! hOrder = ceiling((LL+MM+NN-ii-jj-kk)/2.0_dp)+ii+jj+kk
                    temp = CEILING((LL+MM+NN+ii+jj+kk)/2.0_dp-1.0_dp)
                   ! ceil = temp/2_ikind + mod(temp, 2_ikind) ! integer division with rounding up
                    hOrder = temp
                    ! Barrilero style:
                    ! hOrder = ceiling((LL+MM+NN+ii+jj+kk)/2_ikind)
                    BD(hOrder+1)=BD(hOrder+1) + C3
!                    if (LL==2 .and. MM==0 .and. NN==0) then
!                        if (horder==3 .and. ii==2 .and. jj==2 .and. kk==2) then
!                                print*,ap(2,2)
!
!                            end if
!                    end if


                  !  print*,c3
                end do
            end do
        end do

         BD=0.25_dp * (5.0_dp/pi)**(0.5_dp) * BD


    END SUBROUTINE

    Subroutine bessels2(sum,order,mu,H, h_saved)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=ikind)::i,beta
        integer(kind=ikind), intent(in) :: order
        real (kind=dp), dimension(:), intent(in) :: mu, h_saved
        real (kind=dp),  intent(in) :: H
        real(kind=dp):: qrad,hqr,hfunc2
        real(kind=dp),dimension(size(mu)):: pmu, muOH,h1,h0,h2
        real(kind=dp), dimension(size(mu)), intent(out) :: sum

         sum=0.0_dp
       do i=1,size(mu)
        do beta=0,order+6

            qrad=H * mu(i)
            hqr = mu(i)**(beta-2) / h**beta
            hfunc2 = hqr * spherical_bessel_jn(beta,qrad)
            sum(i)=sum(i)+hfunc2*h_saved(beta+1)

        end do
      enddo
     End Subroutine bessels2



     Subroutine bessels2rr(sum,order,mu,H, h_saved)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer:: order2
        integer(kind=ikind)::i,beta,ra,maxord,ind
        integer(kind=ikind), intent(in) :: order
        real (kind=dp), dimension(:), intent(in) :: mu, h_saved
        real (kind=dp),  intent(in) :: H
        real(kind=dp):: qrad,hfunc2,sumpmu,coeff,bess1,bess2,bessnew,time1,time2
        real(kind=dp),dimension(size(mu)):: pmu, muOH,h_1,h_0,h_r,hqr,bess
        real(kind=dp), dimension(size(mu)), intent(out) :: sum
        real(kind=dp), dimension(0:order+6, size(mu)):: allbessels, allbessels2
        allbessels=0.0_dp
            sum=0.0_dp
            Pmu=H * mu
        print*,'calling new function'
        call cpu_time(time1)
          do i=1,size(mu)
            call SPHERICAL (pmu(i),5,allbessels(:,i),allbessels2(:,i),0.001_dp,order2)

          end do
        call cpu_time(time2)
        print*,time2-time1,order2
        call cpu_time(time1)
                 maxord=(maxval(pack([(ind, ind=1,size(h_saved))],h_saved/=0.0_dp)))
         do i=1,size(mu)
            allbessels(:,i)= spherical_bessel_jn_2(maxord-1,pmu(i))

         enddo
        call cpu_time(time2)
         print*,time2-time1
        stop

        do beta=0,order+6

                hqr = mu**(beta-2) / h**beta
                bess= allbessels(beta,:)

                sum=sum+bess*hqr*h_saved(beta+1)



         enddo





     End Subroutine bessels2rr






    Subroutine bessels2rr2(sum,order,mu,H, h_saved)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=ikind)::i,beta,ra
        integer(kind=ikind), intent(in) :: order
        real (kind=dp), dimension(:), intent(in) :: mu, h_saved
        real (kind=dp),  intent(in) :: H
        real(kind=dp):: qrad,hfunc2,sumpmu,coeff
        real(kind=dp),dimension(size(mu)):: pmu, muOH,h_1,h_0,h_r,hqr,bess1,bess2,bessnew
        real(kind=dp), dimension(size(mu)), intent(out) :: sum

            sum=0.0_dp
            Pmu=H * mu
      do i=1,size(mu)
            beta=order+6
            hqr(i) = mu(i)**(beta-2) / h**beta
            call sphbes(beta,pmu(i),bess2(i))

            sum(i)=sum(i)+bess2(i)*hqr(i)*h_saved(beta+1)

            beta=order+5
            hqr = mu**(beta-2) / h**beta
            call sphbes(beta,pmu(i),bess1(i))
            sum(i)=sum(i)+bess1(i)*hqr(i)*h_saved(beta+1)

      end do

            do ra=order+4,0,-1
                bessnew=bess1*(2*ra+3)/pmu-bess2
                hqr=mu**(ra-2) / h**ra
                sum=sum+bessnew*hqr*h_saved(ra+1)
                bess2=bess1
                bess1=bessnew
            end do





     End Subroutine bessels2rr2






    SUBROUTINE rrdj2a(lmax,x,a)
    ! This subroutine uses a recurrence relation to calculate the
    ! first set of expansion coefficients, a, at given coordinate
    ! x up to an angular momentum of L
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer, parameter    :: dp = selected_real_kind(15)
        integer(kind=ikind), intent(in)                                        :: lmax
        integer(kind=ikind)                                                    :: l, p
        real(kind=dp), intent(in)                                  :: x
        real(kind=dp), dimension(:,:), intent(out)   :: a

        a=0.0_dp
        a(1,3) = 1.0_dp
        a(2,4) = -x

        do l = 2, lmax
            a(l+1,3) = -(l-1) * a(l-1,3)
        enddo

        do l = 2, lmax
            do p = 1, l
                a(l+1,p+3) = -(l-1) * a(l-1,p+3) - x * a(l,p+2)
            enddo
        enddo

    END SUBROUTINE

     SUBROUTINE rrdj2b(lmax,a,b)
    ! This subroutine uses the first set of of expansion coefficients,
    ! a, to calculate the second set, b.
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer, parameter    :: dp = selected_real_kind(15)
        integer, intent(in)                                                   :: lmax
        integer                                                    :: l, p
        real(kind=dp), dimension(:,:), intent(in)                  :: a
        real(kind=dp), dimension(:,:), intent(out)   :: b
        real(kind=dp)                                              :: carg

        b=0.0_dp
       ! lmax = size(a,1)-1



        do l = 1, lmax
            do p = 0, l
                carg = (l+p)/2.0_dp
                b(l+1,p+1) = 2.0_dp * ceiling(carg) * a(l+1,p+3)
            enddo
        enddo

    END SUBROUTINE



    Subroutine bessels(sum,order,mu,H, h_saved)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=ikind)::i
        integer(kind=ikind), intent(in) :: order
        real (kind=dp), dimension(:), intent(in) :: mu, h_saved
        real (kind=dp),  intent(in) :: H
        real(kind=dp),dimension(size(mu)):: pmu, muOH,h1,h0,h2
        real(kind=dp), dimension(size(mu)), intent(out) :: sum



            Pmu=H * mu
            sum=0.0_dp
            h0=sin(Pmu)/Pmu
            sum=sum+h0*h_saved(1)


            if (order>=1) then
               muOH=mu/H
               h1=(sin(Pmu)/Pmu**2-cos(Pmu)/Pmu)
               sum=sum+h1*muOh*h_saved(2)


               do i=2,order

                   h2=(h0*((2*i-1)/Pmu-h1))
                   sum=sum+h2*muOh**(i)*h_saved(i+1)
                   h0=h1
                   h1=h2
               end do

            end if

     End Subroutine bessels


    SUBROUTINE BesselMats(h_sum, mu, H, LLmax, h_saved)
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
        coeff=h_saved(2)
        h_sum=h_sum+coeff*h_1
        do ra = 2, LLmax
            coeff=h_saved(ra+1)
            h_r= ((2*ra-1)/(Pmu)*h_1-h_0*muOH)*muOH
            h_sum=h_sum+h_r*coeff
            h_0=h_1
            h_1=h_r
        end do
    end if

    END SUBROUTINE BesselMats



        SUBROUTINE Hermite_like_coeffs(a, LLmax, Hx)
                  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        REAL(kind=dp), INTENT(out), DIMENSION(13,13) :: a
        INTEGER(kind=ikind), INTENT(in)     :: LLmax
        REAL(kind=dp), INTENT(in)           :: Hx

        ! loop vars
        INTEGER(kind=ikind)     :: L, ka
        a = 0.0_dp
        a(1,1)=1;               ! L=0

        if (LLmax>0) a(2,2)=-Hx ! L=1

        do L = 2, LLmax
            a(L+1,1)=-(L-1)*a(L-1,1) ! L>1, ka = 0

            do ka=1, LLmax           ! L>1, ka > 0
                a(L+1,ka+1)=-Hx*a(L,ka)-(L-1)*a(L-1,ka+1);
            end do
        end do
    END SUBROUTINE Hermite_like_coeffs

    SUBROUTINE rrdj0(lmax,x,a)
    ! This subroutine uses a recurrence relation to calculate the
    ! first set of expansion coefficients, a, at given coordinate
    ! x (or y or z) up to an angular momentum of L

        integer, parameter    :: dp = selected_real_kind(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer, intent(in)                                        :: lmax
        integer(kind=ikind)                                                    :: l,p
        real(kind=dp), intent(in)                                  :: x
        real(kind=dp), dimension(13,13), intent(out)   :: a

        a=0.0_dp

        a(1,1) = 1.0_dp
        a(2,2) = -x

        do l = 2, lmax
            a(l+1,1) = -(l-1) * a(l-1,1)
        enddo

        do l = 2, lmax
            do p = 1, l
                a(l+1,p+1) = -(l-1) * a(l-1,p+1) - x * a(l,p)
            enddo
        enddo

    END SUBROUTINE







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