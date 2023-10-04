Module TSj0groupsfast
    use Bessels_j0
    use p0_cases
    use twordms
    use MD
    use twordms
    use onerdm
    implicit none
contains

    subroutine total_scattering_j0_groups_fast(q, l, m, n, ngto, ng, nq, maxl, typec, state1, &
            state2, ncontr, group, gs, gf, gc, confs, ga, xx, yy, zz, coeffs, mmod, civecs, geom, &
            cutoffmd, cutoffz, cutoffcentre,contrvec, result2)
        implicit none
        !Precision parameters
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

        !Entering variables from readers
        integer(kind = ikind), intent(in) :: ngto, ng, nq, maxl, typec, state1, state2, ncontr
        integer(kind = ikind), intent(in), dimension(:), allocatable :: l, m, n, group, gs, gf, gc,contrvec
        integer(kind = ikind), dimension(:, :), allocatable, intent(in) :: confs
        real(kind = dp), intent(in), dimension(:), allocatable :: ga, xx, yy, zz, q, coeffs
        real(kind = dp), intent(in), dimension(:, :), allocatable :: mmod, civecs, geom
        real(kind = dp), intent(in) :: cutoffmd, cutoffz, cutoffcentre

        !twordm variables
        integer(kind = ikind), dimension(:), allocatable :: m1, m2, m3, m4
        integer(kind = ikind), dimension(:, :), allocatable :: mat, ep3, ndiff2
        integer(kind = ikind) :: nmat, i, j,numberlines
        real(kind = dp), dimension(:), allocatable :: total

        !Variables to create total_variables
        real(kind = dp), DIMENSION(size(q), 4 * maxval(l) + 1, 4 * maxval(l) + 1, 4 * maxval(l) + 1) :: P0matrix
        real(kind = dp), dimension(maxl * 2 + 1, maxl + 1, maxl + 1, ng, ng) :: ddx, ddy, ddz
        real(kind = dp), dimension(:, :, :), allocatable :: z1, z2
        real(kind = dp), dimension(:,:),allocatable :: px, py, pz
        real(kind = dp), dimension(:, :, :),allocatable :: e12
        integer(kind = ikind), dimension(maxval(group)) :: group_start, group_count
        integer(kind = ikind), dimension(size(group)) :: group_sorted
         character(len=60):: File_in,File_out

        !Result out
        real(kind = dp), allocatable, dimension(:), intent(out) :: Result2
        group_sorted=group
        call Bubble_Sort(group_sorted)
        do i = 1, Ng
            do j = 1, Ngto
                if (group_sorted(j) == i) then
                    group_start(i) = j
                    group_count(i) = count(group==i)
                    exit
                end if
            end do
        end do

        P0matrix = 0.0_dp
        call set_P0(P0matrix, 4 * maxval(l), q)

       ! fileout_8='es.dat'
        open(unit=15, file='bitwise.dat')
              numberlines=0
              do while(.true.)
                  read (15, *, end=999) i
                  numberlines=numberlines+1
              end do
              

999 continue
        close(15)
        print*,numberlines
        file_in='bitwise.dat'
        file_out='es.dat'
        call mcci_to_bit(file_in,file_out,numberlines)

        call createtwordm_bit(file_out,numberlines,mat,total)

       ! call maxcoincidence(confs, ep3, ndiff2)

       !call createtwordm(confs, civecs, ndiff2, ep3, mat, total, state1, state2)

        call system('python3 Zcontr.py')


        call variables_total(px, py, pz, ddx, ddy, ddz, &
                e12, maxl, ngto, ng, group_start, group_count, group, ga, l, m, n, xx, yy, zz, &
                mmod,  q, nq)

        call tot_integration(ng, ncontr,px, py, pz, l, m, n, p0matrix, ddx, ddy, ddz, &
                group_start, group_count, group, &
                gs,gf,gc,contrvec,cutoffz, cutoffmd, cutoffcentre, q,coeffs, e12, result2)

    end subroutine


    subroutine variables_total(px, py, pz, ddx, ddy, ddz, &
            e12, maxl, ngto, ng, group_start, group_count, group, ga, l, m, n, xx, yy, zz, &
            mmod, q, nq)

        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind = ikind), intent(in) :: ngto, ng, nq, maxl
        integer(kind = ikind), intent(in), dimension(:), allocatable :: l, m, n, group
        integer(kind = ikind), intent(in), dimension(:) :: group_start, group_count

        real(kind = dp), intent(in), dimension(:), allocatable :: ga, xx, yy, zz, q
        real(kind = dp), dimension(ng) :: ga2, xx2, yy2, zz2
        real(kind = dp), intent(in), dimension(:, :), allocatable :: mmod
        real(kind = dp), intent(out), dimension(:, :), allocatable :: px, py, pz


        real(kind = dp), dimension(ngto, ngto) :: cmat
        !real(kind=dp), intent(out), dimension(ngto,ngto,ngto,ngto) :: zcontr

        real(kind = dp), intent(out), dimension(:,:,:), allocatable :: e12

        real(kind = dp), intent(out), dimension(maxl * 2 + 1, maxl + 1, maxl + 1, ng, ng) :: ddx, ddy, ddz

        integer(kind = ikind) :: max4l, max2l, i, j, k, N1, N2, nn1, nn2, ii, jj, ls, ms, ns, count
        integer(kind = ikind), dimension(ngto) :: ll

        ! real(kind = dp), dimension(:, :, :), intent(out), allocatable :: z11, z22

        real(kind = dp) :: pi, gap, time1, time2, time3, time4

        integer(kind = ikind), dimension(:), allocatable :: iduplicates, jduplicates, m11, m22, m33, m44
        real(kind = dp), dimension(ngto, ngto, maxl * 2 + 1, maxl + 1, maxl + 1) :: dx, dy, dz
        real(kind = dp), dimension(nq) :: preexp



        real(kind = dp) :: obwohl
        integer(kind = ikind) :: counter
        logical :: divided

        pi = acos(-1.0000)
        max4l = maxval(l) * 4
        max2l = maxval(l) * 2 + 1




        allocate(px(ng,ng), py(ng,ng), pz(ng,ng), e12(nq,ng,ng))
        gap = 0.0
        px = 0.0
        py = 0.0
        pz = 0.0
        e12 = 0.0
        ddx = 0.0
        ddy = 0.0
        ddz = 0.0
        dx = 0.0
        dy = 0.0
        dz = 0.0







        call fill_md_table(dx, l, xx, ga)
        call fill_md_table(dy, m, yy, ga)
        call fill_md_table(dz, n, zz, ga)




        do jj = 1, Ng
            ! essentially taking the values for the first GTO in the group.
            ! All gtos in the group have the same prefactors as the prefactors
            ! do not't depend on l, m and n
            j = group_start(jj)
            do ii = 1, Ng
                i = group_start(ii)
                gaP = ga(i) + ga(j)
                Px(ii, jj) = (ga(i) * xx(i) + ga(j) * xx(j)) / gaP
                Py(ii, jj) = (ga(i) * yy(i) + ga(j) * yy(j)) / gaP
                Pz(ii, jj) = (ga(i) * zz(i) + ga(j) * zz(j)) / gaP

                E12(:, ii, jj) = (pi / gaP)**1.5 * exp(-q * q * 0.25 / gaP) &
                        * exp(-ga(i) * ga(j) / gaP * ((xx(i) - xx(j))**2. + (yy(i) - yy(j))**2. + (zz(i) - zz(j))**2.))
            end do
        end do


        do jj = 1, Ngto
            j = group(jj)
            do ii = 1, Ngto
                i = group(ii)
                do ls = 1, l(ii) + l(jj) + 1
                    Ddx(ls, l(ii) + 1, l(jj) + 1, j, i) = Dx(ii, jj, ls, l(ii) + 1, l(jj) + 1)
                end do

                do ms = 1, m(ii) + m(jj) + 1
                    Ddy(ms, m(ii) + 1, m(jj) + 1, j, i) = Dy(ii, jj, ms, m(ii) + 1, m(jj) + 1)
                end do

                do ns = 1, n(ii) + n(jj) + 1
                    Ddz(ns, n(ii) + 1, n(jj) + 1, j, i) = Dz(ii, jj, ns, n(ii) + 1, n(jj) + 1)
                end do
            end do
        end do


        print*, 'leaving variables total'

    end subroutine variables_total

    subroutine tot_integration(ncap,ncontr, px, py, pz, l, m, n, p0matrix, dx, dy, dz,&
            group_start, group_count, group, &
            gs,gf,gc,contrvec,cutoffz, cutoffmd, cutoffcentre, q,coeff, e12, tsi)

        use omp_lib
        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

        real(kind = kind(1.d0)), external :: ddot
        INTEGER(kind = ikind), INTENT(IN) :: ncap,ncontr
        INTEGER(kind = ikind), INTENT(IN), DIMENSION(:) :: l, m, n, group, group_start, group_count
        INTEGER(kind = ikind), INTENT(IN), DIMENSION(:),allocatable:: gs,gf,gc,contrvec
        real(kind=dp), intent(in), dimension(:)::coeff
        real(kind=dp),  dimension(:),allocatable::c1,c2,c3,c4

        REAL(kind = dp), intent(in), dimension(:, :, :, :, :), target :: dx, dy, dz
        real(kind = dp), dimension(:, :, :), pointer :: dx1, dy1, dz1, dx2, dy2, dz2
        REAL(kind = dp), intent(in), dimension(:, :, :, :) :: p0matrix
        REAL(kind = dp), intent(in), dimension(:, :, :),allocatable :: e12
        !real(kind = dp), intent(in), dimension(:, :, :) :: z1, z2
        real(kind = dp), dimension(:, :, :, :), allocatable :: zcontrred, zcontrred2
        REAL(kind = dp), intent(in), dimension(:, :),allocatable :: px, py, pz
        REAL(kind = dp), intent(in), dimension(:), allocatable :: q
        REAL(kind = dp), intent(in) :: cutoffz, cutoffmd, cutoffcentre

        REAL(kind = dp), dimension(size(q)) :: f
        integer(kind = ikind), dimension(:, :), allocatable :: posits,vecgroup,vecgto
        real(kind = dp), dimension(:, :), allocatable :: za, zb, cmat
        integer(kind = ikind), dimension(:), allocatable :: posi, posj, posk, posr
        REAL(kind = dp), intent(out), dimension(:), allocatable :: tsi
        real(kind = dp) :: hx, hy, hz, h
        integer(kind = ikind), dimension(size(l)) :: ll
        integer(kind = ikind) :: nq, i, j, k, r, count1, ng, ii, jj,npos,iii,count2
        integer(kind = ikind) :: spi, spj, spk, spr, szo, nt, ngto,newcap
        real(kind = dp), dimension(:, :, :, :), allocatable :: Zbig
        integer(kind=ikind),dimension(size(l)) :: veccounts
        integer(kind=ikind), dimension(:), allocatable::group_pos,ngtovec,group_sorted,group_pos2
        integer(kind=ikind), dimension(:,:), allocatable:: bigvec

        nq = size(q)
        ng = maxval(group)
        ll = l + m + n
        ngto = size(l)
        allocate(ngtovec(ngto),posits(ng, maxval(group_count)),&
                tsi(nq),vecgroup(ng, maxval(group_count)), group_pos(maxval(group_count)),&
                vecgto(ng,maxval(group_count)), group_sorted(size(group)), group_pos2(maxval(group_count)))

        group_sorted=group
        call Bubble_Sort(group_sorted)
        posits = 1
        count1 = 1
        do i = 1, ncap
            count1 = 1
            do j = group_start(i), group_start(i) + group_count(i) - 1
                posits(i, count1) = j
                count1 = count1 + 1
            end do
        end do

        ! szo = size(z1(:, 1, 1))

        do i=1,ncontr
            veccounts(gs(i):gf(i))=i
        end do

        group_pos=0
        ngtovec=0
        print*,size(group_pos)
        do i=1,ngto
            ngtovec(i)=i

        end do
        print*,group
        print*,group_start

        do j=1,ng

            npos=COUNT(group==j)

            group_pos(1:npos)=pack(ngtovec, group==j)
            group_pos2(1:npos)=pack(ngtovec,group_sorted==j)
            if (group_count(j)/=npos) then
                print*,'ouch'
                stop
            end if
            do i=1,npos
                vecgroup(j,i)=contrvec(group_pos(i))
                vecgto(j,i)=group_pos(i)
            end do

        enddo

        newcap=CEILING(1.d0/24.d0*(ncap-2)*(ncap-1)*(3*ncap-1)*ncap)

        count2=1
        allocate(bigvec(4,newcap))
        do i=1,ncap
            do j=i+1,ncap
                do k=i+1,ncap
                    do r=k+1,ncap
                        bigvec(1,count2)=i
                        bigvec(2,count2)=j
                        bigvec(3,count2)=k
                        bigvec(4,count2)=r
                        count2=count2+1


                    enddo
                enddo
            enddo
        enddo


        print*, OMP_get_num_threads()
        call OMP_set_num_threads(16)
        print*, OMP_get_num_threads()
        !
        !First big loop

        tsi = 0.0_dp

        if (any(isnan(p0matrix))) print*, 'ouch'

        count1 = 0
        allocate(Zbig(ncontr,ncontr,ncontr,ncontr))

        open(40, file='Zcotr.dat', status='old', access='stream', form='unformatted')
        read(40) Zbig
        close(40)
        print*,'starting loop'
        !$OMP PARALLEL do private(posI,posK,posJ,posR,spi,spj,spk,spr,zcontrred), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,iii,i,j,k,r,dx1,dx2,dy1,dy2,dz1,dz2,c1,c2,c3,c4) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz, posits,cutoffmd,group_count,group_start) REDUCTION(+:tsi), &
        !$OMP& schedule(dynamic)
        do iii=1,newcap

            i=bigvec(1,iii)
            j=bigvec(2,iii)
            k=bigvec(3,iii)
            r=bigvec(4,iii)

            hx = px(k, r) - px(i, j)
            hy = py(k, r) - py(i, j)
            hz = pz(k, r) - pz(i, j)
            h = (hx * hx + hy * hy + hz * hz)**0.5

            allocate(posI(size(posits(i, :group_count(i)))), posJ(size(posits(j, :group_count(j)))) &
                    , posK(size(posits(k, :group_count(k)))), posR(size(posits(r, :group_count(r)))))


            posI = vecgroup(i, :group_count(i))
            posJ = vecgroup(j, :group_count(j))
            posK = vecgroup(k, :group_count(k))
            posR = vecgroup(r, :group_count(r))

            spi = size(posI)
            spj = size(posJ)
            spk = size(posK)
            spr = size(posR)


            allocate(zcontrred(spi, spj, spk, spr), c1(spi), c2(spj), c3(spk), c4(spr))
            Zcontrred=Zbig(posI,posJ,posK,posR)
            posI = vecgto(i, :group_count(i))
            posJ = vecgto(j, :group_count(j))
            posK = vecgto(k, :group_count(k))
            posR = vecgto(r, :group_count(r))
            c1=coeff(posI)
            c2=coeff(posJ)
            c3=coeff(posK)
            c4=coeff(posR)


            deallocate(posI, posJ, posK, posR)
            dx1 => dx(:, :, :, j, i)
            dy1 => dy(:, :, :, j, i)
            dz1 => dz(:, :, :, j, i)
            dx2 => dx(:, :, :, r, k)
            dy2 => dy(:, :, :, r, k)
            dz2 => dz(:, :, :, r, k)

            if (h < cutoffcentre) then

                call tot_integral_ijkr_pzero(nq, l, m, n, group_start, group_count, p0matrix, dx1, dy1, dz1, dx2, &
                        dy2, dz2, i, j, k, r, &
                        zcontrred, c1,c2,c3,c4,cutoffz, cutoffmd, f)

            else

                call tot_integral_k_ijkr(q, l, m, n, group_start, group_count, hx, hy, hz, h, dx1, dy1, dz1, dx2, &
                        dy2, dz2, i, j, k, r, &
                        zcontrred,c1,c2,c3,c4, cutoffz, cutoffmd, f)

            end if
            tsi = tsi + 8.000 * f * e12(:, i, j) * e12(:, k, r)
            count1 = count1 + 1

            deallocate(zcontrred,c1,c2,c3,c4)


        end do

        !$OMP END parallel DO

        print*,OMP_get_num_threads()
        print*, 'intermediate step', tsi(1), count1

        !$OMP PARALLEL do private(c1,c2,c4,posI,posJ,posR,spi,spj,spk,spr,zcontrred), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,j,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz,posits, cutoffmd,group_count,group_start) REDUCTION(+:tsi)
        do i = 1, ncap
            do j = i + 1, ncap
                do r = i + 1, ncap
                    hx = px(i, r) - px(i, j)
                    hy = py(i, r) - py(i, j)
                    hz = pz(i, r) - pz(i, j)
                    h = sqrt((hx * hx + hy * hy + hz * hz))
                    allocate(posI(size(posits(i, :group_count(i)))), posJ(size(posits(j, :group_count(j)))) &
                            , posR(size(posits(r, :group_count(r)))))

                    posI = vecgroup(i, :group_count(i))
                    posJ = vecgroup(j, :group_count(j))

                    posR = vecgroup(r, :group_count(r))

                    spi = size(posI)
                    spj = size(posJ)
                    spr = size(posR)

                    allocate(zcontrred(spi, spj, spi, spk), c1(spi), c2(spi), c4(spr))
                    Zcontrred=Zbig(posI,posJ,posI,posR)
                    posI = vecgto(i, :group_count(i))
                    posJ = vecgto(j, :group_count(j))
                    posR = vecgto(r, :group_count(r))
                    c1=coeff(posI)
                    c2=coeff(posJ)

                    c4=coeff(posR)


                    !
                    dx1 => dx(:, :, :, j, i)
                    dy1 => dy(:, :, :, j, i)
                    dz1 => dz(:, :, :, j, i)
                    dx2 => dx(:, :, :, r, i)
                    dy2 => dy(:, :, :, r, i)
                    dz2 => dz(:, :, :, r, i)

                    if (h < cutoffcentre) then
                        call tot_integral_ijkr_pzero(nq, l, m, n, group_start, group_count, p0matrix, dx1, dy1, dz1, dx2, &
                                dy2, dz2, i, j, i, r, &
                                zcontrred, c1,c2,c1,c4, cutoffz, cutoffmd, f)
                    else

                        call tot_integral_k_ijkr(q, l, m, n, group_start, group_count, hx, hy, hz, h, dx1, dy1, dz1, dx2, &
                                dy2, dz2, i, j, i, r, &
                                zcontrred, c1,c2,c1,c4, cutoffz, cutoffmd, f)

                    end if
                    tsi = tsi + 4.000 * f * e12(:, i, j) * e12(:, i, r)
                    count1= count1 + 1
                    deallocate(posI, posJ, posR)
                    deallocate(zcontrred,c1,c2,c4)


                end do
            end do

        end do
        !$OMP END parallel DO
        print*, 'intermediate step', tsi(1)
        !
        !$OMP PARALLEL do private(posI,posK,posR,spi,spj,spk,spr,zcontrred,c1,c3,c4), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,k,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz,posits, cutoffmd,group_count,group_start) REDUCTION(+:tsi)
        do i = 1, ncap
            do k = 1, ncap
                do r = k + 1, ncap
                    hx = px(k, r) - px(i, i)
                    hy = py(k, r) - py(i, i)
                    hz = pz(k, r) - pz(i, i)
                    h = sqrt((hx * hx + hy * hy + hz * hz))

                    allocate(posI(size(posits(i, :group_count(i)))) &
                            , posK(size(posits(k, :group_count(k)))), posR(size(posits(r, :group_count(r)))))


                    posI = vecgroup(i, :group_count(i))
                    posK = vecgroup(k, :group_count(k))
                    posR = vecgroup(r, :group_count(r))
                    spi = size(posI)
                    spk = size(posK)
                    spr = size(posR)

                    allocate(zcontrred(spi, spi, spk, spr),c1(spi),c3(spk),c4(spr) )
                    Zcontrred=Zbig(posI,posI,posK,posR)
                    posI = vecgto(i, :group_count(i))
                    posK = vecgto(k, :group_count(k))
                    posR = vecgto(r, :group_count(r))
                    c1=coeff(posI)
                    c3=coeff(posK)
                    c4=coeff(posR)




                    dx1 => dx(:, :, :, i, i)
                    dy1 => dy(:, :, :, i, i)
                    dz1 => dz(:, :, :, i, i)
                    dx2 => dx(:, :, :, r, k)
                    dy2 => dy(:, :, :, r, k)
                    dz2 => dz(:, :, :, r, k)

                    if (h < cutoffcentre) then
                        call tot_integral_ijkr_pzero(nq, l, m, n, group_start, group_count, p0matrix, dx1, dy1, dz1, dx2, &
                                dy2, dz2, i, i, k, r, &
                                zcontrred, c1,c1,c3,c4, cutoffz, cutoffmd, f)
                    else

                        call tot_integral_k_ijkr(q, l, m, n, group_start, group_count, hx, hy, hz, h, dx1, dy1, dz1, dx2, &
                                dy2, dz2, i, i, k, r, &
                                zcontrred, c1,c1,c3,c4, cutoffz, cutoffmd, f)

                    end if
                    tsi = tsi + 4.000 * f * e12(:, i, i) * e12(:, k, r)
                    count1 = count1 + 1
                    deallocate(posI, posK, posR)
                    deallocate(zcontrred,c1,c3,c4)


                end do
            end do
        end do

        !$OMP END parallel DO
        print*, 'intermediate step', tsi(1)
        !
        !$OMP PARALLEL do private(posI,posK,spi,spk,zcontrred,c1,c3), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,k,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz, cutoffmd,posits,group_count,group_start) REDUCTION(+:tsi)
        do i = 1, ncap
            do k = i + 1, ncap

                hx = px(k, k) - px(i, i)
                hy = py(k, k) - py(i, i)
                hz = pz(k, k) - pz(i, i)
                h = sqrt((hx * hx + hy * hy + hz * hz))

                allocate(posI(size(posits(i, :group_count(i)))), &
                        posK(size(posits(k, :group_count(k)))))
                posI = vecgroup(i, :group_count(i))
                posK = vecgroup(k, :group_count(k))


                spi = size(posI)
                spk = size(posK)

                allocate(zcontrred(spi, spi, spk, spk),c1(spi),c3(spk))
                Zcontrred=Zbig(posI,posI,posK,posK)
                posI = vecgto(i, :group_count(i))
                posK = vecgto(k, :group_count(k))

                c1=coeff(posI)
                c3=coeff(posK)

                dx1 => dx(:, :, :, i, i)
                dy1 => dy(:, :, :, i, i)
                dz1 => dz(:, :, :, i, i)
                dx2 => dx(:, :, :, k, k)
                dy2 => dy(:, :, :, k, k)
                dz2 => dz(:, :, :, k, k)

                if (h < cutoffcentre) then
                    call tot_integral_ijkr_pzero(nq, l, m, n, group_start, group_count, p0matrix, dx1, dy1, dz1, dx2, &
                            dy2, dz2, i, i, k, k, &
                            zcontrred, c1,c1,c3,c3, cutoffz, cutoffmd, f)
                else

                    call tot_integral_k_ijkr(q, l, m, n, group_start, group_count, hx, hy, hz, h, dx1, dy1, dz1, dx2, &
                            dy2, dz2, i, i, k, k, &
                            zcontrred, c1,c1,c3,c3, cutoffz, cutoffmd, f)
                end if
                tsi = tsi + 2.000 * f * e12(:, i, i) * e12(:, k, k)
                deallocate(posI, posK,c1,c3)
                deallocate(zcontrred)


            end do
        end do


        !$OMP END parallel DO
        print*,tsi(1)
        !
        !$OMP PARALLEL do private(posI,spi,zcontrred,c1), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,ll, p0matrix), &
        !$OMP& shared( cutoffz, cutoffmd,posits,group_count,group_start) REDUCTION(+:tsi)
        do i = 1, ncap
            allocate(posI(size(posits(i, :group_count(i)))))

            posI = vecgroup(i, :group_count(i))
            spi = size(posI)
            allocate(zcontrred(spi, spi, spi, spi),c1(spi))
            zcontrred=Zbig(posI,posI,posI,posI)
            posI = vecgto(i, :group_count(i))


            c1=coeff(posI)


            dx1 => dx(:, :, :, i, i)
            dy1 => dy(:, :, :, i, i)
            dz1 => dz(:, :, :, i, i)
            dx2 => dx(:, :, :, i, i)
            dy2 => dy(:, :, :, i, i)
            dz2 => dz(:, :, :, i, i)


            call tot_integral_ijkr_pzero(nq, l, m, n, group_start, group_count, p0matrix, dx1, dy1, dz1, dx2, &
                    dy2, dz2, i, i, i, i, &
                    zcontrred, c1,c1,c1,c1, cutoffz, cutoffmd, f)

            tsi = tsi + f * e12(:, i, i) * e12(:, i, i)
            count1 = count1 + 1
            deallocate(posI,c1)
            deallocate(zcontrred)

        end do
        !$OMP END parallel DO
        print*, tsi(1)

    end subroutine tot_integration


    subroutine tot_integral_ijkr_pzero(nq, l, m, n, gs, gc, p0mat, dx1, dy1, dz1, dx2, dy2, dz2, gi, gj, gk, gr, zcontr, &
            coeff1,coeff2,coeff3,coeff4,cutoff1, cutoff2, f)

        implicit none

        ! definition of input
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(kind = ikind), INTENT(IN) :: nq, gi, gj, gk, gr
        REAL(kind = dp), INTENT(IN) :: cutoff1, cutoff2
        REAL(kind = dp), INTENT(IN), DIMENSION(:, :, :, :) :: p0mat
        real(kind = dp), intent(in), dimension(:, :, :) :: dx1, dy1, dz1, dx2, dy2, dz2
        REAL(kind = dp), INTENT(IN), DIMENSION(:, :, :, :) :: zcontr
        real(kind=dp), intent(in), dimension(:)::coeff1,coeff2,coeff3,coeff4
        INTEGER(kind = ikind), INTENT(IN), DIMENSION(:) :: gs, gc, l, m, n
        ! definition of output
        REAL(kind = dp), INTENT(OUT), DIMENSION(nq) :: f
        ! definition of loop indices
        INTEGER(kind = ikind) :: h1
        INTEGER(kind = ikind) :: ll, mm, nn, llp, mmp, nnp
        ! definition of internal variables
        INTEGER(kind = ikind) :: i, j, k, r, posi, posj, posk, posr
        REAL(kind = dp) :: ztot
        real(dp), external :: ddot
        REAL(kind = dp) :: mdl, mdm, mdn, mdlp, mdmp, mdnp
        real(kind = dp) :: prod6, prod5, prod4, prod3, prod2, prod1
        !        REAL(kind=dp), DIMENSION(size(Z1(:,1,1)))                          :: zij1, zkr1
        !        REAL(kind=dp), DIMENSION(size(Z2(:,1,1)))                          :: zij2, zkr2

        posI = 1
        ! posI=apos(i)

        f = 0.0
        ! loop through all possible ways to get total angular momentum lmax1
        do i = gs(gi), gs(gi) + gc(gi) - 1
            posj = 1
            do j = gs(gj), gs(gj) + gc(gj) - 1
                posk = 1
                do k = gs(gk), gs(gk) + gc(gk) - 1
                    posr = 1
                    do r = gs(gr), gs(gr) + gc(gr) - 1
                        ztot = zcontr(posI, posJ, posK, posR)*coeff1(posI)*&
                                coeff2(posJ)*coeff3(posK)*coeff4(posR)
                        !

                        ! continue only if larger
                        if (abs(ztot) < cutoff1) then
                            posr = posr + 1
                            cycle
                        end if
                        do ll = 0, l(i) + l(j)
                            MDL = Dx1(ll + 1, l(i) + 1, l(j) + 1)
                            if (abs(MDL)<1.0e-30) cycle
                            prod1 = MDL * ztot
                            ! MD coeff 2
                            do mm = 0, m(i) + m(j)
                                MDM = Dy1(mm + 1, m(i) + 1, m(j) + 1)
                                if (abs(MDM)<1.0e-30) cycle
                                prod2 = MDM * prod1
                                ! MD coeff 3
                                do nn = 0, n(i) + n(j)
                                    H1 = (-1.0)**(ll + mm + nn)
                                    MDN = Dz1(nn + 1, n(i) + 1, n(j) + 1)
                                    if (abs(MDN)<1.0e-30) cycle ! check if MD coeff is 0
                                    prod3 = MDN * H1 * prod2
                                    ! MD coeff 4
                                    do llp = 0, l(k) + l(r)
                                        MDLp = Dx2(llp + 1, l(k) + 1, l(r) + 1)
                                        if (abs(MDLp)<1.0e-30) cycle ! check if MD coeff is 0
                                        prod4 = MDLp * prod3
                                        ! MD coeff 5
                                        do mmp = 0, m(k) + m(r)
                                            MDMp = Dy2(mmp + 1, m(k) + 1, m(r) + 1)
                                            if (abs(MDMp)<1.0e-30) cycle ! check if MD coeff is 0
                                            prod5 = MDMp * prod4
                                            ! MD coeff 6
                                            do nnp = 0, n(k) + n(r)
                                                MDNp = Dz2(nnp + 1, n(k) + 1, n(r) + 1)
                                                prod6 = MDNp * prod5
                                                ! cutoff after MD
                                                if (abs(prod6)>cutoff2) then

                                                    ! add the contribution to the total
                                                    F = F + prod6 * P0mat(:, ll + llp + 1, mm + mmp + 1, nn + nnp + 1)
                                                    !  print*,prod6, maxval(abs(P0mat(:,ll+llp+1,mm+mmp+1,nn+nnp+1))), ll+llp+1,mm+mmp+1,nn+nnp+1
                                                    !Int=Int+prodD.*f;
                                                end if
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        posr = posr + 1
                    end do
                    posk = posk + 1
                end do
                posj = posj + 1
            end do
            posi = posi + 1
        end do
        if (maxval(abs(F))>1E20) then
            print*, 'pzero is the problem'

        end if
    end subroutine


    subroutine tot_integral_k_ijkr(mu, l, m, n, group_start, group_count, hx, hy, hz, h, dx1, dy1, dz1, dx2, dy2, dz2, gi, gj, gk, gr, &
            zcontr, coeff1,coeff2,coeff3,coeff4, &
            cutoff1, cutoff2, f)

        implicit none


        !INTEGER, PARAMETER :: dp = selected_real_kind(2*precision(1.0_dp))
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        real(kind=dp), intent(in), dimension(:)::coeff1,coeff2,coeff3,coeff4
        integer(kind = selected_int_kind(8)), intent(in) :: gi, gj, gk, gr
        integer(kind = selected_int_kind(8)), dimension(:), intent(in) :: l, m, n, group_start, group_count
        real(kind = dp), intent(in) :: cutoff1, cutoff2, hx, hy, hz, h
        real(kind = dp), intent(in), dimension(:, :, :, :) :: zcontr
        real(kind = dp), intent(in), dimension(:), allocatable :: mu

        real(kind = dp), intent(in), dimension(:, :, :) :: dx1, dy1, dz1, dx2, dy2, dz2
        !real(kind=dp), pointer,dimension(:,:,:) :: dxx,dyy,dzz, dxx2,dyy2,dzz2

        integer(kind = selected_int_kind(8)) :: ka, posi, posj, posk, posr, ra, i, j, k, r
        integer(kind = selected_int_kind(8)) :: h1
        integer(kind = selected_int_kind(8)) :: ll, mm, nn, llp, mmp, nnp, llmax

        real(kind = dp) ::  prodd, ztot, mdn, mdl, mdm, mdlp, mdmp, mdnp, z11, z22
        INTEGER(kind = ikind), parameter :: dim = 20
        real(kind = dp) :: prod1, prod2, prod3, prod4, prod5, prod6
        REAL(kind = dp), DIMENSION(dim, dim) :: a, b, c
        REAL(kind = dp), DIMENSION(dim) :: h_saved
        REAL(kind = dp), DIMENSION(dim, dim, dim, dim) :: h_pre2
        REAL(kind = dp), DIMENSION(dim) :: BD



        !  real(kind=dp), dimension(:), allocatable :: pmu, h_0, h_1, h_r, muoh,zij, zij2, zkr, zkr2

        real(kind = dp), intent(out), dimension(size(mu)) :: f
        real(kind = dp), allocatable, dimension(:) :: l1vec, l2vec, suml1l2

        LLmax = l(group_start(gi)) + m(group_start(gi)) + n(group_start(gi)) + &
                l(group_start(gj)) + m(group_start(gj)) + n(group_start(gj)) + &
                l(group_start(gk)) + m(group_start(gk)) + n(group_start(gk)) + &
                l(group_start(gr)) + m(group_start(gr)) + n(group_start(gr))

        if (LLmax + 1 > dim) then
            print*, "only s, p, d, f, and g type GTOs are supported"
            stop
        end if

        a = 0.0_dp
        b = 0.0_dp
        c = 0.0_dp

        call rrdj0(LLmax, Hx, a)
        call rrdj0(LLmax, Hy, b)
        call rrdj0(LLmax, Hz, c)

        bd = 0.0_dp

        do k = 0, LLmax
            do j = 0, LLmax - k
                do i = 0, LLmax - k - j
                    call BesselDeriv(BD, k, j, i, a, b, c, LLmax)
                    h_pre2(:, k + 1, j + 1, i + 1) = BD
                end do
            end do
        end do

        ! posi=apos(i)

        h_saved = 0.0_dp

        ! loop through all possible ways to get total angular momentum lmax1

        posI = 1


        ! loop through all possible ways to get total angular momentum lmax1
        do i = group_start(gi), group_start(gi) + group_count(gi) - 1
            posj = 1
            do j = group_start(gj), group_start(gj) + group_count(gj) - 1
                posk = 1
                do k = group_start(gk), group_start(gk) + group_count(gk) - 1
                    posr = 1
                    do r = group_start(gr), group_start(gr) + group_count(gr) - 1

                        ztot = zcontr(posI, posJ, posK, posR)*coeff1(posI)*&
                                coeff2(posJ)*coeff3(posK)*coeff4(posR)


                        ! continue only if larger
                        if (abs(ztot) < cutoff1) then
                            posr = posr + 1
                            cycle
                        end if
                        do ll = 0, l(i) + l(j)
                            MDL = Dx1(ll + 1, l(i) + 1, l(j) + 1)
                            !  if (abs(MDL)<1.0e-30) cycle
                            prod1 = MDL * ztot
                            ! MD coeff 2
                            do mm = 0, m(i) + m(j)
                                MDM = Dy1(mm + 1, m(i) + 1, m(j) + 1)
                                !  if (abs(MDM)<1.0e-30) cycle
                                prod2 = MDM * prod1
                                ! MD coeff 3
                                do nn = 0, n(i) + n(j)
                                    H1 = (-1.0)**(ll + mm + nn)
                                    MDN = Dz1(nn + 1, n(i) + 1, n(j) + 1)
                                    !   if (abs(MDN)<1.0e-30) cycle ! check if MD coeff is 0
                                    prod3 = MDN * H1 * prod2
                                    ! MD coeff 4
                                    do llp = 0, l(k) + l(r)
                                        MDLp = Dx2(llp + 1, l(k) + 1, l(r) + 1)
                                        !      if (abs(MDLp)<1.0e-30) cycle ! check if MD coeff is 0
                                        prod4 = MDLp * prod3
                                        ! MD coeff 5
                                        do mmp = 0, m(k) + m(r)
                                            MDMp = Dy2(mmp + 1, m(k) + 1, m(r) + 1)
                                            !        if (abs(MDMp)<1.0e-30) cycle ! check if MD coeff is 0
                                            prod5 = MDMp * prod4
                                            ! MD coeff 6
                                            do nnp = 0, n(k) + n(r)
                                                MDNp = Dz2(nnp + 1, n(k) + 1, n(r) + 1)
                                                prod6 = MDNp * prod5
                                                ! cutoff after MD
                                                if (abs(prod6)>cutoff2) then
                                                    ! add the contribution to the total
                                                    h_saved = h_saved + h_pre2(:, ll + llp + 1, mm + mmp + 1, nn + nnp + 1) * prod6

                                                end if
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        posr = posr + 1
                    end do
                    posk = posk + 1
                end do
                posj = posj + 1
            end do
            posi = posi + 1
        end do

        CALL bessels0rr(F, llmax, mu, H, h_saved)

        !CALL BesselSum(F, mu, H, LLmax, h_saved)

    end subroutine


End Module TSj0groupsfast