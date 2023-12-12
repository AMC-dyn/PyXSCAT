Module TESfast

    use TSj0groupsfast
    use Bessels_j0
    use p0_cases
    use twordms
    use MD
    use onerdm
    use Zcontr
    use iso_fortran_env, only:  int8, int16, int32, int64
    implicit none
contains

    subroutine total_electron_scattering_fast(Zn,q, l, m, n, ngto, ng, nq, maxl, typec, state1, &
            state2, ncontr, group, gs, gf, gc, confs, ga, xx, yy, zz, coeffs, mmod, civecs, geom, &
            cutoffmd, cutoffz, cutoffcentre,contrvec,norbs,Iee_total,Iee_elastic,Ine,Inn, result2,result3)
        implicit none
        !Precision parameters
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

        !Entering variables from readers
        integer(kind = ikind), intent(in) :: ngto, ng, nq, maxl, typec, state1, state2, ncontr,norbs
        integer(kind = ikind), intent(in), dimension(:), allocatable :: l, m, n, group, gs, gf, gc,contrvec
        integer(kind = ikind), dimension(:, :), allocatable, intent(in) :: confs
        real(kind = dp), intent(in), dimension(:), allocatable :: ga, xx, yy, zz, q, coeffs
        real(kind = dp), intent(in), dimension(:, :), allocatable :: mmod, civecs, geom
        real(kind = dp), intent(in) :: cutoffmd, cutoffz, cutoffcentre
        real(kind=dp), intent(in), dimension(:), allocatable:: Zn

        !twordm variables
        integer(kind = ikind), dimension(:), allocatable :: m1, m2, m3, m4
        integer(kind = ikind), dimension(:, :), allocatable ::  ep3, ndiff2
        integer(kind=int64), dimension(:,:), allocatable:: mat
        integer(kind = ikind) :: nmat, i, j,k
        integer(kind=int64):: numberlines,maxnmo
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
        real(kind = dp), dimension(:, :, :, :), allocatable :: Zbig
        real(kind=dp),  dimension(:,:), allocatable :: onerdm_matrix
        real(kind=dp),dimension(ngto,ngto):: z
        real(kind=dp)::rij


        !Result out
        real(kind = dp), allocatable, dimension(:), intent(out) :: Result2,result3
        real(kind = dp), allocatable, dimension(:),intent(out) :: Iee_total,Iee_elastic, Ine, Inn


        allocate(Iee_total(nq),Iee_elastic(nq), Ine(nq), Inn(nq))
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
        call contraction(Zbig)


        call variables_total(px, py, pz, ddx, ddy, ddz, &
                e12, maxl, ngto, ng, group_start, group_count, group, ga, l, m, n, xx, yy, zz, &
                mmod,  q, nq)

        call tot_integration(ng, ncontr,px, py, pz, l, m, n, p0matrix, ddx, ddy, ddz, &
                group_start, group_count, group, &
                gs,gf,gc,contrvec,cutoffz, cutoffmd, cutoffcentre, q,coeffs, e12, Zbig, Iee_total)


        open(unit=15, file='bitwise.dat')
        numberlines=0
        do while(.true.)

            read (15, *, end=999) i
            numberlines=numberlines+1
        end do


        999 continue
        !numberlines=numberlines-1
        close(15)
        print*,numberlines
        file_in='bitwise.dat'
        file_out='es.dat'
        !call mcci_to_bit(file_in,file_out,numberlines)

        maxnmo=norbs
        print*,'maxumum number or orbitals, ', maxnmo
        call one_rdm_bit(file_out,onerdm_matrix, maxnmo,numberlines)

        call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
                mmod,onerdm_matrix,norbs,q,nq)

        print*,sum( (/ (onerdm_matrix(i,i), i=1, size(onerdm_matrix, 1)) /) )

        call elastic_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,Iee_elastic)
        call nuclei_electron_integration(Zn,geom,ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,Ine)

        print*,Ine(1)
        print*,'Nuclear', Zn
        Inn=0.0_dp
        do i=1,size(Zn)
            do j=1,size(Zn)
                do k=1,nq
                    Rij=sqrt(sum((geom(i,:)-geom(j,:))**2) )
                    Inn(k)=Inn(k)+Zn(i)*Zn(j)* sinc(q(k)*abs(Rij))

                end do
            enddo
        end do
        open(16,file='total_scat_from_el.dat')
        do k=1,nq
            write(16,*)q(k), Iee_total(k)+12
        end do
        close(16)
        open(16,file='mixed_scat_from_el.dat')
        do k=1,nq
            write(16,*)q(k),-2*Ine(k)
        end do

        close(16)

        open(16,file='nuclear_scat_from_el.dat')
        do k=1,nq
            write(16,*)q(k),Inn(k)
        end do
        close(16)
        result2=(sum( (/ (onerdm_matrix(i,i), i=1, size(onerdm_matrix, 1)) /) )+Iee_total-2*Ine+Inn)
        result3=Iee_elastic-2*Ine+Inn
    end subroutine

    subroutine variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl,ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
            mmod,onerdm,maxnmo,q,nq)


        use uniquemodule
        implicit none


        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=ikind), intent(in)::ngto,ng,  nq, maxl,maxnmo
        integer(kind=ikind), intent(in),dimension(:),allocatable :: l, m, n,group
        integer(kind=ikind), intent(in),dimension(:):: group_count,group_start
        real(kind=dp), intent(in),dimension(:),allocatable :: ga, xx, yy, zz, q
        real(kind=dp), intent(in),dimension(:,:),allocatable :: mmod
        real(kind=dp), intent(in),dimension(:,:)::onerdm
        real(kind=dp), intent(out), dimension(ng,ng):: px,py,pz
        real(kind=dp), intent(out), dimension(ngto,ngto) :: z
        real(kind=dp), intent(out), dimension(nq,ng,ng) :: e12
        real(kind=dp), intent(out), dimension(maxl*2+1,maxl+1,maxl+1,ng,ng) :: ddx,ddy,ddz


        integer(kind=ikind) :: max4l, max2l, i, j, k, ii,jj,ls, ms, ns,count,cc,cc2
        real(kind=dp) :: pi, gap
        real(kind=dp),  dimension(ngto,ngto,maxl*2+1,maxl+1,maxl+1):: dx,dy,dz
        real(kind=dp),  dimension(nq):: preexp
        real(kind=dp):: temp



        integer(kind=ikind) :: counter

        pi = acos(-1.0000)
        max4l=maxval(l)*4
        max2l=maxval(l)*2+1


        z=0.0_dp
        print*,maxval(mmod), maxval(onerdm)
        counter=0
        do  ii=1,ngto

            do jj=1,ngto

                temp=0.0_dp

                do cc=1,maxnmo
                    do cc2=1,maxnmo


                        temp=temp+onerdm(cc,cc2)*mmod(cc,ii)*mmod(cc2,jj)
                        temp=temp+onerdm(cc,cc2)*mmod(cc,jj)*mmod(cc2,ii)

                    enddo
                enddo
                Z(ii,jj)=temp/2.0_dp;

            enddo
        enddo
        print*,'Z(1,2)',Z(1,2)


        gap=0.0
        dx=0.0
        dy=0.0
        dz=0.0
        px=0.0
        py=0.0
        pz=0.0
        e12=0.0
        ddx=0.0
        ddy=0.0
        ddz=0.0





        call fill_md_table(dx,l,xx,ga)
        call fill_md_table(dy,m,yy,ga)
        call fill_md_table(dz,n,zz,ga)


        do jj = 1,Ng
            ! essentially taking the values for the first GTO in the group.
            ! All gtos in the group have the same prefactors as the prefactors
            ! do not't depend on l, m and n
            j = group_start(jj)
            do ii = 1, Ng
                i = group_start(ii)
                gaP=ga(i)+ga(j)
                Px(ii,jj)=(ga(i)*xx(i) + ga(j)*xx(j))/gaP
                Py(ii,jj)=(ga(i)*yy(i) + ga(j)*yy(j))/gaP
                Pz(ii,jj)=(ga(i)*zz(i) + ga(j)*zz(j))/gaP
                E12(:,ii,jj) = (pi/gaP)**1.5 * exp(-q*q*0.25/gaP) &
                        * exp(-ga(i)*ga(j)/gaP*((xx(i)-xx(j))**2. + (yy(i)-yy(j))**2. + (zz(i)-zz(j))**2.))
            end do
        end do



        do jj = 1, Ngto
            j = group(jj)
            do ii = 1, Ngto
                i = group(ii)
                do ls = 1, l(ii)+l(jj)+1
                    Ddx(ls, l(ii)+1,l(jj)+1,j,i)=Dx(ii,jj,ls, l(ii)+1,l(jj)+1)
                end do

                do ms = 1, m(ii)+m(jj)+1
                    Ddy(ms, m(ii)+1,m(jj)+1,j,i)=Dy(ii,jj,ms, m(ii)+1,m(jj)+1)
                end do

                do ns = 1, n(ii)+n(jj)+1
                    Ddz(ns, n(ii)+1,n(jj)+1,j,i)=Dz(ii,jj,ns, n(ii)+1,n(jj)+1)
                end do
            end do
        end do

        print*,maxval(e12)
    end subroutine variables_elastic
    subroutine nuclei_electron_integration(Zn,geom,ncap,px,py,pz,l,m,n,p0matrix,dx,dy,dz,z,group_start,group_count,group, &
            cutoffz,cutoffmd, cutoffcentre,q,e12,tsi)
        use omp_lib
        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(kind=ikind), INTENT(IN) :: ncap
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:),allocatable :: l,m,n,group
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:)::  group_start,group_count

        REAL(kind=dp), intent(in), dimension(:,:,:,:,:),target :: dx,dy,dz
        real(kind=dp),dimension(:,:,:), pointer :: dx1,dy1,dz1
        real(kind=dp),dimension(:,:,:,:), intent(in) :: p0matrix
        REAL(kind=dp), intent(in), dimension(:,:,:) ::e12
        real(kind=dp), intent(in), dimension(:,:)::z,geom

        REAL(kind=dp), intent(in), dimension(:,:) :: px,py,pz
        REAL(kind=dp), intent(in), dimension(:),allocatable :: q,Zn
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

        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:)         :: gs,gc
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:),allocatable::l,m,n
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
        real(kind=dp), intent(in), dimension(:),allocatable   ::  mu

        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1
        !real(kind=dp), pointer,dimension(:,:,:) :: dxx,dyy,dzz, dxx2,dyy2,dzz2



        integer(kind=selected_int_kind(8)) :: ka, posi, posj, ra,i,j,k,r
        integer(kind=selected_int_kind(8)) ::  h1
        integer(kind=selected_int_kind(8)) ::  ll, mm, nn, llp, mmp, nnp,llmax

        real(kind=dp)  ::coeff,prodd,ztot,mdn, mdl, mdm, mdlp, mdmp,mdnp,z1
        INTEGER(kind=ikind), parameter :: dim = 20
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
            print*, "only s, p, d, f and g type GTOs are supported"
            stop
        end if



        !call Hermite_like_coeffs(a, LLmax, Hx)
        !call Hermite_like_coeffs(b, LLmax, Hy)
        !call Hermite_like_coeffs(c, LLmax, Hz)


        a = 0.0_dp
        b = 0.0_dp
        c = 0.0_dp

        call rrdj0(LLmax, Hx, a)
        call rrdj0(LLmax, Hy, b)
        call rrdj0(LLmax, Hz, c)



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




        !CALL BesselSum(F, mu, H, LLmax, h_saved)
        CALL bessels0rr(F, llmax, mu, H, h_saved)

    end subroutine
    function sinc (a)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        real(kind=dp):: sinc, a
        if (abs(a) < 1.0d-10) then
            sinc = 1
        else
            sinc = sin(a) / (a)
        end if
    end function



    subroutine elastic_integration(ncap,px,py,pz,l,m,n,p0matrix,dx,dy,dz,z,group_start,group_count,group, &
            cutoffz,cutoffmd, cutoffcentre,q,e12,tsi)


        use omp_lib
        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(kind=ikind), INTENT(IN) :: ncap
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:), allocatable :: l,m,n,group
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:) :: group_start,group_count

        REAL(kind=dp), intent(in), dimension(:,:,:,:,:),target :: dx,dy,dz
        real(kind=dp),dimension(:,:,:), pointer :: dx1,dy1,dz1,dx2,dy2,dz2
        real(kind=dp),dimension(:,:,:,:), intent(in) :: p0matrix
        REAL(kind=dp), intent(in), dimension(:,:,:) ::e12
        real(kind=dp), intent(in), dimension(:,:)::z

        REAL(kind=dp), intent(in), dimension(:,:) :: px,py,pz
        REAL(kind=dp), intent(in), dimension(:),allocatable :: q
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






        nq= size(q)

        print*,OMP_get_num_threads()
        call omp_set_num_threads(16)
        print*,OMP_get_num_threads()
        !First big loop
        print*,'posits allocated'
        tsi=0.0_dp
        !$OMP PARALLEL do private(posI,posK,posJ,posR,spi,spj,spk,spr,za,zb), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,j,k,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz, posits,cutoffmd,group_count,group_start,z) REDUCTION(+:tsi), &
        !$OMP& schedule(dynamic)

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

        print*,'first loop finished', tsi(1)


        !$OMP PARALLEL do private(posI,posJ,posR,spi,spj,spk,spr,za,zb), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,j,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz,posits, cutoffmd,group_count,group_start,z) REDUCTION(+:tsi), &
        !$OMP& schedule(dynamic)
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

        !$OMP PARALLEL do private(posI,posK,posR,spi,spj,spk,spr,za,zb), &
        !$OMP& private(f,ii,jj,h,hx,hy,hz,i,k,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz,posits, cutoffmd,group_count,group_start,z) REDUCTION(+:tsi), &
        !$OMP& schedule(dynamic)
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
        !$OMP& shared( cutoffz, cutoffmd,posits,group_count,group_start) REDUCTION(+:tsi), &
        !$OMP& schedule(dynamic)
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
        !$OMP& shared( cutoffz, cutoffmd,posits,group_count,group_start) REDUCTION(+:tsi), &
        !$OMP& schedule(dynamic)
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


    end subroutine elastic_integration
    subroutine integral_k_ijkr(mu,l,m,n,group_start,group_count,hx,hy,hz,h,dx1,dy1,dz1,dx2,dy2,dz2, gi,gj,gk, gr,&
            za,zb,&
            cutoff1, cutoff2,f)

        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        !INTEGER, PARAMETER :: dp = selected_real_kind(2*precision(1.0_dp))
        integer(kind=selected_int_kind(8)), intent(in)  :: gi,gj,gk,gr
        integer(kind=selected_int_kind(8)), dimension(:), intent(in),allocatable :: l,m,n
        integer(kind=selected_int_kind(8)), dimension(:), intent(in) :: group_start,group_count
        real(kind=dp), intent(in)              :: cutoff1,cutoff2, hx, hy, hz, h
        real(kind=dp), intent(in), dimension(:,:) :: za,zb
        real(kind=dp), intent(in), dimension(:),allocatable   ::  mu

        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1,dx2,dy2,dz2
        !real(kind=dp), pointer,dimension(:,:,:) :: dxx,dyy,dzz, dxx2,dyy2,dzz2



        integer(kind=selected_int_kind(8)) :: ka, posi, posj, posk, posr, ra,i,j,k,r
        integer(kind=selected_int_kind(8)) ::  h1
        integer(kind=selected_int_kind(8)) ::  ll, mm, nn, llp, mmp, nnp,llmax

        real(kind=dp)  ::coeff,prodd,ztot,mdn, mdl, mdm, mdlp, mdmp,mdnp,z1
        INTEGER(kind=ikind), parameter :: dim = 20
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
            print*, "only s, p, d, f, and g type GTOs are supported"
            stop
        end if

        a = 0.0_dp
        b = 0.0_dp
        c = 0.0_dp

        call rrdj0(LLmax, Hx, a)
        call rrdj0(LLmax, Hy, b)
        call rrdj0(LLmax, Hz, c)






        bd=0.0_dp
        !h_pre2=0.0_dp
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

        if (sum(h_saved)==0.0_dp) then
            F=0.0_dp
            return
        end if




        call bessels0rr(F,LLmax,mu,H,h_saved)

        !        do i=1,LLmax+1
        !            F=F+h_saved(i)*bess(i,:)
        !
        !        end do
        !CALL BesselSum(F, mu, H, LLmax, h_saved)



    end subroutine
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

        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:)         :: gs,gc
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:),allocatable::l,m,n
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
                            !   if (abs(MDL)<1.0e-30) cycle
                            prod1 = MDL * ztot
                            ! MD coeff 2
                            do mm = 0, m(i)+m(j)
                                MDM = Dy1(mm+1,m(i)+1,m(j)+1)
                                !   if (abs(MDM)<1.0e-30) cycle
                                prod2 = MDM * prod1
                                ! MD coeff 3
                                do nn =0, n(i)+n(j)
                                    H1=(-1.0)**(ll+mm+nn)
                                    MDN=Dz1(nn+1,n(i)+1,n(j)+1)
                                    !       if (abs(MDN)<1.0e-30) cycle ! check if MD coeff is 0
                                    prod3 = MDN * H1  * prod2
                                    ! MD coeff 4
                                    do llp = 0, l(k)+l(r)
                                        MDLp=Dx2(llp+1,l(k)+1,l(r)+1)
                                        !         if (abs(MDLp)<1.0e-30) cycle ! check if MD coeff is 0
                                        prod4 = MDLp * prod3
                                        ! MD coeff 5
                                        do mmp = 0, m(k)+m(r)
                                            MDMp = Dy2(mmp+1,m(k)+1,m(r)+1)
                                            !     if (abs(MDMp)<1.0e-30) cycle ! check if MD coeff is 0
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



end module