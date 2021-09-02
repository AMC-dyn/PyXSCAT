!-----------------------------------------------------------------------
! 2RDM Fortran Modules for PyXSCAT code
! Andres Moreno Carrascosa and Mats Simmermacher, (2021)
!-----------------------------------------------------------------------

MODULE types

        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

END MODULE types



MODULE uniquemodule

    implicit none 

    contains


    SUBROUTINE unique_integer(matrix,rmat,iuni,irec)
    
        use types
        
        integer(kind=ikind), dimension(:,:), intent(in)                  :: matrix
        integer(kind=ikind), dimension(:,:), allocatable                 :: sprmat
        integer(kind=ikind), dimension(:,:), intent(out), allocatable    :: rmat
        integer(kind=ikind), dimension(:), allocatable                   :: same, iunifull
        integer(kind=ikind), dimension(:), intent(out), allocatable      :: iuni
        integer(kind=ikind), dimension(:), intent(out), allocatable      :: irec
        integer(kind=ikind)                                              :: i, cnt

        allocate(sprmat(size(matrix(:,1)),size(matrix(1,:))))
        allocate(same(size(matrix(:,1))), iunifull(size(matrix(:,1))))
        allocate(irec(size(matrix(:,1))))
        
        irec = 0
        cnt = 0
            
        do i = 1, size(matrix(:,1))   
            if (irec(i) == 0) then
                cnt = cnt + 1
                iunifull(cnt) = i
                sprmat = spread(matrix(i,:), 1, size(matrix(:,i)))
                same = sum(abs(matrix - sprmat), dim=2)
                same = (-1)*same + 1
                same = (abs(same) + same)/2
                irec = irec + same*cnt
            endif
        enddo

        allocate(iuni(cnt))
        allocate(rmat(cnt,size(matrix(1,:))))

        iuni = iunifull(1:cnt)
        rmat = matrix(iuni,:)

    END SUBROUTINE


    SUBROUTINE unique_real(matrix,rmat,iuni,irec)
    
        use types
        
        real(kind=dp), dimension(:,:), intent(in)                      :: matrix
        real(kind=dp), dimension(:,:), allocatable                     :: sprmat
        real(kind=dp), dimension(:,:), intent(out), allocatable        :: rmat
        integer(kind=ikind), dimension(:), allocatable                 :: same, iunifull
        integer(kind=ikind), dimension(:), intent(out), allocatable    :: iuni
        integer(kind=ikind), dimension(:), intent(out), allocatable    :: irec
        integer(kind=ikind)                                            :: i, cnt

        allocate(sprmat(size(matrix(:,1)),size(matrix(1,:))))
        allocate(same(size(matrix(:,1))), iunifull(size(matrix(:,1))))
        allocate(irec(size(matrix(:,1))))
        
        irec = 0
        cnt = 0
            
        do i = 1, size(matrix(:,1))   
            if (irec(i) == 0) then
                cnt = cnt + 1
                iunifull(cnt) = i
                sprmat = spread(matrix(i,:), 1, size(matrix(:,i)))
                same = ceiling(sum(abs(matrix - sprmat), dim=2))
                same = (-1)*same + 1
                same = (abs(same) + same)/2
                irec = irec + same*cnt
            endif
        enddo

        allocate(iuni(cnt))
        allocate(rmat(cnt,size(matrix(1,:))))

        iuni = iunifull(1:cnt)
        rmat = matrix(iuni,:)

    END SUBROUTINE


    SUBROUTINE unique_total(matrix,total,rmat,rtot)
    
        use types
        
        integer(kind=ikind), dimension(:,:), intent(in)                  :: matrix
        real(kind=dp), dimension(:), intent(in)                          :: total
        integer(kind=ikind), dimension(:,:), allocatable                 :: sprmat
        integer(kind=ikind), dimension(:,:), intent(out), allocatable    :: rmat
        real(kind=dp), dimension(:), intent(out), allocatable            :: rtot
        integer(kind=ikind), dimension(:), allocatable                   :: same, iunifull
        real(kind=dp), dimension(:), allocatable                         :: rtotfull
        integer(kind=ikind), dimension(:), allocatable                   :: iuni
        integer(kind=ikind), dimension(:), allocatable                   :: irec
        integer(kind=ikind)                                              :: i, cnt, dim1, dim2


        dim1 = size(matrix(:,1))
        dim2 = size(matrix(1,:))
        print*,dim1,dim2
        allocate(sprmat(dim1,dim2))
        print*, 'first allocation'
        allocate(same(dim1), iunifull(dim1))
        print*,'second allocation'
        allocate(irec(dim1))
        print*,'third allocation'
        allocate(rtotfull(dim1))
        print*,'fourth allocation'
        print*,dim1
        irec = 0
        cnt = 0
            
        do i = 1, dim1
            if (irec(i) == 0) then
                cnt = cnt + 1
                iunifull(cnt) = i
                sprmat = spread(matrix(i,:), 1, size(matrix(:,i)))
                same = sum(abs(matrix - sprmat), dim=2)
                same = (-1)*same + 1
                same = (abs(same) + same)/2
                rtotfull(cnt) = sum(same*total)
                irec = irec + same*cnt

            endif
        enddo

        allocate(iuni(cnt))
        allocate(rmat(cnt,size(matrix(1,:))))
        allocate(rtot(cnt))
        
        iuni = iunifull(1:cnt)
        rmat = matrix(iuni,:)
        rtot = rtotfull(1:cnt)

        print*, 'reduced finished'
    END SUBROUTINE


END MODULE uniquemodule



MODULE membermodule

    implicit none 

    contains


    SUBROUTINE ismember(col,matrix,memval,colnum)

        use types
    
        integer(kind=ikind), intent(in), dimension(4)        :: col
        integer(kind=ikind), intent(in), dimension(:,:)      :: matrix
        integer(kind=ikind), dimension(size(matrix(:,1)))    :: subcol
        logical, intent(out)                                 :: memval
        integer(kind=ikind), intent(out)                     :: colnum

        subcol = sum(abs(matrix - spread(col, 1, size(matrix(:,1)))), dim=2)
        colnum = minloc(subcol, dim=1)

        if (subcol(colnum) == 0) then
           memval = .true.
        else
           colnum = -1
           memval = .false.
        endif

    END SUBROUTINE


END MODULE membermodule



MODULE twordmreader

    implicit none

    contains


    SUBROUTINE BesselDeriv(BD, LL, MM,NN,a,b,c,LLmax)
        use types
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


    SUBROUTINE reduce_density(matst,totalst,m1,m2,m3,m4,total2)

        use types
        use membermodule
        use uniquemodule

        real(kind=dp), intent(in), dimension(:)                    :: totalst
        integer(kind=ikind), intent(in), dimension(:,:) :: matst
        integer(kind=ikind), dimension(:), intent(out), allocatable      :: m1, m2, m3, m4
        real(kind=dp), dimension(:), intent(out),allocatable             :: total2
        integer(kind=ikind)                                              :: i, j, k, l, cc
        integer(kind=ikind)                                              :: n, sdr
        integer(kind=ikind)                                              :: cnt, num
        logical                                                          :: memb

        integer(kind=ikind), dimension(:,:), allocatable                 :: newmat
        real(kind=dp), dimension(:), allocatable      :: stotal,newtotal
        integer(kind=ikind), dimension(:,:), allocatable    :: smat
        integer(kind=ikind)                                              :: mo1, mo2, mo3, mo4
        integer(kind=ikind), dimension(4)                                :: b

        integer(kind=ikind), dimension(:,:), allocatable                 :: mat2

        integer(kind=ikind), dimension(:), allocatable :: num1




        call unique_total(matst, totalst, mat2, total2 )
        
        sdr = size(mat2(:,1))

        allocate(m1(sdr), m2(sdr), m3(sdr), m4(sdr))


        m1 = mat2(:,1)
        m2 = mat2(:,2)
        m3 = mat2(:,3)
        m4 = mat2(:,4)
        


   
    END SUBROUTINE reduce_density


     SUBROUTINE integral_ijkr_pzero(nq,l,m,n,gs,gc,p0mat,dx1,dy1,dz1,dx2,dy2,dz2,gi,gj,gk,gr,zcontr,zcontr2, &
             cutoff1,cutoff2,f)

        use types

        implicit none

        ! definition of input

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

        use types

        implicit none


        !INTEGER, PARAMETER :: dp = selected_real_kind(2*precision(1.0_dp))
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


        SUBROUTINE BesselSum(h_sum, mu, H, LLmax, h_saved)
            use types
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
            use types
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


    subroutine maxcoincidence(confs, ep2,ndiff)
        use types
        implicit none
        integer(kind=ikind), intent(in), dimension(:,:) :: confs
        integer(kind=ikind), intent(out), dimension(:,:), allocatable :: ep2, ndiff
        integer(kind=ikind), dimension(size(confs(:,1)), size(confs(1,:))):: matdum
        integer(kind=ikind), dimension(:,:), allocatable :: mat1
        integer(kind=ikind) :: i,j, c1,c2, count


        allocate(ep2(size(confs(:,1)),size(confs(:,1))))
        ep2=1
        count=0
        matdum=0
        print*,'holaa'
            do i=1,size(confs(:,1))
                count=0
                do j=1,size(confs(1,:))


                    if (confs(i,j)/=0) then
                        count=count+1
                        matdum(i,count)=j
                        end if
                end do
            end do
           print*, 'matdum constructed'
            allocate(mat1(size(confs(:,1)),count), ndiff(size(confs(:,1)),size(confs(:,1))))
            ndiff=0
            mat1=matdum(:,1:count)
            print*, 'mat1 constructed'
            do c1=1, size(confs(:,1))
                do c2=c1+1,size(confs(:,1))
                    do i=1,size(mat1(1,:))
                        if  (mat1(c1,i) /= mat1(c2,i)) then


                            do j=1,size(mat1(1,:))
                                if (mat1(c1,i) /= mat1(c2,j)) then
                                    ep2(c1,c2)=-ep2(c1,c2)
                                    ep2(c2,c1)=ep2(c1,c2)

                                end if
                            end do
                        end if
                    end do

                    do j=1,size(confs(1,:))
                        if (confs(c1,j)/=confs(c2,j)) then
                           ndiff(c1,c2)= ndiff(c1,c2) +1
                           ndiff(c2,c1)= ndiff(c2,c1) +1
                        end if
                    end do

                end do
            end do

    print*,ndiff(1,2)
    end subroutine maxcoincidence

    subroutine createtwordm(confs,civs,ndiff,ep2,mat,total)

    implicit none

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

    integer(kind=ikind), intent(in),dimension(:,:) :: confs
    integer(kind=ikind), intent(in), dimension(:,:) :: ep2, ndiff
    real(kind=dp), intent(in), dimension(:,:) :: civs
    real(kind=dp), intent(out), dimension(:), allocatable :: total
    integer(kind=ikind), intent(out), dimension(:,:), allocatable :: mat

    real(kind=dp), dimension(:,:,:,:), allocatable :: twordm
    integer(kind=ikind) :: ep,nc1,lconfs,norbs,nc2,sorb,rorb,qorb,porb,p,q,r,s,c1,c2,count1,count

    integer(kind=ikind) :: i,i1,i2,n,count2,eg, ndiff1

    integer(kind=ikind), dimension(:), allocatable :: mat1,mat2
    logical(4) :: sdef, rdef, pdef, qdef
    logical(4), dimension(:), allocatable :: logicaltwordms
    real(kind=dp) :: cutoff
    integer(kind=ikind), dimension(:), allocatable :: diffs1,diffs2
    integer(kind=ikind) :: spins,spinr,spinq, spinp
    integer(kind=ikind), dimension(:), allocatable :: spin1,spin2
    integer(kind=ikind),  dimension(:,:), allocatable :: matdum
    real(kind=dp), dimension(:), allocatable :: totaldum

    lconfs=size(confs(1,:))
    norbs=lconfs/2
    allocate(twordm(norbs,norbs,norbs,norbs))
    twordm=0.0_dp
    nc1=1
    nc2=1
    ep=1

    do c1=1,size(confs(:,1))
        do c2=1,size(confs(:,1))

            ndiff1=ndiff(c1,c2)
            ep=ep2(c1,c2)

            if (ndiff(c1,c2)/=0 .and. ndiff1 <=4) then

                sdef = .False.
                rdef = .False.
                pdef = .False.
                qdef = .False.
                if (allocated(diffs1)) deallocate(diffs1)
                if (allocated(diffs2)) deallocate(diffs2)
                if (allocated(spin1)) deallocate(spin1)
                if (allocated(spin2)) deallocate(spin2)
                allocate(diffs1(ndiff(c1,c2)/2), diffs2(ndiff(c1,c2)/2), spin1(ndiff(c1,c2)/2))
                allocate(spin2(ndiff(c1,c2)/2))

                count1=1
                count2=1

                do n=1,size(confs(c1,:))
                    if (confs(c1,n) /= confs(c2,n)) then
                        if (confs(c1,n) /= 0) THEN
                            diffs1(count1)=nint((n) / 2.0_dp + 0.1)
                            spin1(count1)=confs(c1,n)
                            count1=count1+1


                        elseif (confs(c2,n) /= 0) THEN
                            diffs2(count2)=nint((n) / 2.0_dp + 0.1)
                            spin2(count2)=confs(c2,n)
                            count2=count2+1
                        end if
                    end if
                enddo

                if (ndiff(c1,c2) == 4) then

                    sorb = diffs2(1)
                    qorb = diffs2(2)
                    porb = diffs1(2)
                    rorb = diffs1(1)

                    eg = 1.00

                    if (spin2(1) == spin1(1) .and. spin2(2) == spin1(2)) THEN
                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
                                civs(c1,1) * civs(c2,1) * ep * eg

                    end if

                    sorb = diffs2(1)
                    qorb = diffs2(2)
                    rorb = diffs1(2)
                    porb = diffs1(1)
                    eg = -1.00

                    if (spin2(1) == spin1(2) .and. spin2(2) == spin1(1)) then
                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ civs(c1,1) * civs(c2,1) * ep * eg
                    endif
                    qorb = diffs2(1)
                    sorb = diffs2(2)
                    porb = diffs1(2)
                    rorb = diffs1(1)
                    eg = -1.00

                    if (spin1(1) == spin2(2) .and. spin2(1) == spin1(2) ) THEN
                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
                                civs(c1,1) * civs(c2,1) * ep * eg


                    end if

                    qorb = diffs2(1)
                    sorb = diffs2(2)
                    rorb = diffs1(2)
                    porb = diffs1(1)
                    eg = 1.00

                    if (spin1(2) == spin2(2) .and. spin2(1) == spin1(1)) then
                        twordm(porb , rorb , sorb , qorb )= twordm(porb , rorb , sorb , qorb ) &
                                + civs(c1,1) * civs(c2,1) * ep * eg

                    end if


                elseif (ndiff(c1,c2) == 2) THEN

                    qorb = diffs2(1)
                    porb = diffs1(1)
                    eg = 1.00

                    do i=1,size(confs(c2,:))

                        if (confs(c2,i) /=0) then
                            sdef = .True.

                            sorb =nint(i / 2.0_dp + 0.1)

                            spins = confs(c2,i)
                        end if

                        if (confs(c1,i) /= 0) then

                            rdef = .True.
                            rorb = nint((i) / 2.0_dp + 0.1)
                            spinr = confs(c1,i)
                        end if

                        if (sdef .and. rdef .and. spins == spinr .and. spin2(1) == spin1(1)) then
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.

                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) &
                                    + civs(c1,1) * civs(c2,1) * ep * eg



                        else
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.

                        end if
                    enddo
                    qorb = diffs2(1)
                    rorb = diffs1(1)
                    eg = -1.00

                    do  i= 1, size(confs(c2,:))

                        if (confs(c2,i) /= 0) then
                            sdef = .True.
                            sorb = nint((i) / 2.0_dp + 0.1)
                            spins = confs(c2,i)
                        endif
                        if (confs(c1,i) /= 0) then
                            pdef = .True.
                            porb = nint((i) / 2.0_dp + 0.1)
                            spinp = confs(c1,i)
                        end if

                        if (sdef .and. pdef .and. spin1(1) == spins .and. spin2(1) == spinp) then
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.
                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) + &
                                    civs(c1,1) * civs(c2,1) * ep * eg


                        else
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.

                        end if
                    end do



                    sorb = diffs2(1)
                    porb = diffs1(1)
                    eg = -1.00

                    do i=1,size(confs(c2,:))

                        if (confs(c2,i) /= 0) then
                            qdef = .True.
                            qorb = nint((i) / 2.0_dp + 0.1)
                            spinq = confs(c2,i)
                        endif
                        if (confs(c1,i) /= 0) then
                            rdef = .True.
                            rorb = nint((i) / 2.0_dp + 0.1)
                            spinr = confs(c1,i)
                        endif
                        if (rdef .and. qdef .and. spins == spinr .and. spin1(1) == spin2(1)) then
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.
                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) +&
                                    civs(c1,1) * civs(c2,1) * ep * eg


                        else
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.

                        end if
                    enddo

                    sorb = diffs2(1)
                    rorb = diffs1(1)
                    eg = 1.00

                    do i=1,size(confs(c2,:))
                        if (confs(c2,i) /= 0) then
                            qdef = .True.
                            qorb = nint((i) / 2.0_dp + 0.1)
                            spinq = confs(c2,i)

                        end if

                        if (confs(c1,i) /= 0) then
                            pdef = .True.
                            porb = nint((i) / 2.0_dp + 0.1)
                            spinp = confs(c1,i)

                        end if

                        if (qdef .and. pdef .and. spinq == spinp .and. spin2(1) == spin1(1)) then
                            qdef = .False.
                            pdef = .False.

                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) &
                                    + civs(c1,1) * civs(c2,1) * ep * eg

                        else

                            pdef = .False.
                            qdef = .False.

                        end if
                    enddo
                endif

            elseif (ndiff(c1,c2) == 0) then

                ep = 1

                do i1=1,size(confs(c1,:))

                    do i2=1,size(confs(c2,:))

                        if (i1/=i2) then

                            if (confs(c1,i1) /= 0 .and. confs(c1,i2) /= 0) then

                                sorb = nint(i1 / 2.0_dp + 0.1)
                                qorb = nint(i2 / 2.0_dp + 0.1)

                                if (confs(c1,i1) == confs(c1,i2)) then
                                    porb = sorb
                                    rorb = qorb

                                    eg = -1.00

                                    twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) &
                                            + civs(c1,1) * civs(c2,1) * ep * eg

                                end if
                                porb = qorb
                                rorb = sorb
                                eg = 1.00

                                twordm(porb , rorb , sorb , qorb )  = twordm(porb , rorb , sorb , qorb ) &
                                        + civs(c1,1) * civs(c2,1) * ep * eg
                            endif
                        end if

                    enddo
                end do
            end if


        end do

    end do



    cutoff = 1E-09
    count=0
    count2=1

    allocate(logicaltwordms(norbs**4))
    allocate(totaldum(norbs**4), matdum(norbs**4,4))
    logicaltwordms(:)=.False.
    do p=1,norbs
        do q=1,norbs
            do r=1,norbs
                do s=1,norbs
                    totaldum(count2)=twordm(p,q,r,s)
                    matdum(count2,:)=(/p,s,q,r/)
                    if (abs(twordm(p,q,r,s))>=cutoff) then
                        count=count+1
                        logicaltwordms(count2)=.True.

                    end if
                    count2=count2+1
                end do
            end do
        end do
    end do


    allocate(mat(count, 4), total(count))
    count=1
    do i=1,count2
        if (logicaltwordms(i)) then
            mat(count,:)=matdum(i,:)
            total(count)=totaldum(i)
            count=count+1
        end if
    end do

        print*, 'twordm calculated'



    end subroutine createtwordm

subroutine variables_total(px,py,pz,ddx,ddy,ddz,z11,z22,e12,maxl,ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q,nq)

    use types
    use uniquemodule
        implicit none



        integer(kind=ikind), intent(in)::ngto,ng, nmat, nq, maxl
        integer(kind=ikind), intent(in),dimension(:) :: l, m, n, m1, m2, m3, m4,group_start,group,group_count


        real(kind=dp), intent(in),dimension(:) :: ga, xx, yy, zz, total,q
        real(kind=dp), dimension(ng):: ga2, xx2, yy2, zz2
        real(kind=dp), intent(in),dimension(:,:) :: mmod
        real(kind=dp), intent(out), dimension(ng,ng):: px,py,pz
        real(kind=dp), dimension(nmat,ngto,ngto) :: z1
        real(kind=dp), dimension(:,:),allocatable :: za
        real(kind=dp), dimension(:,:),allocatable :: zb
        real(kind=dp), dimension(nmat,ngto,ngto) ::  z2
        real(kind=dp), dimension(ngto,ngto) :: cmat
        !real(kind=dp), intent(out), dimension(ngto,ngto,ngto,ngto) :: zcontr

        real(kind=dp), intent(out), dimension(nq,ngto,ngto) :: e12

        real(kind=dp), intent(out), dimension(maxl*2+1,maxl+1,maxl+1,ng,ng) :: ddx,ddy,ddz


        integer(kind=ikind) :: max4l, max2l, i, j, k, N1, N2, nn1, nn2,ii,jj,ls, ms, ns,count
        integer(kind=ikind), dimension(ngto) :: ll
        integer(kind=ikind), dimension(nmat) :: indexing1,indexing2
        real(kind=dp), dimension(:,:,:),intent(out), allocatable :: z11, z22

        real(kind=dp) :: pi, gap,time1,time2,time3,time4

        integer(kind=ikind), dimension(:), allocatable :: iduplicates,jduplicates,m11,m22,m33,m44
        real(kind=dp),  dimension(ngto,ngto,maxl*2+1,maxl+1,maxl+1):: dx,dy,dz
        real(kind=dp),  dimension(nq):: preexp
        real(kind=dp),  dimension(:), allocatable:: temp1,temp2
        integer(kind=ikind),dimension(nmat,2) :: vec1, vec2
        integer(kind=ikind),dimension(nmat,4) :: mat1
        integer(kind=ikind),dimension(:,:),allocatable :: matfin
        real(kind=dp),dimension(:), allocatable :: totalfin

        pi = acos(-1.0000)
        max4l=maxval(l)*4
        max2l=maxval(l)*2+1


       ! N1=size(ipos)
        do i=1,nmat

            vec1(i,:)=(/m1(i), m2(i)/)
            vec2(i,:)=(/m3(i), m4(i)/)
            call bubble_sort(vec1(i,:))
            call bubble_sort(vec2(i,:))
            mat1(i,1:2)=vec1(i,:)
            mat1(i,3:4)=vec2(i,:)

        end do

        call unique_total(mat1,total,matfin,totalfin)
        allocate(m11(size(totalfin)),m22(size(totalfin)),m33(size(totalfin)),m44(size(totalfin)))
        allocate(z11(size(totalfin),ngto,ngto), z22(size(totalfin),ngto,ngto), &
               temp1(size(totalfin)),temp2(size(totalfin)))
        m11 = matfin(:,1)
        m22 = matfin(:,2)
        m33 = matfin(:,3)
        m44 = matfin(:,4)

        write(*,*)'shape of unique mat',shape(totalfin)
      !  print*,shape(z11), nmat
        z11=0.0_dp
        z22=0.0_dp
        print*,z11(1,1,1),m11(1)


        temp1=0.0
        temp2=0.0
        print*,ngto,shape(z11),shape(z22),shape(mmod)


         do  ii=1,ngto

            do jj=1,ngto


                        temp1 = totalfin * (mmod(m11, ii) * mmod(m22, jj) + mmod(m11, jj) * mmod(m22, ii))
                        temp2 = mmod(m33, ii) * mmod(m44, jj) + mmod(m33, jj) * mmod(m44, ii)

                        z11(:,ii, jj) =  temp1
                        z22(:, ii, jj) = temp2



            enddo
        enddo



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


        call cpu_time(time3)



        call fill_md_table(dx,l,xx,ga)
        call fill_md_table(dy,m,yy,ga)
        call fill_md_table(dz,n,zz,ga)


       ! allocate( ga2(size(apos)), xx2(size(apos)), yy2(size(apos)), zz2(size(apos)) )

!        ll = l + m + n
!        ga2 = ga(apos)
!        xx2 = xx(apos)
!        yy2 = yy(apos)
!        zz2 = zz(apos)
!        ll2 = ll(apos)

     !   N2=size(apos)

       ! allocate(ddx(n2,n2,max2l,maxl,maxl),ddy(n2,n2,max2l,maxl,maxl),ddz(n2,n2,max2l,maxl,maxl))
        !allocate(preexp(nq))

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
        print*,sum(Pz)


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



    print*, 'leaving variables'
    end subroutine variables_total


    SUBROUTINE fill_md_table(D, l, x, ga)
        use types
        implicit none

        ! The MD table to be populated D(Ngto, Ngto, 2maxl+1,maxl+1,maxl+1)
        REAL(kind=dp), INTENT(OUT), DIMENSION(:,:,:,:,:)    :: D
        ! vectors that contains the x coordinates for all GTOs, and all gammas
        REAl(kind=dp), INTENT(IN), DIMENSION(:)             :: ga, x
        ! vector that contains all the angular momenta (lx) for the GTOs
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:)       :: l


        ! loop and temp vars
        INTEGER(kind=ikind) :: i, j, ii, jj, N, maxl, l1, l2
        REAL(kind=dp)       :: gaP, Px, PA, PB, a

        ! number of GTOs
        N = size(l)
        ! maximum angular momentum
        maxl = maxval(l)


        ! the offdiagonal elements
        do i = 1,N
            do j = i+1, N
                gaP=ga(i)+ga(j)
                Px=(ga(i)*x(i)+ga(j)*x(j))/gaP
                ! order the angular momenta so that l1>=l2
                if (l(i)<l(j)) then
                    ii=j
                    jj=i
                else
                    ii=i
                    jj=j
                end if
                PA=Px-x(ii)
                PB=Px-x(jj)

                l1=l(ii)
                l2=l(jj)

                if (l1==0 .and. l2==0) then
                    D(ii,jj,1,1,1)=1
                    D(jj,ii,1,1,1)=1

                elseif (l1==1 .and. l2==0) then
                   D(ii,jj,1,2,1)=PA
                   D(ii,jj,2,2,1)=0.5/gaP
                   D(jj,ii,1,1,2)=D(ii,jj,1,2,1)
                   D(jj,ii,2,1,2)=D(ii,jj,2,2,1)

               elseif (l1==1 .and. l2==1) then
                    a=0.5/gaP
                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,2,1,2)=a
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a

                    D(ii,jj,1,2,2)=PB*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,2,2)=a*D(ii,jj,1,2,1)+PB*D(ii,jj,2,2,1)
                    D(ii,jj,3,2,2)=a*D(ii,jj,2,2,1)

                    D(jj,ii,1,2,2)=D(ii,jj,1,2,2)
                    D(jj,ii,2,2,2)=D(ii,jj,2,2,2)
                    D(jj,ii,3,2,2)=D(ii,jj,3,2,2)

                elseif (l1==2 .and. l2==0) then
                    a=0.5/gaP

                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a

                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)

                    D(jj,ii,1,1,3)=D(ii,jj,1,3,1)
                    D(jj,ii,2,1,3)=D(ii,jj,2,3,1)
                    D(jj,ii,3,1,3)=D(ii,jj,3,3,1)

                elseif (l1==2 .and. l2==1 ) then
                    a=0.5/gaP
                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,2,1,2)=a
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a
                    D(ii,jj,1,2,2)=PB*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,2,2)=a*D(ii,jj,1,2,1)+PB*D(ii,jj,2,2,1)
                    D(ii,jj,3,2,2)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)

                    D(ii,jj,1,3,2)=PB*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,3,2)=a*D(ii,jj,1,3,1)+PB*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,3,2)=a*D(ii,jj,2,3,1)+PB*D(ii,jj,3,3,1)
                    D(ii,jj,4,3,2)=a*D(ii,jj,3,3,1)

                    D(jj,ii,1,2,3)=D(ii,jj,1,3,2)
                    D(jj,ii,2,2,3)=D(ii,jj,2,3,2)
                    D(jj,ii,3,2,3)=D(ii,jj,3,3,2)
                    D(jj,ii,4,2,3)=D(ii,jj,4,3,2)

                elseif (l1==2 .and. l2==2) then
                    a=0.5/gaP

                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,2,1,2)=a
                    D(ii,jj,1,1,3)=PB*D(ii,jj,1,1,2)+D(ii,jj,2,1,2)
                    D(ii,jj,2,1,3)=a*D(ii,jj,1,1,2)+PB*D(ii,jj,2,1,2)
                    D(ii,jj,3,1,3)=a*D(ii,jj,2,1,2)
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a
                    D(ii,jj,1,2,2)=PB*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,2,2)=a*D(ii,jj,1,2,1)+PB*D(ii,jj,2,2,1)
                    D(ii,jj,3,2,2)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,2,3)=PB*D(ii,jj,1,2,2)+D(ii,jj,2,2,2)
                    D(ii,jj,2,2,3)=a*D(ii,jj,1,2,2)+PB*D(ii,jj,2,2,2)+2.*D(ii,jj,3,2,2)
                    D(ii,jj,3,2,3)=a*D(ii,jj,2,2,2)+PB*D(ii,jj,3,2,2)
                    D(ii,jj,4,2,3)=a*D(ii,jj,3,2,2)

                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,3,2)=PB*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,3,2)=a*D(ii,jj,1,3,1)+PB*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,3,2)=a*D(ii,jj,2,3,1)+PB*D(ii,jj,3,3,1)
                    D(ii,jj,4,3,2)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,3,3)=PB*D(ii,jj,1,3,2)+D(ii,jj,2,3,2)
                    D(ii,jj,2,3,3)=a*D(ii,jj,1,3,2)+PB*D(ii,jj,2,3,2)+2.*D(ii,jj,3,3,2)
                    D(ii,jj,3,3,3)=a*D(ii,jj,2,3,2)+PB*D(ii,jj,3,3,2)+3.*D(ii,jj,4,3,2)
                    D(ii,jj,4,3,3)=a*D(ii,jj,3,3,2)+PB*D(ii,jj,4,3,2)
                    D(ii,jj,5,3,3)=a*D(ii,jj,4,3,2)

                    D(jj,ii,1,3,3)=D(ii,jj,1,3,3)
                    D(jj,ii,2,3,3)=D(ii,jj,2,3,3)
                    D(jj,ii,3,3,3)=D(ii,jj,3,3,3)
                    D(jj,ii,4,3,3)=D(ii,jj,4,3,3)
                    D(jj,ii,5,3,3)=D(ii,jj,5,3,3)

                elseif (l1==3 .and. l2==0) then
                    a=0.5/gaP
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a
                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,4,1)=PA*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,4,1)=a*D(ii,jj,1,3,1)+PA*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,4,1)=a*D(ii,jj,2,3,1)+PA*D(ii,jj,3,3,1)
                    D(ii,jj,4,4,1)=a*D(ii,jj,3,3,1)

                    D(jj,ii,1,1,4)=D(ii,jj,1,4,1)
                    D(jj,ii,2,1,4)=D(ii,jj,2,4,1)
                    D(jj,ii,3,1,4)=D(ii,jj,3,4,1)
                    D(jj,ii,4,1,4)=D(ii,jj,4,4,1)

                elseif (l1==3 .and. l2==1) then
                    a=0.5/gaP
                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,2,1,2)=a
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a
                    D(ii,jj,1,2,2)=PB*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,2,2)=a*D(ii,jj,1,2,1)+PB*D(ii,jj,2,2,1)
                    D(ii,jj,3,2,2)=a*D(ii,jj,2,2,1)

                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,3,2)=PB*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,3,2)=a*D(ii,jj,1,3,1)+PB*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,3,2)=a*D(ii,jj,2,3,1)+PB*D(ii,jj,3,3,1)
                    D(ii,jj,4,3,2)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,4,1)=PA*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,4,1)=a*D(ii,jj,1,3,1)+PA*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,4,1)=a*D(ii,jj,2,3,1)+PA*D(ii,jj,3,3,1)
                    D(ii,jj,4,4,1)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,4,2)=PB*D(ii,jj,1,4,1)+D(ii,jj,2,4,1)
                    D(ii,jj,2,4,2)=a*D(ii,jj,1,4,1)+PB*D(ii,jj,2,4,1)+2.*D(ii,jj,3,4,1)
                    D(ii,jj,3,4,2)=a*D(ii,jj,2,4,1)+PB*D(ii,jj,3,4,1)+3.*D(ii,jj,4,4,1)
                    D(ii,jj,4,4,2)=a*D(ii,jj,3,4,1)+PB*D(ii,jj,4,4,1)
                    D(ii,jj,5,4,2)=a*D(ii,jj,4,4,1)


                    D(jj,ii,1,2,4)=D(ii,jj,1,4,2)
                    D(jj,ii,2,2,4)=D(ii,jj,2,4,2)
                    D(jj,ii,3,2,4)=D(ii,jj,3,4,2)
                    D(jj,ii,4,2,4)=D(ii,jj,4,4,2)
                    D(jj,ii,5,2,4)=D(ii,jj,5,4,2)
                    ! need to add symmetric combination

                elseif (l1==3 .and. l2==2) then
                    a=0.5/gaP
                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,2,1,2)=a
                    D(ii,jj,1,1,3)=PB*D(ii,jj,1,1,2)+D(ii,jj,2,1,2)
                    D(ii,jj,2,1,3)=a*D(ii,jj,1,1,2)+PB*D(ii,jj,2,1,2)
                    D(ii,jj,3,1,3)=a*D(ii,jj,2,1,2)
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a
                    D(ii,jj,1,2,2)=PB*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,2,2)=a*D(ii,jj,1,2,1)+PB*D(ii,jj,2,2,1)
                    D(ii,jj,3,2,2)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,2,3)=PB*D(ii,jj,1,2,2)+D(ii,jj,2,2,2)
                    D(ii,jj,2,2,3)=a*D(ii,jj,1,2,2)+PB*D(ii,jj,2,2,2)+2.*D(ii,jj,3,2,2)
                    D(ii,jj,3,2,3)=a*D(ii,jj,2,2,2)+PB*D(ii,jj,3,2,2)
                    D(ii,jj,4,2,3)=a*D(ii,jj,3,2,2)
                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,3,2)=PB*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,3,2)=a*D(ii,jj,1,3,1)+PB*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,3,2)=a*D(ii,jj,2,3,1)+PB*D(ii,jj,3,3,1)
                    D(ii,jj,4,3,2)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,3,3)=PB*D(ii,jj,1,3,2)+D(ii,jj,2,3,2)
                    D(ii,jj,2,3,3)=a*D(ii,jj,1,3,2)+PB*D(ii,jj,2,3,2)+2.*D(ii,jj,3,3,2)
                    D(ii,jj,3,3,3)=a*D(ii,jj,2,3,2)+PB*D(ii,jj,3,3,2)+3.*D(ii,jj,4,3,2)
                    D(ii,jj,4,3,3)=a*D(ii,jj,3,3,2)+PB*D(ii,jj,4,3,2)
                    D(ii,jj,5,3,3)=a*D(ii,jj,4,3,2)
                    D(ii,jj,1,4,1)=PA*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,4,1)=a*D(ii,jj,1,3,1)+PA*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,4,1)=a*D(ii,jj,2,3,1)+PA*D(ii,jj,3,3,1)
                    D(ii,jj,4,4,1)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,4,2)=PB*D(ii,jj,1,4,1)+D(ii,jj,2,4,1)
                    D(ii,jj,2,4,2)=a*D(ii,jj,1,4,1)+PB*D(ii,jj,2,4,1)+2.*D(ii,jj,3,4,1)
                    D(ii,jj,3,4,2)=a*D(ii,jj,2,4,1)+PB*D(ii,jj,3,4,1)+3.*D(ii,jj,4,4,1)
                    D(ii,jj,4,4,2)=a*D(ii,jj,3,4,1)+PB*D(ii,jj,4,4,1)
                    D(ii,jj,5,4,2)=a*D(ii,jj,4,4,1)
                    D(ii,jj,1,4,3)=PB*D(ii,jj,1,4,2)+D(ii,jj,2,4,2)
                    D(ii,jj,2,4,3)=a*D(ii,jj,1,4,2)+PB*D(ii,jj,2,4,2)+2.*D(ii,jj,3,4,2)
                    D(ii,jj,3,4,3)=a*D(ii,jj,2,4,2)+PB*D(ii,jj,3,4,2)+3.*D(ii,jj,4,4,2)
                    D(ii,jj,4,4,3)=a*D(ii,jj,3,4,2)+PB*D(ii,jj,4,4,2)+4.*D(ii,jj,5,4,2)
                    D(ii,jj,5,4,3)=a*D(ii,jj,4,4,2)+PB*D(ii,jj,5,4,2)
                    D(ii,jj,6,4,3)=a*D(ii,jj,5,4,2)


                    D(jj,ii,1,3,4)=D(ii,jj,1,4,3)
                    D(jj,ii,2,3,4)=D(ii,jj,2,4,3)
                    D(jj,ii,3,3,4)=D(ii,jj,3,4,3)
                    D(jj,ii,4,3,4)=D(ii,jj,4,4,3)
                    D(jj,ii,5,3,4)=D(ii,jj,5,4,3)
                    D(jj,ii,6,3,4)=D(ii,jj,6,4,3)

                    ! need to add symmetric combination
                elseif (l1==3 .and. l2==3) then
                    a=0.5/gaP
                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,2,1,2)=a
                    D(ii,jj,1,1,3)=PB*D(ii,jj,1,1,2)+D(ii,jj,2,1,2)
                    D(ii,jj,2,1,3)=a*D(ii,jj,1,1,2)+PB*D(ii,jj,2,1,2)
                    D(ii,jj,3,1,3)=a*D(ii,jj,2,1,2)
                    D(ii,jj,1,1,4)=PB*D(ii,jj,1,1,3)+D(ii,jj,2,1,3)
                    D(ii,jj,2,1,4)=a*D(ii,jj,1,1,3)+PB*D(ii,jj,2,1,3)+2.*D(ii,jj,3,1,3)
                    D(ii,jj,3,1,4)=a*D(ii,jj,2,1,3)+PB*D(ii,jj,3,1,3)
                    D(ii,jj,4,1,4)=a*D(ii,jj,3,1,3)
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a
                    D(ii,jj,1,2,2)=PB*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,2,2)=a*D(ii,jj,1,2,1)+PB*D(ii,jj,2,2,1)
                    D(ii,jj,3,2,2)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,2,3)=PB*D(ii,jj,1,2,2)+D(ii,jj,2,2,2)
                    D(ii,jj,2,2,3)=a*D(ii,jj,1,2,2)+PB*D(ii,jj,2,2,2)+2.*D(ii,jj,3,2,2)
                    D(ii,jj,3,2,3)=a*D(ii,jj,2,2,2)+PB*D(ii,jj,3,2,2)
                    D(ii,jj,4,2,3)=a*D(ii,jj,3,2,2)

                    D(ii,jj,1,2,4)=PB*D(ii,jj,1,2,3)+D(ii,jj,2,2,3)
                    D(ii,jj,2,2,4)=a*D(ii,jj,1,2,3)+PB*D(ii,jj,2,2,3)+2.*D(ii,jj,3,2,3)
                    D(ii,jj,3,2,4)=a*D(ii,jj,2,2,3)+PB*D(ii,jj,3,2,3)+3.*D(ii,jj,4,2,3)
                    D(ii,jj,4,2,4)=a*D(ii,jj,3,2,3)+PB*D(ii,jj,4,2,3)
                    D(ii,jj,5,2,4)=a*D(ii,jj,4,2,3)
                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,3,2)=PB*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,3,2)=a*D(ii,jj,1,3,1)+PB*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,3,2)=a*D(ii,jj,2,3,1)+PB*D(ii,jj,3,3,1)
                    D(ii,jj,4,3,2)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,3,3)=PB*D(ii,jj,1,3,2)+D(ii,jj,2,3,2)
                    D(ii,jj,2,3,3)=a*D(ii,jj,1,3,2)+PB*D(ii,jj,2,3,2)+2.*D(ii,jj,3,3,2)
                    D(ii,jj,3,3,3)=a*D(ii,jj,2,3,2)+PB*D(ii,jj,3,3,2)+3.*D(ii,jj,4,3,2)
                    D(ii,jj,4,3,3)=a*D(ii,jj,3,3,2)+PB*D(ii,jj,4,3,2)
                    D(ii,jj,5,3,3)=a*D(ii,jj,4,3,2)

                    D(ii,jj,1,3,4)=PB*D(ii,jj,1,3,3)+D(ii,jj,2,3,3)
                    D(ii,jj,2,3,4)=a*D(ii,jj,1,3,3)+PB*D(ii,jj,2,3,3)+2.*D(ii,jj,3,3,3)
                    D(ii,jj,3,3,4)=a*D(ii,jj,2,3,3)+PB*D(ii,jj,3,3,3)+3.*D(ii,jj,4,3,3)
                    D(ii,jj,4,3,4)=a*D(ii,jj,3,3,3)+PB*D(ii,jj,4,3,3)+4.*D(ii,jj,5,3,3)
                    D(ii,jj,5,3,4)=a*D(ii,jj,4,3,3)+PB*D(ii,jj,5,3,3)
                    D(ii,jj,6,3,4)=a*D(ii,jj,5,3,3)
                    D(ii,jj,1,4,1)=PA*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,4,1)=a*D(ii,jj,1,3,1)+PA*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,4,1)=a*D(ii,jj,2,3,1)+PA*D(ii,jj,3,3,1)
                    D(ii,jj,4,4,1)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,4,2)=PB*D(ii,jj,1,4,1)+D(ii,jj,2,4,1)
                    D(ii,jj,2,4,2)=a*D(ii,jj,1,4,1)+PB*D(ii,jj,2,4,1)+2.*D(ii,jj,3,4,1)
                    D(ii,jj,3,4,2)=a*D(ii,jj,2,4,1)+PB*D(ii,jj,3,4,1)+3.*D(ii,jj,4,4,1)
                    D(ii,jj,4,4,2)=a*D(ii,jj,3,4,1)+PB*D(ii,jj,4,4,1)
                    D(ii,jj,5,4,2)=a*D(ii,jj,4,4,1)

                    D(ii,jj,1,4,3)=PB*D(ii,jj,1,4,2)+D(ii,jj,2,4,2)
                    D(ii,jj,2,4,3)=a*D(ii,jj,1,4,2)+PB*D(ii,jj,2,4,2)+2.*D(ii,jj,3,4,2)
                    D(ii,jj,3,4,3)=a*D(ii,jj,2,4,2)+PB*D(ii,jj,3,4,2)+3.*D(ii,jj,4,4,2)
                    D(ii,jj,4,4,3)=a*D(ii,jj,3,4,2)+PB*D(ii,jj,4,4,2)+4.*D(ii,jj,5,4,2)
                    D(ii,jj,5,4,3)=a*D(ii,jj,4,4,2)+PB*D(ii,jj,5,4,2)
                    D(ii,jj,6,4,3)=a*D(ii,jj,5,4,2)


                    D(ii,jj,1,4,4)=PB*D(ii,jj,1,4,3)+D(ii,jj,2,4,3)
                    D(ii,jj,2,4,4)=a*D(ii,jj,1,4,3)+PB*D(ii,jj,2,4,3)+2.*D(ii,jj,3,4,3)
                    D(ii,jj,3,4,4)=a*D(ii,jj,2,4,3)+PB*D(ii,jj,3,4,3)+3.*D(ii,jj,4,4,3)
                    D(ii,jj,4,4,4)=a*D(ii,jj,3,4,3)+PB*D(ii,jj,4,4,3)+4.*D(ii,jj,5,4,3)
                    D(ii,jj,5,4,4)=a*D(ii,jj,4,4,3)+PB*D(ii,jj,5,4,3)+5.*D(ii,jj,6,4,3)
                    D(ii,jj,6,4,4)=a*D(ii,jj,5,4,3)+PB*D(ii,jj,6,4,3)
                    D(ii,jj,7,4,4)=a*D(ii,jj,6,4,3)

                    D(jj,ii,1,4,4)=D(ii,jj,1,4,4)
                    D(jj,ii,2,4,4)=D(ii,jj,2,4,4)
                    D(jj,ii,3,4,4)=D(ii,jj,3,4,4)
                    D(jj,ii,4,4,4)=D(ii,jj,4,4,4)
                    D(jj,ii,5,4,4)=D(ii,jj,5,4,4)
                    D(jj,ii,6,4,4)=D(ii,jj,6,4,4)
                    D(jj,ii,7,4,4)=D(ii,jj,7,4,4)
                    ! need to add symmetric combination
                else
                    print*, "case not programmed yet: l1/2= ", l1, l2
                    stop
                end if

            end do
        end do


        ! diagonal case

        do i = 1,N
            j=i

            gaP=ga(i)+ga(j)
            Px=(ga(i)*x(i)+ga(j)*x(j))/gaP
            PA=Px-x(i)
            PB=Px-x(j)

            l1=l(i)
            l2=l(j)

            if (l1==0 .and. l2==0) then
                D(i,j,1,1,1)=1

            elseif (l1==1) then
                a=0.5/gaP
                D(i,j,1,1,2)=PB
                D(i,j,2,1,2)=a
                D(i,j,1,2,1)=PA
                D(i,j,2,2,1)=a

                D(i,j,1,2,2)=PB*D(i,j,1,2,1)+D(i,j,2,2,1)
                D(i,j,2,2,2)=a*D(i,j,1,2,1)+PB*D(i,j,2,2,1)
                D(i,j,3,2,2)=a*D(i,j,2,2,1)

            elseif (l1==2 ) then
                a=0.5/gaP

                D(i,j,1,1,2)=PB
                D(i,j,2,1,2)=a
                D(i,j,1,1,3)=PB*D(i,j,1,1,2)+D(i,j,2,1,2)
                D(i,j,2,1,3)=a*D(i,j,1,1,2)+PB*D(i,j,2,1,2)
                D(i,j,3,1,3)=a*D(i,j,2,1,2)
                D(i,j,1,2,1)=PA
                D(i,j,2,2,1)=a
                D(i,j,1,2,2)=PB*D(i,j,1,2,1)+D(i,j,2,2,1)
                D(i,j,2,2,2)=a*D(i,j,1,2,1)+PB*D(i,j,2,2,1)
                D(i,j,3,2,2)=a*D(i,j,2,2,1)
                D(i,j,1,2,3)=PB*D(i,j,1,2,2)+D(i,j,2,2,2)
                D(i,j,2,2,3)=a*D(i,j,1,2,2)+PB*D(i,j,2,2,2)+2.*D(i,j,3,2,2)
                D(i,j,3,2,3)=a*D(i,j,2,2,2)+PB*D(i,j,3,2,2)
                D(i,j,4,2,3)=a*D(i,j,3,2,2)

                D(i,j,1,3,1)=PA*D(i,j,1,2,1)+D(i,j,2,2,1)
                D(i,j,2,3,1)=a*D(i,j,1,2,1)+PA*D(i,j,2,2,1)
                D(i,j,3,3,1)=a*D(i,j,2,2,1)
                D(i,j,1,3,2)=PB*D(i,j,1,3,1)+D(i,j,2,3,1)
                D(i,j,2,3,2)=a*D(i,j,1,3,1)+PB*D(i,j,2,3,1)+2.*D(i,j,3,3,1)
                D(i,j,3,3,2)=a*D(i,j,2,3,1)+PB*D(i,j,3,3,1)
                D(i,j,4,3,2)=a*D(i,j,3,3,1)
                D(i,j,1,3,3)=PB*D(i,j,1,3,2)+D(i,j,2,3,2)
                D(i,j,2,3,3)=a*D(i,j,1,3,2)+PB*D(i,j,2,3,2)+2.*D(i,j,3,3,2)
                D(i,j,3,3,3)=a*D(i,j,2,3,2)+PB*D(i,j,3,3,2)+3.*D(i,j,4,3,2)
                D(i,j,4,3,3)=a*D(i,j,3,3,2)+PB*D(i,j,4,3,2)
                D(i,j,5,3,3)=a*D(i,j,4,3,2)

            elseif (l1==3 ) then
                    a=0.5/gaP
                    D(i,j,1,1,2)=PB
                    D(i,j,1,1,2)=PB
                    D(i,j,2,1,2)=a
                    D(i,j,1,1,3)=PB*D(i,j,1,1,2)+D(i,j,2,1,2)
                    D(i,j,2,1,3)=a*D(i,j,1,1,2)+PB*D(i,j,2,1,2)
                    D(i,j,3,1,3)=a*D(i,j,2,1,2)
                    D(i,j,1,1,4)=PB*D(i,j,1,1,3)+D(i,j,2,1,3)
                    D(i,j,2,1,4)=a*D(i,j,1,1,3)+PB*D(i,j,2,1,3)+2.*D(i,j,3,1,3)
                    D(i,j,3,1,4)=a*D(i,j,2,1,3)+PB*D(i,j,3,1,3)
                    D(i,j,4,1,4)=a*D(i,j,3,1,3)
                    D(i,j,1,2,1)=PA
                    D(i,j,2,2,1)=a
                    D(i,j,1,2,2)=PB*D(i,j,1,2,1)+D(i,j,2,2,1)
                    D(i,j,2,2,2)=a*D(i,j,1,2,1)+PB*D(i,j,2,2,1)
                    D(i,j,3,2,2)=a*D(i,j,2,2,1)
                    D(i,j,1,2,3)=PB*D(i,j,1,2,2)+D(i,j,2,2,2)
                    D(i,j,2,2,3)=a*D(i,j,1,2,2)+PB*D(i,j,2,2,2)+2.*D(i,j,3,2,2)
                    D(i,j,3,2,3)=a*D(i,j,2,2,2)+PB*D(i,j,3,2,2)
                    D(i,j,4,2,3)=a*D(i,j,3,2,2)

                    D(i,j,1,2,4)=PB*D(i,j,1,2,3)+D(i,j,2,2,3)
                    D(i,j,2,2,4)=a*D(i,j,1,2,3)+PB*D(i,j,2,2,3)+2.*D(i,j,3,2,3)
                    D(i,j,3,2,4)=a*D(i,j,2,2,3)+PB*D(i,j,3,2,3)+3.*D(i,j,4,2,3)
                    D(i,j,4,2,4)=a*D(i,j,3,2,3)+PB*D(i,j,4,2,3)
                    D(i,j,5,2,4)=a*D(i,j,4,2,3)
                    D(i,j,1,3,1)=PA*D(i,j,1,2,1)+D(i,j,2,2,1)
                    D(i,j,2,3,1)=a*D(i,j,1,2,1)+PA*D(i,j,2,2,1)
                    D(i,j,3,3,1)=a*D(i,j,2,2,1)
                    D(i,j,1,3,2)=PB*D(i,j,1,3,1)+D(i,j,2,3,1)
                    D(i,j,2,3,2)=a*D(i,j,1,3,1)+PB*D(i,j,2,3,1)+2.*D(i,j,3,3,1)
                    D(i,j,3,3,2)=a*D(i,j,2,3,1)+PB*D(i,j,3,3,1)
                    D(i,j,4,3,2)=a*D(i,j,3,3,1)
                    D(i,j,1,3,3)=PB*D(i,j,1,3,2)+D(i,j,2,3,2)
                    D(i,j,2,3,3)=a*D(i,j,1,3,2)+PB*D(i,j,2,3,2)+2.*D(i,j,3,3,2)
                    D(i,j,3,3,3)=a*D(i,j,2,3,2)+PB*D(i,j,3,3,2)+3.*D(i,j,4,3,2)
                    D(i,j,4,3,3)=a*D(i,j,3,3,2)+PB*D(i,j,4,3,2)
                    D(i,j,5,3,3)=a*D(i,j,4,3,2)

                    D(i,j,1,3,4)=PB*D(i,j,1,3,3)+D(i,j,2,3,3)
                    D(i,j,2,3,4)=a*D(i,j,1,3,3)+PB*D(i,j,2,3,3)+2.*D(i,j,3,3,3)
                    D(i,j,3,3,4)=a*D(i,j,2,3,3)+PB*D(i,j,3,3,3)+3.*D(i,j,4,3,3)
                    D(i,j,4,3,4)=a*D(i,j,3,3,3)+PB*D(i,j,4,3,3)+4.*D(i,j,5,3,3)
                    D(i,j,5,3,4)=a*D(i,j,4,3,3)+PB*D(i,j,5,3,3)
                    D(i,j,6,3,4)=a*D(i,j,5,3,3)
                    D(i,j,1,4,1)=PA*D(i,j,1,3,1)+D(i,j,2,3,1)
                    D(i,j,2,4,1)=a*D(i,j,1,3,1)+PA*D(i,j,2,3,1)+2.*D(i,j,3,3,1)
                    D(i,j,3,4,1)=a*D(i,j,2,3,1)+PA*D(i,j,3,3,1)
                    D(i,j,4,4,1)=a*D(i,j,3,3,1)
                    D(i,j,1,4,2)=PB*D(i,j,1,4,1)+D(i,j,2,4,1)
                    D(i,j,2,4,2)=a*D(i,j,1,4,1)+PB*D(i,j,2,4,1)+2.*D(i,j,3,4,1)
                    D(i,j,3,4,2)=a*D(i,j,2,4,1)+PB*D(i,j,3,4,1)+3.*D(i,j,4,4,1)
                    D(i,j,4,4,2)=a*D(i,j,3,4,1)+PB*D(i,j,4,4,1)
                    D(i,j,5,4,2)=a*D(i,j,4,4,1)

                    D(i,j,1,4,3)=PB*D(i,j,1,4,2)+D(i,j,2,4,2)
                    D(i,j,2,4,3)=a*D(i,j,1,4,2)+PB*D(i,j,2,4,2)+2.*D(i,j,3,4,2)
                    D(i,j,3,4,3)=a*D(i,j,2,4,2)+PB*D(i,j,3,4,2)+3.*D(i,j,4,4,2)
                    D(i,j,4,4,3)=a*D(i,j,3,4,2)+PB*D(i,j,4,4,2)+4.*D(i,j,5,4,2)
                    D(i,j,5,4,3)=a*D(i,j,4,4,2)+PB*D(i,j,5,4,2)
                    D(i,j,6,4,3)=a*D(i,j,5,4,2)

                    D(i,j,1,4,4)=PB*D(i,j,1,4,3)+D(i,j,2,4,3)
                    D(i,j,2,4,4)=a*D(i,j,1,4,3)+PB*D(i,j,2,4,3)+2.*D(i,j,3,4,3)
                    D(i,j,3,4,4)=a*D(i,j,2,4,3)+PB*D(i,j,3,4,3)+3.*D(i,j,4,4,3)
                    D(i,j,4,4,4)=a*D(i,j,3,4,3)+PB*D(i,j,4,4,3)+4.*D(i,j,5,4,3)
                    D(i,j,5,4,4)=a*D(i,j,4,4,3)+PB*D(i,j,5,4,3)+5.*D(i,j,6,4,3)
                    D(i,j,6,4,4)=a*D(i,j,5,4,3)+PB*D(i,j,6,4,3)
                    D(i,j,7,4,4)=a*D(i,j,6,4,3)

            else

                    print*, "case not programmed yet: l1/2= " , l1, l2
                    stop
            end if
        end do

    END SUBROUTINE fill_md_table

    SUBROUTINE Bubble_Sort(a)
        use types
        implicit none
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

        use types
       IMPLICIT NONE
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

END MODULE twordmreader


