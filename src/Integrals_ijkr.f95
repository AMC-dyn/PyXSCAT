module integrals_ijkr

    implicit none 

    contains 


    subroutine tot_integral_k_ijkr(mu,lmax1,lmax2,lmax3,lmax4,hx,hy,hz,h,dx, dy, dz, i,j, k, r, z, z2, apos, cutoffz, &
                cutoffmd,int_res)

        implicit none 

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        integer(kind=selected_int_kind(8)), intent(in)  :: lmax1,lmax2,lmax3,lmax4,i,j,k,r
        integer(kind=selected_int_kind(8)), dimension(:), intent(in) :: apos
        real(kind=dp), intent(in)              :: cutoffz,cutoffmd, hx, hy, hz, h
        real(kind=dp), intent(in), dimension(:,:,:) :: z, z2
        real(kind=dp), intent(in), dimension(:)   ::  mu
        real(kind=dp), intent(in),  dimension(:,:,:,:,:) :: dx,dy,dz



        integer(kind=selected_int_kind(8)) :: llmax, l, m, n, ka, posi, posj, posk, posr, ra
        integer(kind=selected_int_kind(8)) ::  l1, l2, l3, l4, m1, m2, m3, m4, n1, n2, n3, n4, h1
        integer(kind=selected_int_kind(8)) ::  ll, mm, nn, lp, mp, np

        real(kind=dp)  ::coeff,prodd,ztot,mdn, mdl, mdm, mdlp, mdmp
        real(kind=dp), dimension(:,:), allocatable :: a, b, c
        real(kind=dp), dimension(:,:,:,:), allocatable :: h_pre2
        real(kind=dp), dimension(:), allocatable :: h_saved, bd, pmu, h_sum, h_0, h_1, h_r, muoh,zij, zij2, zkr, zkr2

        real(kind=dp), intent(out), dimension(size(mu)) :: int_res





        llmax=lmax1+lmax2+lmax3+lmax4

        allocate(zij(size(Z(:,1,1))), zij2(size(Z(:,1,1))), zkr(size(Z(:,1,1))), zkr2(size(Z(:,1,1))))



        allocate(a(llmax+1,llmax+1),b(llmax+1,llmax+1), c(llmax+1,llmax+1),bd(llmax+1))
        a(1,1)=1
        b(1,1)=1
        c(1,1)=1

        if (llmax==0) then 
            a(2,2)=-hx
            b(2,2)=-hy
            c(2,2)=-hz
        endif

        do l=2,llmax
            a(l+1,1)=-(l-1)*a(l-1,1)
            b(l+1,1)=-(l-1)*b(l-1,1)
            c(l+1,1)=-(l-1)*c(l-1,1)

            do ka=1,llmax
                a(l+1,ka+1)=-hx*a(l,ka)-(l-1)*a(l-1,ka+1)
                b(l+1,ka+1)=-hy*b(l,ka)-(l-1)*b(l-1,ka+1)
                c(l+1,ka+1)=-hz*c(l,ka)-(l-1)*c(l-1,ka+1)
            enddo
        enddo

        allocate(h_saved(llmax+1),h_pre2(llmax+1,llmax+1,llmax+1,llmax+1))

        do l=0,llmax
            do m=0,llmax-1
                do n=0,llmax-l-m
                    call besselderiv(bd,l, m,n,a,b,c,llmax)
               
                    h_pre2(:,l+1,m+1,n+1) = bd
                enddo
            enddo
        enddo
        posi=apos(i+1)


! loop through all possible ways to get total angular momentum lmax1
        do l1=0,lmax1
            do m1=0,(lmax1-l1)
                n1=lmax1-l1-m1
                posj=apos(j+1)
        !loop through all possible ways to get total angular momentum lmax2
                do l2=0,lmax2
                    do m2=0,(lmax2-l2)
                        n2=lmax2-l2-m2


                        zij(:)=z(:,posi,posj)
                        zij2(:)=z2(:,posi,posj)

                        posk=apos(k+1)
                        do l3=0,lmax3
                            do m3=0,(lmax3-l3)
                                n3=lmax3-l3-m3
                                posr=apos(r+1)

                        ! loop through all possible ways to get total angular momentum lmax4
                                do l4=0,lmax4
                                    do m4=0,(lmax4-l4)
                                        n4=lmax4-l4-m4;

                                        zkr=z(:,posk,posr)
                                        zkr2=z2(:,posk,posr)
                                ! total prefactor

                                        ztot=sum(zij*zkr2+zij2*zkr)/8.d0

                                        if (abs(ztot)<cutoffz) then 
                                            posr=posr+1
                                            cycle
                                        endif


                                        do l=0,(l1+l2)
                                            mdl=dx(i+1,j+1,l+1,l1+1,l2+1)*ztot
                                   
                                            if (mdl==0) cycle
                                                 
                                         
                                            do m=0,(m1+m2)
                                                mdm=dy(i+1,j+1,m+1,m1+1,m2+1)*mdl
                                                if (mdm==0) cycle
                                         
                                                do n=0,(n1+n2)
                                                    h1=(-1)**(l+m+n)
                                                    mdn=dz(i+1,j+1,n+1,n1+1,n2+1)*mdm*h1
                                                    if (mdn==0) cycle  
                                                    
                                                    do lp=0,(l3+l4)
                                                        mdlp=dx(k+1,r+1,lp+1,l3+1,l4+1)*mdn
                                                        if (mdlp==0) cycle 

                                                    
                                                        ll=l+lp+1;
                                                        do mp=0,(m3+m4)
                                                            mdmp=dy(k+1,r+1,mp+1,m3+1,m4+1)*mdlp
                                                            if (mdmp==0) cycle
                                                     
                                                            mm=m+mp+1
                                                            do np=0,(n3+n4)
                                                                prodd=dz(k+1,r+1,np+1,n3+1,n4+1)*mdmp
                                                                if (abs(prodd)<cutoffmd) cycle

                                                                nn=n+np+1
                                                                h_saved=h_saved+h_pre2(:,ll,mm,nn)*prodd
                                                            enddo
                                                        enddo
                                                    enddo
                                                enddo
                                            enddo
                                        enddo
                                        posr=posr+1
                                    enddo
                                enddo
                                posk=posk+1
                            enddo
                        enddo
                        posj=posj+1
                    enddo
                enddo
                posi=posi+1
            enddo
        enddo


        allocate(pmu(size(mu)), muoh(size(mu)), h_0(size(mu)),h_1(size(mu)), h_r(size(mu)))
        pmu=h*mu


        h_0=sin(pmu)/pmu
        coeff=h_saved(1)
        h_sum=h_0*coeff

        if (llmax==1) then
            h_1=(sin(pmu)/pmu**2-cos(pmu)/pmu)*mu/h 
            coeff=h_saved(2)
            h_sum=h_sum+coeff*h_1
        elseif (llmax>1) then
            muoh=mu/h
            h_1=(sin(pmu)/pmu**2-cos(pmu)/pmu)*muoh
            coeff=h_saved(2)
            h_sum=h_sum+coeff*h_1
            do ra=2,llmax
                coeff=h_saved(ra+1)
                h_r= ((2*ra-1)/(pmu)*h_1-h_0*muoh)*muoh
                h_sum=h_sum+h_r*coeff
                h_0=h_1
                h_1=h_r
            enddo
        endif
        int_res=h_sum
    end subroutine




    subroutine besselderiv(bd, ll, mm,nn,a,b,c,llmax)
        ! the three nested loops give the formula 
        ! in the form of coefficients multiplying the h functions
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        real(kind=dp), intent(out), dimension(llmax+1)  :: bd
        integer(kind=selected_int_kind(8)), intent(in)                 :: ll, mm, nn, llmax
        real(kind=dp), intent(in), dimension(llmax+1,llmax+1)  :: a, b, c

        ! loop and temp variables
        integer(kind=selected_int_kind(8)) :: ii, jj, kk, horder, temp, ceil
        real(kind=dp)       :: c1, c2, c3, ct2, ct3
        ! set this to 0 initially and accumulate
        bd = 0 
        !tempbesselpre=zeros(llmax+1,1);
        do ii = 0, ll
            c1=a(ll+1,ii+1)
            if (abs(c1)<1.0e-30) cycle

            do jj = 0, mm
                ct2=b(mm+1,jj+1)
                if (abs(ct2)<1.0e-30) cycle
                c2 = c1 * ct2
        
                do kk = 0, nn
                    ct3=c(nn+1,kk+1)
                    if (abs(ct3)<1.0e-30) cycle
                    c3 = c2 * ct3
                    ! horder = ceiling((ll+mm+nn-ii-jj-kk)/2.0_dp)+ii+jj+kk
                    temp = ll+mm+nn-ii-jj-kk
                    ceil = temp/2 + mod(temp, 2) ! integer division with rounding up
                    horder = ceil+ii+jj+kk
                    bd(horder+1)=bd(horder+1) + c3
                end do
            end do
        end do
        end subroutine
    end module integrals_ijkr