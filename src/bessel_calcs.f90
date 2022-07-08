module bessel_calcs
    implicit none
contains
    SUBROUTINE BesselDeriv(BD, LL, MM,NN,a,b,c,LLmax)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        ! The three nested loops give the formula
        ! in the form of coefficients multiplying the h functions
        REAL(kind=dp), INTENT(out), DIMENSION(LLmax+1)  :: BD
        INTEGER(kind=ikind), INTENT(in)                 :: LL, MM, NN, LLmax
        REAL(kind=dp), INTENT(in), DIMENSION(13,13)  :: a, b, c

        ! loop and temp variables
        INTEGER(kind=ikind) :: ii, jj, kk, hOrder, temp, ceil
        REAL(kind=dp)       :: C1, C2, C3, Ct2, Ct3
        ! set this to 0 initially and accumulate
        BD = 0.0_dp
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

    SUBROUTINE BesselDeriv1j(BD, hz,h,LL,MM,NN, a, b, c,LLmax)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        ! The three nested loops give the formula
        ! in the form of coefficients multiplying the h functions
        REAL(kind=dp), INTENT(out), DIMENSION(13)  :: BD
        INTEGER(kind=ikind), INTENT(in)                 :: LL, MM, NN, LLmax
        REAL(kind=dp), INTENT(in), DIMENSION(13,13)  :: a, b, c
        real(kind=dp),intent(in) :: h,hz
        ! loop and temp variables
        INTEGER(kind=ikind) :: ii, jj, kk, hOrder, temp, ceil
        REAL(kind=dp)       :: C1, C2, C3,pi,eta



        pi=dacos(-1.0_dp)
        eta=3.0_dp*hz**2-h**2

        BD = 0.0_dp

        do ii = 0, LL+2
            C1=a(LL+1,ii+1)




            do jj = 0, MM+2
                C2=b(MM+1,jj+1)




                do kk = 0, NN+2
                    C3=c(NN+1,kk+1)


                    C3 = c1*c2*c3

                    temp = CEILING((LL+MM+NN+ii+jj+kk)/2.0_dp-1.0_dp)

                    hOrder = temp

                    BD(hOrder+1)=BD(hOrder+1) + C3

                end do
            end do
        end do

        BD=0.5_dp * (3.0_dp/pi)**(0.5_dp)  * BD


    END SUBROUTINE

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
        real(kind=dp), dimension(0:18, size(mu)):: allbessels, allbessels2
        allbessels=0.0_dp
        sum=0.0_dp
        Pmu=H * mu



        do i=1,size(mu)
            if (abs(pmu(i))<0.05) then
                call van(allbessels(0:6,i),6,pmu(i))
                allbessels(7:18,i)=0.0_dp
            else
                call van(allbessels(0:18,i),18,pmu(i))
            end if

        enddo

        !  print*,(allbessels(0:18,1))
        ! stop
        do beta=0,18

            hqr = mu**(beta-2) / h**beta
            bess= allbessels(beta,:)

            sum=sum+bess*hqr*h_saved(beta+1)



        enddo





    End Subroutine bessels2rr

    Subroutine bessels1rr(sum,order,mu,H, h_saved)
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
        real(kind=dp), dimension(0:16, size(mu)):: allbessels, allbessels2
        allbessels=0.0_dp
        sum=0.0_dp
        Pmu=H * mu



        do i=1,size(mu)
            if (abs(pmu(i))<0.05) then
                call van(allbessels(0:6,i),6,pmu(i))
                allbessels(7:16,i)=0.0_dp
            else
                call van(allbessels(0:16,i),16,pmu(i))
            end if

        enddo

        !  print*,(allbessels(0:18,1))
        ! stop
        do beta=0,16

            hqr = mu**(beta-1) / h**beta
            bess= allbessels(beta,:)

            sum=sum+bess*hqr*h_saved(beta+1)



        enddo





    End Subroutine bessels1rr

    Subroutine bessels0rr(sum,order,mu,H, h_saved)
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
        real(kind=dp), dimension(0:16, size(mu)):: allbessels, allbessels2
        allbessels=0.0_dp
        sum=0.0_dp
        Pmu=H * mu



        do i=1,size(mu)
            if (abs(pmu(i))<0.05) then
                call van(allbessels(0:6,i),6,pmu(i))
                allbessels(7:18,i)=0.0_dp
            else
                call van(allbessels(0:16,i),16,pmu(i))
            end if

        enddo

        !  print*,(allbessels(0:18,1))
        ! stop
        do beta=0,order
            hqr = (mu / h)**beta
            bess= allbessels(beta,:)
            sum=sum+bess*hqr*h_saved(beta+1)
        enddo





    End Subroutine bessels0rr

    subroutine van (jbwd,n, x)
        implicit none
        double precision, intent(in):: x
        integer, intent(in) :: n
        integer i, nmax
        double precision, intent(out),dimension(0:N):: jbwd
        double precision ::  ix, jmx1, jmx, j0, j1


        nmax= max(n, ceiling(1.142476370122814*n-4.027048776987268))


        ix = 1 / x
        j0 = sin(x) * ix
        j1 = (j0 - cos(x)) * ix
        !
        ! Backward recursion.
        !
        jmx1 = 0
        jmx = 1
        do i = n + nmax, n + 2, -1
            jbwd(0) = (2 * i + 1) * (ix * jmx - jmx1 / (2 * i + 1))
            jmx1 = jmx
            jmx = jbwd(0)
        end do
        do i = n + 1, 1, -1
            jbwd(i - 1) = (2 * i + 1) * ix * jmx - jmx1
            jmx1 = jmx
            jmx = jbwd(i - 1)
        end do
        if (abs(jmx) >= abs(jmx1)) then
            j0 = j0 / jbwd(0)
        else
            j0 = j1 / jbwd(1)
        end if
        jbwd = j0 * jbwd
    end

    SUBROUTINE rrdj1xy(lmax,x,a)
        ! This subroutine uses a recurrence relation to calculate the
        ! first set of expansion coefficients, a, at given coordinate
        ! x (or y) up to an angular momentum of L

        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer, parameter    :: dp = selected_real_kind(15)
        integer, intent(in)                                        :: lmax
        integer                                                    :: l,p
        real(kind=dp), intent(in)                                  :: x
        real(kind=dp), dimension(:,:), intent(out) :: a


        a(1,2) = 1.0_dp
        a(2,3) = -x

        do l = 2, lmax
            a(l+1,2) = -(l-1) * a(l-1,2)
        enddo

        do l = 2, lmax
            do p = 1, l
                a(l+1,p+2) = -(l-1) * a(l-1,p+2) - x * a(l,p+1)
            enddo
        enddo

    END SUBROUTINE

    SUBROUTINE rrdj1z(nmax,z,a)
        ! This subroutine uses a recurrence relation to calculate the
        ! first set of expansion coefficients, a, at given coordinate
        ! z up to an angular momentum of L

        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer, parameter    :: dp = selected_real_kind(15)
        integer, intent(in)                                        :: nmax
        integer                                                    :: n,s
        real(kind=dp), intent(in)                                  :: z
        real(kind=dp), dimension(:,:), allocatable                 :: a0
        real(kind=dp), dimension(:,:), intent(out)   :: a

        allocate( a0(nmax+2,nmax+3) )


        a0(1,2) = -1.0_dp
        a0(2,3) = z

        do n = 2, nmax+1
            a0(n+1,2) = -(n-1) * a0(n-1,2)
        enddo

        do n = 2, nmax+1
            do s = 1, n
                a0(n+1,s+2) = -(n-1) * a0(n-1,s+2) - z * a0(n,s+1)
            enddo
        enddo

        do n = 1, nmax+1
            do s = 1, n+1
                a(n,s) = a0(n+1,s+1)
            enddo
        enddo

    END SUBROUTINE

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

    SUBROUTINE Hermite_like_coeffs(a, LLmax, Hx)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        REAL(kind=dp), INTENT(out), DIMENSION(13,13) :: a
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


end module bessel_calcs