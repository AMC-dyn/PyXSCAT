! Created by  on 26/09/2023.

module Bessels_j0
implicit none
    contains
    subroutine bessels0rr(h_sum, order, mu, H, h_saved)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer :: order2
        integer(kind = ikind) :: i, beta, ra, maxord, ind
        integer(kind = ikind), intent(in) :: order
        real (kind = dp), dimension(:), intent(in) :: h_saved
        real (kind = dp), dimension(:), intent(in), allocatable :: mu
        real (kind = dp), intent(in) :: H
        real(kind = dp) :: qrad, hfunc2, sumpmu, coeff, bess1, bess2, bessnew, time1, time2
        real(kind = dp), dimension(size(mu)) :: pmu, muOH, h_2, h_1, h_0, h_r, hqr, bess, sinpmu, cospmu, invpmu
        real(kind = dp), dimension(size(mu)) :: h_3, h_4, h_5
        real(kind = dp), dimension(size(mu)), intent(out) :: h_sum
        real(kind = dp), dimension(100, size(mu)) :: allbessels, allbessels2
        allbessels = 0.0_dp

        Pmu = H * mu
        invpmu = mu / H
        sinpmu = sin(pmu)
        cospmu = cos(pmu)
        h_sum = sinpmu / pmu * h_saved(1)

        if (order==1) then
            h_1 = (sinpmu / Pmu**2 - cospmu / Pmu) * invpmu ! could be written in more efficient way by deleting mu

            h_sum = h_sum + h_saved(2) * h_1

        elseif (order==2) then
            h_2 = ((3.d0 / (Pmu)**3 - 1.d0 / Pmu) * sinpmu - 3.d0 / (pmu)**2 * cospmu) * invpmu**2
            h_1 = (sinpmu / Pmu**2 - cospmu / Pmu) * invpmu
            h_sum = h_sum + h_saved(2) * h_1 + h_saved(3) * h_2
        elseif(order==3) then

            h_3 = (pmu * (-15.0d0 + pmu**2.0d0) * cospmu + 3.0d0 * (5.0d0 - 2.0d0 * pmu**2.0d0) * sinpmu) / pmu**4.0d0 * invpmu**3.0d0
            h_2 = ((3.d0 / (Pmu)**3 - 1.d0 / Pmu) * sinpmu - 3.d0 / (pmu)**2.0d0 * cospmu) * invpmu**2.0d0
            h_1 = (sinpmu / Pmu**2 - cospmu / Pmu) * invpmu
            h_sum = h_sum + h_saved(2) * h_1 + h_saved(3) * h_2 + h_saved(4) * h_3
            !        elseif(order==4) then
            !            h_4=(5.0d0*pmu*(-21.0d0 + 2.0d0*pmu**2.0d0)*cospmu + (105.0d0 - 45.0d0*pmu**2.0d0 &
            !&             + pmu**4)*sinpmu)/pmu**5*invpmu**4
            !            h_3 = (pmu*(-15.0d0 + pmu**2.0d0)*cospmu + 3.0d0*(5.0d0 - 2.0d0*pmu**2.0d0)*sinpmu)/pmu**4.0d0*invpmu**3.0d0
            !            h_2=((3.d0/(Pmu)**3-1.d0/Pmu)*sinpmu-3.d0/(pmu)**2*cospmu)*invpmu**2
            !            h_1=(sinpmu/Pmu**2-cospmu/Pmu)*invpmu
            !            h_sum=h_sum+h_saved(2)*h_1+h_saved(3)*h_2+h_saved(4)*h_3+h_saved(5)*h_4
            !        elseif (order==5) then
            !             h_5 = (-(pmu*(945.0d0 - 105.0d0*pmu**2.0d0 + pmu**4.0d0)*cospmu) + 15.0d0*(63.0d0 - 28.0d0*pmu**2.0d0 &
            !&             + pmu**4.0d0)*sinpmu)/pmu**6.0d0
            !            h_4=(5.0d0*pmu*(-21.0d0 + 2.0d0*pmu**2.0d0)*cospmu + (105.0d0 - 45.0d0*pmu**2.0d0 &
            !&             + pmu**4.0d0)*sinpmu)/pmu**5.0d0*invpmu**4.0d0
            !            h_3 = (pmu*(-15.0d0 + pmu**2.0d0)*cospmu + 3.0d0*(5.0d0 - 2.0d0*pmu**2.0d0)*sinpmu)/pmu**4.0d0*invpmu**3.0d0
            !            h_2=((3.d0/(Pmu)**3-1.d0/Pmu)*sinpmu-3.d0/(pmu)**2.0d0*cospmu)*invpmu**2.0d0
            !            h_1=(sinpmu/Pmu**2-cospmu/Pmu)*invpmu
            !            h_sum=h_sum+h_saved(2)*h_1+h_saved(3)*h_2+h_saved(4)*h_3+h_saved(5)*h_4+h_5*h_saved(6)
        elseif(order>=4) then
            !            h_5 = (-(pmu*(945.0d0 - 105.0d0*pmu**2.0d0 + pmu**4.0d0)*cospmu) + 15.0d0*(63.0d0 - 28.0d0*pmu**2.0d0 &
            !&             + pmu**4.0d0)*sinpmu)/pmu**6.0d0
            !            h_4=(5.0d0*pmu*(-21.0d0 + 2.0d0*pmu**2.0d0)*cospmu + (105.0d0 - 45.0d0*pmu**2.0d0 &
            !&             + pmu**4.0d0)*sinpmu)/pmu**5.0d0*invpmu**4.0d0
            !            h_3 = (pmu*(-15.0d0 + pmu**2.0d0)*cospmu + 3.0d0*(5.0d0 - 2.0d0*pmu**2.0d0)*sinpmu)/pmu**4.0d0*invpmu**3.0d0
            !            h_2=((3.d0/(Pmu)**3-1.d0/Pmu)*sinpmu-3.d0/(pmu)**2.0d0*cospmu)*invpmu**2.0d0
            !            h_1=(sinpmu/Pmu**2-cospmu/Pmu)*invpmu
            !            h_sum=h_sum+h_saved(2)*h_1+h_saved(3)*h_2+h_saved(4)*h_3
            call dphrec(pmu, allbessels, 100, order + 1, size(mu))
            do beta = 1, order
                h_sum = h_sum + allbessels(beta + 1, :) * invpmu**beta * h_saved(beta + 1)
            enddo

        end if


        !         do i=1,size(mu)
        !                if (abs(pmu(i))<0.05) then
        !                    call van(allbessels(0:6,i),6,pmu(i))
        !                    allbessels(7:40,i)=0.0_dp
        !                elseif(abs(pmu(i))>100) then
        !                    print*,'oh my goood'
        !                else
        !                    call van(allbessels(0:40,i),20,pmu(i))
        !                 end if
        !
        !         enddo
        !
        !       !  print*,(allbessels(0:18,1))
        !       ! stop
        !        h_sum=0.0_dp
        !        do beta=0,order
        !            hqr = (mu / h)**beta
        !            bess= allbessels(beta,:)
        !            h_sum=h_sum+bess*hqr*h_saved(beta+1)
        !        enddo

    end subroutine bessels0rr
    subroutine rrdj0(lmax, x, a)
        ! This subroutine uses a recurrence relation to calculate the
        ! first set of expansion coefficients, a, at given coordinate
        ! x (or y or z) up to an angular momentum of L

        integer, parameter :: dp = selected_real_kind(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer, intent(in) :: lmax
        integer(kind = ikind) :: l, p
        real(kind = dp), intent(in) :: x
        real(kind = dp), dimension(20, 20), intent(out) :: a

        a = 0.0_dp

        a(1, 1) = 1.0_dp
        a(2, 2) = -x

        do l = 2, lmax
            a(l + 1, 1) = -(l - 1) * a(l - 1, 1)
        enddo

        do l = 2, lmax
            do p = 1, l
                a(l + 1, p + 1) = -(l - 1) * a(l - 1, p + 1) - x * a(l, p)
            enddo
        enddo

    end subroutine
    subroutine dphrec(x, rbes, leps1, lmax1, nq)

        IMPLICIT REAL(kind = 8)(A-H, o-Z)
        implicit integer(kind = 8)(i-n)
        DIMENSION RBES(LEPS1, nq), RNEU(LEPS1, nq)
        DIMENSION X(nq), XX(nq), CX(nq), DX(nq), CU(nq), A(nq), SX(nq)

        XX = X * X
        !LEPS=0.5D0*dsqrt(XX/EPSL**(1.D0/3.00)+9.00) + 0.5D0
        !IF(LEPS>=LEPS1) GOTO 101

        rbes = 0.0d0
        rneu = 0.0d0
        RBES(LEPS1, :) = 1.0D0
        RBES(LEPS1 - 1, :) = 1.0D0
        CX = DCOS(X)
        SX = DSIN(X)
        RNEU(1, :) = CX
        RNEU(2, :) = CX + X * SX

        DO K = 3, LEPS1
            LU = LEPS1 - K + 2
            RBES(LU - 1, :) = RBES(LU, :) - XX / (4.D0 * LU * LU - 1.D0) * RBES(LU + 1, :)
        end do
        A = RBES(1, :) * RNEU(2, :) - XX / 3.0d0 * RBES(2, :) * CX
        DO K = 1, LEPS1
            RBES(K, :) = RBES(K, :) / A
        end do
        CU = 1.0D0 / X

        DO K = 1, LMAX1
            W = 2.D0 * (K - 1)
            CU = X / (W + 1.D0) * CU
            RBES(K, :) = CU * RBES(K, :)
        end do

        message = 0
    END
    subroutine BesselDeriv(BD, LL, MM, NN, a, b, c, LLmax)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        ! The three nested loops give the formula
        ! in the form of coefficients multiplying the h functions
        REAL(kind = dp), INTENT(out), DIMENSION(LLmax + 1) :: BD
        INTEGER(kind = ikind), INTENT(in) :: LL, MM, NN, LLmax
        REAL(kind = dp), INTENT(in), DIMENSION(20, 20) :: a, b, c

        ! loop and temp variables
        INTEGER(kind = ikind) :: ii, jj, kk, hOrder, temp, ceil
        REAL(kind = dp) :: C1, C2, C3, Ct2, Ct3
        ! set this to 0 initially and accumulate
        BD = 0.0_dp
        !TempBesselPre=zeros(LLmax+1,1);
        do ii = 0, LL
            C1 = a(LL + 1, ii + 1)
            if (abs(C1)<1.0E-30) cycle

            do jj = 0, MM
                Ct2 = b(MM + 1, jj + 1)
                if (abs(Ct2)<1.0E-30) cycle
                C2 = C1 * Ct2

                do kk = 0, NN
                    Ct3 = c(NN + 1, kk + 1)
                    if (abs(Ct3)<1.0E-30) cycle
                    C3 = C2 * Ct3
                    ! hOrder = ceiling((LL+MM+NN-ii-jj-kk)/2.0_dp)+ii+jj+kk
                    temp = LL + MM + NN - ii - jj - kk
                    ceil = temp / 2_ikind + mod(temp, 2_ikind) ! integer division with rounding up
                    hOrder = ceil + ii + jj + kk
                    ! Barrilero style:
                    ! hOrder = ceiling((LL+MM+NN+ii+jj+kk)/2_ikind)
                    BD(hOrder + 1) = BD(hOrder + 1) + C3
                end do
            end do
        end do
    END SUBROUTINE

    Subroutine bessel0sum(h_sum,order,BesselQ,BesselQ2,mu,H,h_saved)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer:: order2
        real (kind=dp):: fact
        integer(kind=ikind)::i,beta,ra,maxord,ind
        integer(kind=ikind), intent(in) :: order
        real (kind=dp), dimension(:), intent(in) :: h_saved
        real (kind=dp), dimension(:), intent(in), allocatable :: mu
        real(kind=dp), dimension(:,:), intent(in)::BesselQ,BesselQ2
        real (kind=dp),  intent(in) :: H
        real(kind=dp):: qrad,hfunc2,sumpmu,coeff,bess1,bess2,bessnew,time1,time2
        real(kind=dp),dimension(size(mu)):: pmu, muOH,h_2,h_1,h_0,h_r,hqr,bess,sinpmu,cospmu,invpmu,ratan
        real(kind=dp),dimension(size(mu)):: pmu2,pmu3
        real(kind=dp),dimension(size(mu))::h_3,h_4,h_5
        real(kind=dp), dimension(size(mu)), intent(out) :: h_sum
        real(kind=dp), dimension(100, size(mu)):: allbessels, allbessels2
        allbessels=0.0_dp

            Pmu=H * mu
            Pmu2=Pmu*Pmu
            Pmu3=Pmu*Pmu*Pmu
            invpmu=mu/H
            sinpmu=sin(pmu)
            cospmu=cos(pmu)
            h_sum=sinpmu/pmu*h_saved(1)

        if (order==1) then
        h_1=(sinpmu/pmu2-cospmu/Pmu)*invpmu ! could be written in more efficient way by deleting mu

        h_sum=h_sum+h_saved(2)*h_1

        elseif (order==2) then
            h_2=((3.d0/pmu3-1.d0/Pmu)*sinpmu-3.d0/pmu2*cospmu)*invpmu**2
            h_1=(sinpmu/Pmu2-cospmu/Pmu)*invpmu
            h_sum=h_sum+h_saved(2)*h_1+h_saved(3)*h_2
        elseif(order==3) then


            h_3 = (pmu*(-15.0d0 + pmu**2.0d0)*cospmu + 3.0d0*(5.0d0 - 2.0d0*pmu**2.0d0)*sinpmu)/pmu**4.0d0*invpmu**3.0d0
            h_2=((3.d0/(Pmu)**3-1.d0/Pmu)*sinpmu-3.d0/(pmu)**2.0d0*cospmu)*invpmu**2.0d0
            h_1=(sinpmu/Pmu**2-cospmu/Pmu)*invpmu
            h_sum=h_sum+h_saved(2)*h_1+h_saved(3)*h_2+h_saved(4)*h_3
        elseif(order==4) then
            h_4=(5.0d0*pmu*(-21.0d0 + 2.0d0*pmu**2.0d0)*cospmu + (105.0d0 - 45.0d0*pmu**2.0d0 &
&             + pmu**4)*sinpmu)/pmu**5*invpmu**4
            h_3 = (pmu*(-15.0d0 + pmu**2.0d0)*cospmu + 3.0d0*(5.0d0 - 2.0d0*pmu**2.0d0)*sinpmu)/pmu**4.0d0*invpmu**3.0d0
            h_2=((3.d0/(Pmu)**3-1.d0/Pmu)*sinpmu-3.d0/(pmu)**2*cospmu)*invpmu**2
            h_1=(sinpmu/Pmu**2-cospmu/Pmu)*invpmu
            h_sum=h_sum+h_saved(2)*h_1+h_saved(3)*h_2+h_saved(4)*h_3+h_saved(5)*h_4
        elseif (order==5) then
             h_5 = (-(pmu*(945.0d0 - 105.0d0*pmu**2.0d0 + pmu**4.0d0)*cospmu) + 15.0d0*(63.0d0 - 28.0d0*pmu**2.0d0 &
&             + pmu**4.0d0)*sinpmu)/pmu**6.0d0
            h_4=(5.0d0*pmu*(-21.0d0 + 2.0d0*pmu**2.0d0)*cospmu + (105.0d0 - 45.0d0*pmu**2.0d0 &
&             + pmu**4.0d0)*sinpmu)/pmu**5.0d0*invpmu**4.0d0
            h_3 = (pmu*(-15.0d0 + pmu**2.0d0)*cospmu + 3.0d0*(5.0d0 - 2.0d0*pmu**2.0d0)*sinpmu)/pmu**4.0d0*invpmu**3.0d0
            h_2=((3.d0/(Pmu)**3-1.d0/Pmu)*sinpmu-3.d0/(pmu)**2.0d0*cospmu)*invpmu**2.0d0
            h_1=(sinpmu/Pmu**2-cospmu/Pmu)*invpmu
            h_sum=h_sum+h_saved(2)*h_1+h_saved(3)*h_2+h_saved(4)*h_3+h_saved(5)*h_4+h_5*h_saved(6)
       elseif(order>=6) then
!            h_5 = (-(pmu*(945.0d0 - 105.0d0*pmu**2.0d0 + pmu**4.0d0)*cospmu) + 15.0d0*(63.0d0 - 28.0d0*pmu**2.0d0 &
!&             + pmu**4.0d0)*sinpmu)/pmu**6.0d0
!            h_4=(5.0d0*pmu*(-21.0d0 + 2.0d0*pmu**2.0d0)*cospmu + (105.0d0 - 45.0d0*pmu**2.0d0 &
!&             + pmu**4.0d0)*sinpmu)/pmu**5.0d0*invpmu**4.0d0
!            h_3 = (pmu*(-15.0d0 + pmu**2.0d0)*cospmu + 3.0d0*(5.0d0 - 2.0d0*pmu**2.0d0)*sinpmu)/pmu**4.0d0*invpmu**3.0d0
!            h_2=((3.d0/(Pmu)**3-1.d0/Pmu)*sinpmu-3.d0/(pmu)**2.0d0*cospmu)*invpmu**2.0d0
!            h_1=(sinpmu/Pmu**2-cospmu/Pmu)*invpmu
!            h_sum=h_sum+h_saved(2)*h_1+h_saved(3)*h_2+h_saved(4)*h_3
            !call dphrec(pmu,allbessels,100,order+1,size(mu))

            do beta=1,order
                ratan=0.0_dp

                !omp$ simd reduction(+:ratan)
                do i=0,40

                    ratan=ratan+(H**2.0_dp-1.0_dp)**dble(i)*BesselQ2(i+1,:)*BesselQ(i+beta+1,:)

                end do

                h_sum=h_sum+H**beta*ratan*invpmu**beta*h_saved(beta+1)
          enddo

        end if


!         do i=1,size(mu)
!                if (abs(pmu(i))<0.05) then
!                    call van(allbessels(0:6,i),6,pmu(i))
!                    allbessels(7:40,i)=0.0_dp
!                elseif(abs(pmu(i))>100) then
!                    print*,'oh my goood'
!                else
!                    call van(allbessels(0:40,i),20,pmu(i))
!                 end if
!
!         enddo
!
!       !  print*,(allbessels(0:18,1))
!       ! stop
!        h_sum=0.0_dp
!        do beta=0,order
!            hqr = (mu / h)**beta
!            bess= allbessels(beta,:)
!            h_sum=h_sum+bess*hqr*h_saved(beta+1)
!        enddo



     End Subroutine bessel0sum

end module Bessels_j0