!-----------------------------------------------------------------------
! Fortran Modules for PyXSCAT code
! Andres Moreno Carrascosa and Mats Simmermacher, 2021
!-----------------------------------------------------------------------
module types
    implicit none
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
end module types

module integrals_ijkr

    use twordmreader


    implicit none 

    contains



     SUBROUTINE P_LMN(P0, L, M, N, q)
         use types

        ! the three values of the angular momentum
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
        call Bubble_Sort(A)

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

    subroutine integration(ncap,px,py,pz,ll,p0matrix,dx,dy,dz,z1,z2,apos,cutoffz,cutoffmd, cutoffcentre,q,e12,tsi)

        use types

        implicit none

        
        INTEGER(kind=ikind), INTENT(IN) :: ncap
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:) :: ll,apos
        
        REAL(kind=dp), intent(in), dimension(:,:,:,:,:) :: dx,dy,dz
        REAL(kind=dp), intent(in), dimension(:,:,:,:) :: p0matrix
        REAL(kind=dp), intent(in), dimension(:,:,:) :: z1,z2,e12
        REAL(kind=dp), intent(in), dimension(:,:) :: px,py,pz
        REAL(kind=dp), intent(in), dimension(:) :: q
        REAL(kind=dp), intent(in) :: cutoffz, cutoffmd,cutoffcentre
        
        
        REAL(kind=dp), dimension(size(q)) :: f
        REAL(kind=dp), intent(out), dimension(size(q)) :: tsi
        real(kind=dp) :: hx,hy,hz,h
        integer(kind=ikind) :: nq,i,j,k,r
        
        
        
        
        nq= size(q)
        tsi=0.0_dp

        !First big loop 
        do i=1,ncap
            do j=i+1,ncap
                do k=i+1,ncap
                    do r=k+1,ncap
                        hx = px(k, r) - px(i, j)
                        hy = py(k, r) - py(i, j)
                        hz = pz(k, r) - pz(i, j)
                        h = sqrt((hx * hx + hy * hy + hz * hz))
                        if (h < cutoffcentre) then

                            call integral_ijkr_pzero(nq, ll(i), ll(j), ll(k), ll(r), p0matrix, dx, dy, &
                                    dz, i, j, k, r, &
                                    z1, z2, apos, cutoffz, cutoffmd, f)

                        else

                            call tot_integral_k_ijkr(q, ll(i), ll(j), ll(k), ll(r), hx, hy, hz, h, dx, &
                                    dy, dz, &
                                    i, j, k, r, &
                                    z1, z2, apos, cutoffz, cutoffmd, f)



                        end if
                        tsi = tsi + 8.000 * f * e12(:, i, j) * e12(:, k, r)
                    end do
                end do
            end do
        end do

        do i=1,ncap
            do j=i+1,ncap
                do r=i+1,ncap
                    hx = px(i, r) - px(i, j)
                    hy = py(i, r) - py(i, j)
                    hz = pz(i, r) - pz(i, j)
                    h = sqrt((hx * hx + hy * hy + hz * hz))
                    if (h < cutoffcentre) then
                        call integral_ijkr_pzero(nq, ll(i), ll(j), ll(i), ll(r), p0matrix, dx, dy, &
                                dz, i, j, i, r, &
                                z1, z2, apos, cutoffz, cutoffmd, f)
                    else

                        call tot_integral_k_ijkr(q, ll(i), ll(j), ll(i), ll(r), hx, hy, hz, h, dx, &
                                dy, dz, &
                                i, j, i, r, &
                                z1, z2, apos, cutoffz, cutoffmd, f)



                    end if
                    tsi = tsi + 4.000 * f * e12(:, i, j) * e12(:, i, r)
                end do
            end do

        end do

        do i=1,ncap
            do k=1,ncap
                do r=k+1,ncap
                    hx = px(k, r) - px(i, i)
                    hy = py(k, r) - py(i, i)
                    hz = pz(k, r) - pz(i, i)
                    h = sqrt((hx * hx + hy * hy + hz * hz))
                    if (h < cutoffcentre) then
                        call integral_ijkr_pzero(nq, ll(i), ll(i), ll(k), ll(r), p0matrix, dx, dy, &
                                dz, i, i, k, r, &
                                z1, z2, apos, cutoffz, cutoffmd, f)
                    else

                        call tot_integral_k_ijkr(q, ll(i), ll(i), ll(k), ll(r), hx, hy, hz, h, dx, &
                                dy, dz, &
                                i, i, k, r, &
                                z1, z2, apos, cutoffz, cutoffmd, f)



                    end if
                    tsi = tsi+ 4.000 * f * e12(:, i, i) * e12(:, k, r)
                end do
            end do
        end do


        do i=1,ncap
            do k=i+1,ncap

                hx = px(k, k) - px(i, i)
                hy = py(k, k) - py(i, i)
                hz = pz(k, k) - pz(i, i)
                h = sqrt((hx * hx + hy * hy + hz * hz))
                if (h < cutoffcentre) then
                    call integral_ijkr_pzero(nq, ll(i), ll(i), ll(k), ll(k), p0matrix, dx, dy, &
                            dz, i, i, k, k, &
                            z1, z2, apos, cutoffz, cutoffmd, f)
                else

                    call tot_integral_k_ijkr(q, ll(i), ll(i), ll(k), ll(k), hx, hy, hz, h, dx, &
                            dy, dz, &
                            i, i, k, k, &
                            z1, z2, apos, cutoffz, cutoffmd, f)
                end if
                tsi = tsi+ 2.000 * f * e12(:, i, i) * e12(:, k, k)
            end do
        end do
        do i=1,ncap
            call integral_ijkr_pzero(nq, ll(i), ll(i), ll(i), ll(i), p0matrix, dx, dy, dz, i, i, i, i, &
            z1, z2, apos, cutoffz, cutoffmd,f)

            tsi = tsi + f * e12(:, i, i) * e12(:, i, i)

        end do

    print*, maxval(tsi)

    end subroutine integration

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



subroutine variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,e12,ll2,maxl, ipos,nipos,apos,napos,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q,nq,list1,listN1,list2,listN2)

    use types
        implicit none



        integer(kind=ikind), intent(in):: nipos, napos, nmat, nq, maxl
        integer(kind=ikind), intent(in),dimension(:) :: l, m, n, apos, ipos, m1, m2, m3, m4
        integer(kind=ikind), intent(in),dimension(:) :: listN1,listN2
        integer(kind=ikind), intent(in),dimension(:,:) :: list1, list2
        integer(kind=ikind), intent(out) :: ll2(napos)


        real(kind=dp), intent(in),dimension(:) :: ga, xx, yy, zz, total,q
        real(kind=dp), dimension(napos):: ga2, xx2, yy2, zz2
        real(kind=dp), intent(in),dimension(:,:) :: mmod
        real(kind=dp), intent(out), dimension(napos,napos):: px,py,pz
        real(kind=dp), intent(out), dimension(nmat,nipos,nipos) :: z1, z2
        real(kind=dp), intent(out), dimension(nq,nipos,nipos) :: e12

        real(kind=dp), intent(out), dimension(maxl*2+1,maxl+1,maxl+1,napos,napos) :: ddx,ddy,ddz


        integer(kind=ikind) :: max4l, max2l, i, j, k, N1, N2, nn1, nn2,ii,jj,ls, ms, ns
        integer(kind=ikind), dimension(nipos) :: ll
        real(kind=dp) :: pi, gap

        integer(kind=ikind), dimension(:), allocatable :: iduplicates,jduplicates
        real(kind=dp),  dimension(nipos,nipos,maxl*2+1,maxl+1,maxl+1):: dx,dy,dz
        real(kind=dp),  dimension(nq):: preexp
        real(kind=dp),  dimension(nmat):: temp1,temp2

        pi = acos(-1.0000)
        max4l=maxval(l)*4
        max2l=maxval(l)*2+1


        N1=size(ipos)


        write(*,*)'P0CALCULATED'
        print*,shape(z1), nmat
        z1=0.0_dp
        z2=0.0_dp
        print*,z1(1,1,1),m1(1)


        temp1=0.0
        temp2=0.0
        do  i=1,Nipos
            if (allocated(iduplicates)) deallocate(iduplicates)

            allocate(iduplicates(listN1(i)))
            iduplicates = list1(i,1:listN1(i))

            do j=1,nipos
                 if (allocated(jduplicates)) deallocate(jduplicates)

                 allocate(jduplicates(listN1(j)))
                 jduplicates = list1(j,1:listN1(j))

                 do nn1=1,size(iduplicates)
                     ii=iduplicates(nn1)
                     do nn2=1,size(jduplicates)
                        jj=jduplicates(nn2)

                        temp1 = total * (mmod(m1, ii) * mmod(m2, jj) + mmod(m1, jj) * mmod(m2, ii))
                        temp2 = mmod(m3, ii) * mmod(m4, jj) + mmod(m3, jj) * mmod(m4, ii)

                        z1(:, i, j) =  z1(:, i, j) + temp1
                        z2(:, i, j) =  z2(:, i, j) + temp2

                     enddo
                 enddo
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




        call fill_md_table(dx,l,xx,ga)
        call fill_md_table(dy,m,yy,ga)
        call fill_md_table(dz,n,zz,ga)


       ! allocate( ga2(size(apos)), xx2(size(apos)), yy2(size(apos)), zz2(size(apos)) )

        ll = l + m + n
        ga2 = ga(apos)
        xx2 = xx(apos)
        yy2 = yy(apos)
        zz2 = zz(apos)
        ll2 = ll(apos)

        N2=size(apos)

       ! allocate(ddx(n2,n2,max2l,maxl,maxl),ddy(n2,n2,max2l,maxl,maxl),ddz(n2,n2,max2l,maxl,maxl))
        !allocate(preexp(nq))

        do i=1,napos
            if (allocated(iduplicates)) deallocate(iduplicates)


            allocate(iduplicates(listN2(i)))
            iduplicates = list2(i,1:listN2(i))
            do j=1,napos
                 if (allocated(jduplicates)) deallocate(jduplicates)


                 allocate(jduplicates(listN2(j)))
                 jduplicates = list2(j,1:listN2(j))


                 gap = ga2(i) + ga2(j)
                 px(i,j)=(ga2(i)*xx2(i)+ga2(j)*xx2(j))/(gap)
                 py(i,j)=(ga2(i)*yy2(i)+ga2(j)*yy2(j))/(gap)

                 pz(i,j)=(ga2(i)*zz2(i)+ga2(j)*zz2(j))/(gap)
                 preexp= (pi/(gap))**1.5 * dexp(-q*q*0.25*(1/gap))


                 e12(:,i,j)=preexp*dexp(-(ga2(i)*ga2(j)/gap*((xx2(i)-xx2(j))**2+(yy2(i)-yy2(j))**2+(zz2(i)-zz2(j))**2)))

                 do nn1=1,size(iduplicates)
                     ii=iduplicates(nn1)
                     do nn2=1,size(jduplicates)
                        jj=jduplicates(nn2)

                        do ls=1,l(ii) + l(jj) + 1
                            ddx(ls, l(jj)+1, l(ii)+1, j,i) = dx(ii, jj, ls, l(ii)+1, l(jj)+1)
                        end do

                        do ms=1,m(ii) + m(jj) + 1
                            ddy(ms, m(jj)+1, m(ii)+1, j,i) = dy(ii, jj, ms, m(ii)+1, m(jj)+1)
                        end do

                        do ns=1, n(ii) + n(jj) + 1
                            ddz(ns, n(jj)+1,n(ii)+1, j,i) = dz(ii, jj, ns, n(ii)+1, n(jj)+1)
                        end do

                     end do
                 end do
            end do
        end do

    print*, 'leaving variables', size(ddx(1,1,1,:,1))
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
                    ! do this in a loop: for L=1:l1+l2+1
                    !D(i,j,1,5,5)= MD(1,5, 5, PA, PB, gaP)
                    !D(i,j,2,5,5)= MD(2,5, 5, PA, PB, gaP)
                    print*, "case not programmed yet: l1/2= " , l1, l2
                    stop
            end if
        end do

    END SUBROUTINE fill_md_table







subroutine total_scattering_calculation(maxl, ipos,nipos,apos,napos,ga,l,m,n,xx,yy,zz, &
        mmod,q,nq,list1,listN1,list2,listN2,p0matrix, &
        cutoffz,cutoffmd,cutoffcentre,confs,civecs,result)


    use types
    implicit none



        integer(kind=ikind), intent(in):: nipos, napos,  nq, maxl
        integer(kind=ikind), intent(in),dimension(:) :: l, m, n, apos, ipos
        integer(kind=ikind), intent(in),dimension(:) :: listN1,listN2
        integer(kind=ikind), intent(in),dimension(:,:) :: list1, list2
        integer(kind=ikind), dimension(:,:), intent(in):: confs
        real(kind=dp), intent(in),dimension(:) :: ga, xx, yy, zz, q
        real(kind=dp), intent(in),dimension(:,:) :: mmod, civecs
        real(kind=dp), intent(in), dimension(:,:,:,:) :: p0matrix
        real(kind=dp), intent(in) :: cutoffmd, cutoffz,cutoffcentre

        real(kind=dp), intent(out), dimension(nq):: result

         real(kind=dp),  dimension(maxl*2+1,maxl+1,maxl+1,napos,napos) :: ddx,ddy,ddz
        real(kind=dp), dimension(napos,napos) :: px,py,pz
        real(kind=dp),  dimension(:,:,:), allocatable :: z1, z2
        real(kind=dp),  dimension(:), allocatable :: total,newtotal
        real(kind=dp),  dimension(nq,nipos,nipos) :: e12
        integer(kind=ikind), dimension(napos) :: ll
        integer(kind=ikind), dimension(:), allocatable :: m1, m2, m3, m4
        integer(kind=ikind), dimension(:,:), allocatable :: mat,ep3,ndiff2
        integer(kind=ikind):: nmat,i,j
        real(kind=dp) :: start,time1,time2,time3,time4


        call cpu_time(start)
        call maxcoincidence(confs,ep3,ndiff2)
        call cpu_time(time1)
       ! allocate(ep3(size(confs(:,1)),size(confs(:,1))),ndiff2(size(confs(:,1)),size(confs(:,1))) )
        print*,'Time maxcoincidence',time1-start


        call createtwordm(confs,civecs,ndiff2,ep3,mat,total)

        call cpu_time(time2)
         print*,'Time maxcoincidence',time2-time1


        allocate(m1(size(mat(:,1))), m2(size(mat(:,1))), m3(size(mat(:,1))), m4(size(mat(:,1))))
!
!
!
        m1 = mat(:,1)
        m2 = mat(:,2)
        m3 = mat(:,3)
        m4 = mat(:,4)
        !call reduce_density(mat,total,m1,m2,m3,m4,newtotal)

        print*,'Reduced matrix'
        allocate(z1(size(m1), nipos, nipos), z2(size(m1), nipos, nipos))

        nmat=size(m1)


        call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,e12,ll,maxl, ipos,nipos,apos,napos,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q,nq,list1,listN1,list2,listN2)

        call cpu_time(time3)
        print*,'Time variables', time3-time2
        print*,shape(ddx),shape(ddy),shape(ddz)
        call integration(napos,px,py,pz,ll,p0matrix,ddx,ddy,ddz,z1,z2,apos,cutoffz,cutoffmd, cutoffcentre,q,e12,result)

        call cpu_time(time4)

         print*,'Time calc', time4-time3
        end subroutine total_scattering_calculation

        end module integrals_ijkr