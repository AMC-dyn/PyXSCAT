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

        real(kind=dp), intent(out), dimension(napos,napos,maxl*2+1,maxl+1,maxl+1) :: ddx,ddy,ddz


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
                            ddx(i, j, ls, l(ii)+1, l(jj)+1) = dx(ii, jj, ls, l(ii)+1, l(jj)+1)
                        end do

                        do ms=1,m(ii) + m(jj) + 1
                            ddy(i, j, ms, m(ii)+1, m(jj)+1) = dy(ii, jj, ms, m(ii)+1, m(jj)+1)
                        end do

                        do ns=1, n(ii) + n(jj) + 1
                            ddz(i, j, ns, n(ii)+1, n(jj)+1) = dz(ii, jj, ns, n(ii)+1, n(jj)+1)
                        end do

                     end do
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
                    ! do this in a loop: for L=1:l1+l2+1
                    !D(i,j,1,5,5)= MD(1,5, 5, PA, PB, gaP)
                    !D(i,j,2,5,5)= MD(2,5, 5, PA, PB, gaP)
                    print*, "case not programmed yet: l1/2= " , l1, l2
                    stop
            end if
        end do

    END SUBROUTINE fill_md_table



    SUBROUTINE integral_ijkr_pzero(nq,lmax1,lmax2,lmax3,lmax4,p0mat,dx,dy,dz,i,j,k,r,z1,z2,apos,cutoffz,cutoffmd,itgr)

        use types
        implicit none
        ! definition of input

        INTEGER(kind=ikind), INTENT(IN)                       :: nq, lmax1, lmax2, lmax3, lmax4, i, j, k, r
        REAL(kind=dp), INTENT(IN)                             :: cutoffz, cutoffmd
        REAL(kind=dp), INTENT(IN), DIMENSION(:,:,:,:)               :: p0mat
        REAL(kind=dp), INTENT(IN), DIMENSION(:,:,:,:,:)               :: dx, dy, dz
        REAL(kind=dp), INTENT(IN), DIMENSION(:,:,:)               :: z1, z2
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:)         :: apos
        ! definition of output
        REAL(kind=dp), INTENT(OUT), DIMENSION(nq)             :: itgr
        ! definition of loop indices
        INTEGER(kind=ikind)                                   :: l1, m1, l2, m2, l3, m3, l4, m4, h1
        INTEGER(kind=ikind)                                   :: l, m, n, lp, mp, np
        ! definition of internal variables
        INTEGER(kind=ikind)                                   :: n1, n2, n3, n4, posi, posj, posk, posr
        REAL(kind=dp)                                         :: ztot
        REAL(kind=dp)                                         :: mdl, mdm, mdn, mdlp, mdmp, prodd
        REAL(kind=dp), DIMENSION(size(Z1(:,1,1)))                          :: zij1, zkr1
        REAL(kind=dp), DIMENSION(size(Z2(:,1,1)))                          :: zij2, zkr2
        REAL(kind=dp), DIMENSION(nq)                          :: f


        posI=apos(i)

        itgr=0.0
        ! loop through all possible ways to get total angular momentum lmax1
        do l1 = 0, lmax1
            do m1 = 0, lmax1-l1
                n1 = lmax1-l1-m1
                posj = apos(j)
                ! loop through all possible ways to get total angular momentum lmax2
                do l2 = 0, lmax2
                    do m2 = 0, lmax2-l2
                        n2 = lmax2-l2-m2
                        zij1 = z1(:,posi,posj)
                        zij2 = z2(:,posi,posj)
                        posk = apos(k)
                        ! loop through all possible ways to get total angular momentum lmax3
                        do l3 = 0, lmax3
                            do m3 = 0, lmax3-l3
                                n3 = lmax3-l3-m3
                                posr = apos(r)
                                ! loop through all possible ways to get total angular momentum lmax4
                                do l4 = 0, lmax4
                                    do m4 = 0, lmax4-l4
                                        n4 = lmax4-l4-m4

                                        zkr1 = z1(:,posk,posr)
                                        zkr2 = z2(:,posk,posr)

                                        ! total prefactor
                                        ztot = sum(zij1*zkr2 + zij2*zkr1) / 8

                                        ! continue only if larger
                                        if (abs(ztot) < cutoffz) then
                                            posr = posr+1
                                            cycle
                                        end if
                                        ! the 6-dimensional sum over MD coefficents
                                        do l = 0, l1+l2
                                            mdl = dx(i,j,l+1,l1+1,l2+1) * ztot
                                            if (mdl == 0) cycle
                                            do m = 0, m1+m2
                                                mdm = dy(i,j,m+1,m1+1,m2+1) * mdl
                                                if (mdm == 0) cycle
                                                do n = 0, n1+n2
                                                    h1 = (-1)**(l+m+n)
                                                    mdn = dz(i,j,n+1,n1+1,n2+1) * mdm * h1
                                                    if (mdn == 0) cycle
                                                    do lp = 0, l3+l4
                                                        mdlp = dx(k,r,lp+1,l3+1,l4+1) * mdn
                                                        if (mdlp == 0) cycle
                                                        do mp = 0, m3+m4
                                                            mdmp = dy(k,r,mp+1,m3+1,m4+1) * mdlp
                                                            if (mdmp == 0) cycle
                                                            do np = 0, n3+n4
                                                                prodd = dz(k,r,np+1,n3+1,n4+1) * mdmp
                                                                ! cutoff after md
                                                                if (abs(prodd) < cutoffmd) cycle

                                                                ! add the contribution to the total

                                                                f = p0mat(:,l+lp+1,m+mp+1,n+np+1)
                                                                itgr = itgr + prodd*f
                                                            end do
                                                        end do
                                                    end do
                                                end do
                                            end do
                                        end do
                                        posr = posr+1
                                    end do
                                end do
                                posk = posk+1
                            end do
                        end do
                        posj = posj+1
                    end do
                end do
                posi = posi+1
            end do
        end do

    END SUBROUTINE

    subroutine tot_integral_k_ijkr(mu,lmax1,lmax2,lmax3,lmax4,hx,hy,hz,h,dx, dy, dz, i,j, k, r, z, z2, apos, cutoffz, &
                cutoffmd,int_res)

        use types
        implicit none


        !INTEGER, PARAMETER :: dp = selected_real_kind(2*precision(1.0_dp))
        integer(kind=selected_int_kind(8)), intent(in)  :: lmax1,lmax2,lmax3,lmax4,i,j,k,r
        integer(kind=selected_int_kind(8)), dimension(:), intent(in) :: apos
        real(kind=dp), intent(in)              :: cutoffz,cutoffmd, hx, hy, hz, h
        real(kind=dp), intent(in), dimension(:,:,:) :: z, z2
        real(kind=dp), intent(in), dimension(:)   ::  mu

        real(kind=dp), intent(in),  dimension(:,:,:,:,:) :: dx,dy,dz



        integer(kind=selected_int_kind(8)) :: llmax, l, m, n, ka, posi, posj, posk, posr, ra
        integer(kind=selected_int_kind(8)) ::  l1, l2, l3, l4, m1, m2, m3, m4, n1, n2, n3, n4, h1
        integer(kind=selected_int_kind(8)) ::  ll2, mm2, nn2, lp, mp, np

        real(kind=dp)  ::coeff,prodd,ztot,mdn, mdl, mdm, mdlp, mdmp
        real(kind=dp), dimension(lmax1+lmax2+lmax3+lmax4+1)  :: bd
        real(kind=dp), dimension(:,:), allocatable :: a, b, c
        real(kind=dp), dimension(:,:,:,:), allocatable :: h_pre2
        real(kind=dp), dimension(:), allocatable :: h_saved, pmu, h_sum, h_0, h_1, h_r, muoh,zij, zij2, zkr, zkr2

        real(kind=dp), intent(out), dimension(size(mu)) :: int_res





        llmax=lmax1+lmax2+lmax3+lmax4

        allocate(zij(size(Z(:,1,1))), zij2(size(Z(:,1,1))), zkr(size(Z(:,1,1))), zkr2(size(Z(:,1,1))))
        zij=0.0_dp
        zij2=0.0_dp
        zkr=0.0_dp
        zkr2=0.0_dp


        allocate(a(llmax+1,llmax+1),b(llmax+1,llmax+1), c(llmax+1,llmax+1))
        a=0.0_dp
        b=0.0_dp
        c=0.0_dp
        a(1,1)=1.0_dp
        b(1,1)=1.0_dp
        c(1,1)=1.0_dp

        if (llmax>0) then
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
        h_pre2=0.0_dp
        bd=0.0_dp
        do l=0,llmax
            do m=0,llmax-l
                do n=0,(llmax-l-m)
                    call besselderiv(bd,l,m,n,a,b,c,llmax)
                    h_pre2(:,l+1,m+1,n+1) = bd
                enddo
            enddo
        enddo
        posi=apos(i)


        h_saved=0.0_dp
! loop through all possible ways to get total angular momentum lmax1
        do l1=0,lmax1
            do m1=0,(lmax1-l1)
                n1=lmax1-l1-m1
                posj=apos(j)
        !loop through all possible ways to get total angular momentum lmax2
                do l2=0,lmax2
                    do m2=0,(lmax2-l2)
                        n2=lmax2-l2-m2


                        zij(:)=z(:,posi,posj)
                        zij2(:)=z2(:,posi,posj)

                        posk=apos(k)
                        do l3=0,lmax3
                            do m3=0,(lmax3-l3)
                                n3=lmax3-l3-m3
                                posr=apos(r)

                        ! loop through all possible ways to get total angular momentum lmax4
                                do l4=0,lmax4
                                    do m4=0,(lmax4-l4)
                                        n4=lmax4-l4-m4

                                        zkr=z(:,posk,posr)
                                        zkr2=z2(:,posk,posr)
                                ! total prefactor

                                        ztot=sum(zij*zkr2+zij2*zkr)/8.0_dp

                                        if (abs(ztot)<cutoffz) then 
                                            posr=posr+1
                                            cycle
                                        endif


                                        do l=0,(l1+l2)
                                            mdl=dx(i,j,l+1,l1+1,l2+1)*ztot
                                   
                                            if (mdl==0.0_dp) cycle
                                                 
                                         
                                            do m=0,(m1+m2)
                                                mdm=dy(i,j,m+1,m1+1,m2+1)*mdl
                                                if (mdm==0.0_dp) cycle
                                         
                                                do n=0,(n1+n2)
                                                    h1=(-1)**(l+m+n)
                                                    mdn=dz(i,j,n+1,n1+1,n2+1)*mdm*h1
                                                    if (mdn==0.0_dp) cycle
                                                    
                                                    do lp=0,(l3+l4)
                                                        mdlp=dx(k,r,lp+1,l3+1,l4+1)*mdn
                                                        if (mdlp==0.0_dp) cycle

                                                    
                                                        ll2=l+lp+1
                                                        do mp=0,(m3+m4)
                                                            mdmp=dy(k,r,mp+1,m3+1,m4+1)*mdlp
                                                            if (mdmp==0.0_dp) cycle
                                                     
                                                            mm2=m+mp+1
                                                            do np=0,(n3+n4)
                                                                prodd=dz(k,r,np+1,n3+1,n4+1)*mdmp
                                                                if (abs(prodd)<cutoffmd) cycle

                                                                nn2=n+np+1

                                                                h_saved=h_saved+h_pre2(:,ll2,mm2,nn2)*prodd

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


        h_0=dsin(pmu)/pmu
        coeff=h_saved(1)
        h_sum=h_0*coeff

        if (llmax==1) then
            h_1=(dsin(pmu)/pmu**2.0_dp-dcos(pmu)/pmu)*mu/h
            coeff=h_saved(2)
            h_sum=h_sum+coeff*h_1
        elseif (llmax>1) then
            muoh=mu/h
            h_1=(dsin(pmu)/pmu**2.0_dp-dcos(pmu)/pmu)*muoh
            coeff=h_saved(2)
            h_sum=h_sum+coeff*h_1
            do ra=2,llmax
                coeff=h_saved(ra+1)
                h_r= ((2.0d0*ra-1.0d0)/(pmu)*h_1-h_0*muoh)*muoh
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
        use types
        implicit none

        !INTEGER, PARAMETER :: dp = kind(1.0d0)
        real(kind=dp), intent(out), dimension(llmax+1)  :: bd
        integer(kind=selected_int_kind(8)), intent(in)                 :: ll, mm, nn, llmax
        real(kind=dp), intent(in), dimension(llmax+1,llmax+1)  :: a, b, c

        ! loop and temp variables
        integer(kind=selected_int_kind(8)) :: ii, jj, kk, horder, temp, ceil
        real(kind=dp)       :: c1, c2, c3, ct2, ct3
        ! set this to 0 initially and accumulate
        bd = 0.0_dp
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


subroutine total_scattering_calculation(maxl, ipos,nipos,apos,napos,ga,l,m,n,xx,yy,zz, &
        mmod,q,nq,list1,listN1,list2,listN2,p0matrix, &
        cutoffz,cutoffmd,cutoffcentre,confs,civecs,ndiff,ep2,result)


    use types
    implicit none



        integer(kind=ikind), intent(in):: nipos, napos,  nq, maxl
        integer(kind=ikind), intent(in),dimension(:) :: l, m, n, apos, ipos
        integer(kind=ikind), intent(in),dimension(:) :: listN1,listN2
        integer(kind=ikind), intent(in),dimension(:,:) :: list1, list2, ep2, ndiff
        integer(kind=ikind), dimension(:,:), intent(in):: confs
        real(kind=dp), intent(in),dimension(:) :: ga, xx, yy, zz, q
        real(kind=dp), intent(in),dimension(:,:) :: mmod, civecs
        real(kind=dp), intent(in), dimension(:,:,:,:) :: p0matrix
        real(kind=dp), intent(in) :: cutoffmd, cutoffz,cutoffcentre

        real(kind=dp), intent(out), dimension(nq):: result

         real(kind=dp),  dimension(napos,napos,maxl*2+1,maxl+1,maxl+1) :: ddx,ddy,ddz
        real(kind=dp), dimension(napos,napos) :: px,py,pz
        real(kind=dp),  dimension(:,:,:), allocatable :: z1, z2
        real(kind=dp),  dimension(:), allocatable :: total,newtotal
        real(kind=dp),  dimension(nq,nipos,nipos) :: e12
        integer(kind=ikind), dimension(napos) :: ll
        integer(kind=ikind), dimension(:), allocatable :: m1, m2, m3, m4
        integer(kind=ikind), dimension(:,:), allocatable :: mat
        integer(kind=ikind):: nmat

        print*,confs(1,:)
        call createtwordm(confs,civecs,ndiff,ep2,mat,total)
        call reduce_density(mat,total,m1,m2,m3,m4,newtotal)

        allocate(z1(size(m1(:)), nipos, nipos), z2(size(m1(:)), nipos, nipos))

        nmat=size(m1(:))

        print*,nmat, size(newtotal), newtotal(1)
        call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,e12,ll,maxl, ipos,nipos,apos,napos,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q,nq,list1,listN1,list2,listN2)

        print*,'In between variable and integration'

        call integration(napos,px,py,pz,ll,p0matrix,ddx,ddy,ddz,z1,z2,apos,cutoffz,cutoffmd, cutoffcentre,q,e12,result)


        end subroutine total_scattering_calculation

        end module integrals_ijkr