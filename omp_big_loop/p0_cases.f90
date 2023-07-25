Module p0_cases

implicit none


contains

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

                    if ((i+j+k)<16) then
                        ! Fill the current LMN
                        CALL P_LMN_j2(P0, k, j, i, q)

                    end if

                end do
            end do
        end do

        P0=P0*sqrt(5.0_dp/dacos(-1.0_dp))
    END SUBROUTINE set_P0_j2

    SUBROUTINE set_P0_j1(P0, lmax4, q)

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

                    if ((i+j+k)<16) then
                        ! Fill the current LMN
                        CALL P_LMN_j1(P0, k, j, i, q)

                    end if

                end do
            end do
        end do

        P0=P0*sqrt(3.0_dp/dacos(-1.0_dp))
    END SUBROUTINE set_P0_j1

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

    SUBROUTINE P_LMN_j1(P0, L, M, N, q)

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
      if (mod(L+M+N,2)==0 .or. mod(N,2)==0) then
    P0(:,L+1,M+1,N+1)=0.0_dp
elseif (mod(L,2)/=0 .and. mod(M,2)/=0 .and. mod(N,2)/=0) then
    P0(:,L+1,M+1,N+1)=0.0_dp
elseif (L==0 .and. M==0 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/6.0_dp)*q
elseif (L==2 .and. M==0 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/30.0_dp)*q**3.0_dp
elseif (L==0 .and. M==2 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/30.0_dp)*q**3.0_dp
elseif (L==0 .and. M==0 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/10.0_dp)*q**3.0_dp
elseif (L==2 .and. M==2 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/210.0_dp)*q**5.0_dp
elseif (L==2 .and. M==0 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/70.0_dp)*q**5.0_dp
elseif (L==0 .and. M==2 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/70.0_dp)*q**5.0_dp
elseif (L==4 .and. M==0 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/70.0_dp)*q**5.0_dp
elseif (L==0 .and. M==4 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/70.0_dp)*q**5.0_dp
elseif (L==0 .and. M==0 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/14.0_dp)*q**5.0_dp
elseif (L==2 .and. M==2 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/630.0_dp)*q**7.0_dp
elseif (L==4 .and. M==2 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/630.0_dp)*q**7.0_dp
elseif (L==2 .and. M==4 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/630.0_dp)*q**7.0_dp
elseif (L==4 .and. M==0 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/210.0_dp)*q**7.0_dp
elseif (L==0 .and. M==4 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/210.0_dp)*q**7.0_dp
elseif (L==2 .and. M==0 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/126.0_dp)*q**7.0_dp
elseif (L==0 .and. M==2 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/126.0_dp)*q**7.0_dp
elseif (L==6 .and. M==0 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/126.0_dp)*q**7.0_dp
elseif (L==0 .and. M==6 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/126.0_dp)*q**7.0_dp
elseif (L==0 .and. M==0 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/18.0_dp)*q**7.0_dp
elseif (L==4 .and. M==2 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/2310.0_dp)*q**9.0_dp
elseif (L==2 .and. M==4 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/2310.0_dp)*q**9.0_dp
elseif (L==4 .and. M==4 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/2310.0_dp)*q**9.0_dp
elseif (L==2 .and. M==2 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/1386.0_dp)*q**9.0_dp
elseif (L==6 .and. M==2 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/1386.0_dp)*q**9.0_dp
elseif (L==2 .and. M==6 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/1386.0_dp)*q**9.0_dp
elseif (L==4 .and. M==0 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/462.0_dp)*q**9.0_dp
elseif (L==0 .and. M==4 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/462.0_dp)*q**9.0_dp
elseif (L==6 .and. M==0 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/462.0_dp)*q**9.0_dp
elseif (L==0 .and. M==6 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/462.0_dp)*q**9.0_dp
elseif (L==2 .and. M==0 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/198.0_dp)*q**9.0_dp
elseif (L==0 .and. M==2 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/198.0_dp)*q**9.0_dp
elseif (L==8 .and. M==0 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/198.0_dp)*q**9.0_dp
elseif (L==0 .and. M==8 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/198.0_dp)*q**9.0_dp
elseif (L==0 .and. M==0 .and. N==9) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/22.0_dp)*q**9.0_dp
elseif (L==4 .and. M==4 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/10010.0_dp)*q**11.0_dp
elseif (L==4 .and. M==2 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/6006.0_dp)*q**11.0_dp
elseif (L==2 .and. M==4 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/6006.0_dp)*q**11.0_dp
elseif (L==6 .and. M==2 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/6006.0_dp)*q**11.0_dp
elseif (L==2 .and. M==6 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/6006.0_dp)*q**11.0_dp
elseif (L==6 .and. M==4 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/6006.0_dp)*q**11.0_dp
elseif (L==4 .and. M==6 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/6006.0_dp)*q**11.0_dp
elseif (L==2 .and. M==2 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/2574.0_dp)*q**11.0_dp
elseif (L==6 .and. M==0 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(-5.0_dp/6006.0_dp)*q**11.0_dp
elseif (L==0 .and. M==6 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(-5.0_dp/6006.0_dp)*q**11.0_dp
elseif (L==4 .and. M==0 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/858.0_dp)*q**11.0_dp
elseif (L==0 .and. M==4 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/858.0_dp)*q**11.0_dp
elseif (L==8 .and. M==2 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/2574.0_dp)*q**11.0_dp
elseif (L==2 .and. M==8 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/2574.0_dp)*q**11.0_dp
elseif (L==8 .and. M==0 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/858.0_dp)*q**11.0_dp
elseif (L==0 .and. M==8 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/858.0_dp)*q**11.0_dp
elseif (L==2 .and. M==0 .and. N==9) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/286.0_dp)*q**11.0_dp
elseif (L==0 .and. M==2 .and. N==9) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/286.0_dp)*q**11.0_dp
elseif (L==10 .and. M==0 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/286.0_dp)*q**11.0_dp
elseif (L==0 .and. M==10 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/286.0_dp)*q**11.0_dp
elseif (L==0 .and. M==0 .and. N==11) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/26.0_dp)*q**11.0_dp
elseif (L==4 .and. M==4 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/30030.0_dp)*q**13.0_dp
elseif (L==6 .and. M==4 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/30030.0_dp)*q**13.0_dp
elseif (L==4 .and. M==6 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/30030.0_dp)*q**13.0_dp
elseif (L==6 .and. M==2 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/18018.0_dp)*q**13.0_dp
elseif (L==2 .and. M==6 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/18018.0_dp)*q**13.0_dp
elseif (L==4 .and. M==2 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/12870.0_dp)*q**13.0_dp
elseif (L==2 .and. M==4 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/12870.0_dp)*q**13.0_dp
elseif (L==6 .and. M==6 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/18018.0_dp)*q**13.0_dp
elseif (L==8 .and. M==2 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/12870.0_dp)*q**13.0_dp
elseif (L==2 .and. M==8 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/12870.0_dp)*q**13.0_dp
elseif (L==8 .and. M==4 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/12870.0_dp)*q**13.0_dp
elseif (L==4 .and. M==8 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/12870.0_dp)*q**13.0_dp
elseif (L==6 .and. M==0 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/2574.0_dp)*q**13.0_dp
elseif (L==0 .and. M==6 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/2574.0_dp)*q**13.0_dp
elseif (L==8 .and. M==0 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/2574.0_dp)*q**13.0_dp
elseif (L==0 .and. M==8 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/2574.0_dp)*q**13.0_dp
elseif (L==2 .and. M==2 .and. N==9) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/4290.0_dp)*q**13.0_dp
elseif (L==4 .and. M==0 .and. N==9) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/1430.0_dp)*q**13.0_dp
elseif (L==0 .and. M==4 .and. N==9) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/1430.0_dp)*q**13.0_dp
elseif (L==10 .and. M==2 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/4290.0_dp)*q**13.0_dp
elseif (L==2 .and. M==10 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/4290.0_dp)*q**13.0_dp
elseif (L==10 .and. M==0 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/1430.0_dp)*q**13.0_dp
elseif (L==0 .and. M==10 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/1430.0_dp)*q**13.0_dp
elseif (L==2 .and. M==0 .and. N==11) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/390.0_dp)*q**13.0_dp
elseif (L==0 .and. M==2 .and. N==11) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/390.0_dp)*q**13.0_dp
elseif (L==12 .and. M==0 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/390.0_dp)*q**13.0_dp
elseif (L==0 .and. M==12 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/390.0_dp)*q**13.0_dp
elseif (L==0 .and. M==0 .and. N==13) then
    P0(:,L+1,M+1,N+1)=(1.0_dp/30.0_dp)*q**13.0_dp
elseif (L==6 .and. M==4 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/102102.0_dp)*q**15.0_dp
elseif (L==4 .and. M==6 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/102102.0_dp)*q**15.0_dp
elseif (L==6 .and. M==6 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/102102.0_dp)*q**15.0_dp
elseif (L==4 .and. M==4 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/72930.0_dp)*q**15.0_dp
elseif (L==8 .and. M==4 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/72930.0_dp)*q**15.0_dp
elseif (L==4 .and. M==8 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/72930.0_dp)*q**15.0_dp
elseif (L==6 .and. M==2 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/43758.0_dp)*q**15.0_dp
elseif (L==2 .and. M==6 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/43758.0_dp)*q**15.0_dp
elseif (L==8 .and. M==2 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/43758.0_dp)*q**15.0_dp
elseif (L==2 .and. M==8 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/43758.0_dp)*q**15.0_dp
elseif (L==8 .and. M==6 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/43758.0_dp)*q**15.0_dp
elseif (L==6 .and. M==8 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/43758.0_dp)*q**15.0_dp
elseif (L==4 .and. M==2 .and. N==9) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/24310.0_dp)*q**15.0_dp
elseif (L==2 .and. M==4 .and. N==9) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/24310.0_dp)*q**15.0_dp
elseif (L==10 .and. M==2 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/24310.0_dp)*q**15.0_dp
elseif (L==2 .and. M==10 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/24310.0_dp)*q**15.0_dp
elseif (L==8 .and. M==0 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/43758.0_dp)*q**15.0_dp
elseif (L==0 .and. M==8 .and. N==7) then
    P0(:,L+1,M+1,N+1)=(-7.0_dp/43758.0_dp)*q**15.0_dp
elseif (L==10 .and. M==4 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/24310.0_dp)*q**15.0_dp
elseif (L==4 .and. M==10 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/24310.0_dp)*q**15.0_dp
elseif (L==6 .and. M==0 .and. N==9) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/4862.0_dp)*q**15.0_dp
elseif (L==0 .and. M==6 .and. N==9) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/4862.0_dp)*q**15.0_dp
elseif (L==10 .and. M==0 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/4862.0_dp)*q**15.0_dp
elseif (L==0 .and. M==10 .and. N==5) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/4862.0_dp)*q**15.0_dp
elseif (L==2 .and. M==2 .and. N==11) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/6630.0_dp)*q**15.0_dp
elseif (L==4 .and. M==0 .and. N==11) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/2210.0_dp)*q**15.0_dp
elseif (L==0 .and. M==4 .and. N==11) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/2210.0_dp)*q**15.0_dp
elseif (L==12 .and. M==2 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/6630.0_dp)*q**15.0_dp
elseif (L==2 .and. M==12 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/6630.0_dp)*q**15.0_dp
elseif (L==12 .and. M==0 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/2210.0_dp)*q**15.0_dp
elseif (L==0 .and. M==12 .and. N==3) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/2210.0_dp)*q**15.0_dp
elseif (L==2 .and. M==0 .and. N==13) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/510.0_dp)*q**15.0_dp
elseif (L==0 .and. M==2 .and. N==13) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/510.0_dp)*q**15.0_dp
elseif (L==14 .and. M==0 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/510.0_dp)*q**15.0_dp
elseif (L==0 .and. M==14 .and. N==1) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/510.0_dp)*q**15.0_dp
elseif (L==0 .and. M==0 .and. N==15) then
    P0(:,L+1,M+1,N+1)=(-1.0_dp/34.0_dp)*q**15.0_dp
        else

            print*, "CASE NOT PROGRAMED YET", L, M, N
            P0(:,L+1,M+1,N+1) = 0.0_dp
            !STOP

        end if
        !P0=P0*sqrt(5.0_dp/pi)
    END SUBROUTINE P_LMN_j1

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





End Module p0_cases