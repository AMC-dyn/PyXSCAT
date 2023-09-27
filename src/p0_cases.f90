!This module contains the P0 cases for every legendre decomposition term
!j0,j1, and j2 and the table constructor

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
        REAL(KIND=dp), INTENT(IN), DIMENSION(:),allocatable         :: q


        ! loop variables
        INTEGER(KIND=ikind) :: i, j, k


        ! loop for all combinations of angular momenta
        do k = 0, lmax4
            do j = 0, lmax4
                do i = 0, lmax4
                    ! Fill the current LMN
                    if ((i+j+k)<16) then
                        CALL P_LMN(P0, i, j, k, q)
                    endif

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
        REAL(KIND=dp), INTENT(IN), DIMENSION(:),allocatable         :: q

        ! loop variables
        INTEGER(KIND=ikind) :: i, j, k

        p0=0.0_dp
        ! loop for all combinations of angular momenta
        do k = 0, lmax4
            do j = 0, lmax4
                do i = 0, lmax4

                    if ((i+j+k)<=16) then
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
        REAL(KIND=dp), INTENT(IN), DIMENSION(:),allocatable    :: q


        ! For ordering
        INTEGER(KIND=ikind), DIMENSION(1:3)     :: A
        ! helper variable
        !INTEGER(KIND=ikind) :: i

        ! Order L, M and N
        A(1) = L
        A(2) = M
        A(3) = N
        CALL Bubble_Sort(A)

        ! These are analytical solutions to the limits of the derivatives of j_0.
        ! They were generated in Wolfram Mathematica 13.1 by Mats Simmermacher in July 2023.

        if (mod(L,2)/=0 .or. mod(M,2)/=0 .or. mod(N,2)/=0) then
            P0(:,L+1,M+1,N+1)=0.0_dp

        elseif (sum(A)==0) then
            P0(:,L+1,M+1,N+1)=1.0_dp

        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==2) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/3.0_dp)*q**2.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/15.0_dp)*q**4.0_dp
        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/5.0_dp)*q**4.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==2) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/105.0_dp)*q**6.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/35.0_dp)*q**6.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/315.0_dp)*q**8.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. A(3)==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/105.0_dp)*q**8.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. A(3)==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1155.0_dp)*q**10.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. A(3)==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/5005.0_dp)*q**12.0_dp

        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/7.0_dp)*q**6.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/63.0_dp)*q**8.0_dp
        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/9.0_dp)*q**8.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/693.0_dp)*q**10.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/231.0_dp)*q**10.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/99.0_dp)*q**10.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/3003.0_dp)*q**12.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=(5.0_dp/3003.0_dp)*q**12.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/1287.0_dp)*q**12.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/429.0_dp)*q**12.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/15015.0_dp)*q**14.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/9009.0_dp)*q**14.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/6435.0_dp)*q**14.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1287.0_dp)*q**14.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/51051.0_dp)*q**16.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/36465.0_dp)*q**16.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/21879.0_dp)*q**16.0_dp
        elseif (A(1)==0 .and. A(2)==8 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/21879.0_dp)*q**16.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. A(3)==6) then
            P0(:,L+1,M+1,N+1)=(-5.0_dp/969969.0_dp)*q**18.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/138567.0_dp)*q**18.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/415701.0_dp)*q**18.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(5.0_dp/2909907.0_dp)*q**20.0_dp
        elseif (A(1)==4 .and. A(2)==8 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/415701.0_dp)*q**20.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(-5.0_dp/9561123.0_dp)*q**22.0_dp
        elseif (A(1)==8 .and. A(2)==8 .and. A(3)==8) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/47805615.0_dp)*q**24.0_dp

        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/11.0_dp)*q**10.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/143.0_dp)*q**12.0_dp
        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/13.0_dp)*q**12.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/2145.0_dp)*q**14.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/715.0_dp)*q**14.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/195.0_dp)*q**14.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/12155.0_dp)*q**16.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/2431.0_dp)*q**16.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/3315.0_dp)*q**16.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/1105.0_dp)*q**16.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(-3.0_dp/230945.0_dp)*q**18.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/46189.0_dp)*q**18.0_dp
        elseif (A(1)==0 .and. A(2)==8 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/46189.0_dp)*q**18.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/20995.0_dp)*q**18.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/4199.0_dp)*q**18.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/323323.0_dp)*q**20.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/138567.0_dp)*q**20.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/146965.0_dp)*q**20.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/88179.0_dp)*q**20.0_dp
        elseif (A(1)==0 .and. A(2)==10 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(3.0_dp/46189.0_dp)*q**20.0_dp
        elseif (A(1)==0 .and. A(2)==8 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/12597.0_dp)*q**20.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(-5.0_dp/7436429.0_dp)*q**22.0_dp
        elseif (A(1)==4 .and. A(2)==8 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1062347.0_dp)*q**22.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/676039.0_dp)*q**22.0_dp
        elseif (A(1)==2 .and. A(2)==10 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(-3.0_dp/1062347.0_dp)*q**22.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/289731.0_dp)*q**22.0_dp
        elseif (A(1)==0 .and. A(2)==10 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(-3.0_dp/96577.0_dp)*q**22.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/5311735.0_dp)*q**24.0_dp
        elseif (A(1)==4 .and. A(2)==10 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(9.0_dp/26558675.0_dp)*q**24.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/3380195.0_dp)*q**24.0_dp
        elseif (A(1)==4 .and. A(2)==8 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/2414425.0_dp)*q**24.0_dp
        elseif (A(1)==2 .and. A(2)==10 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(3.0_dp/2414425.0_dp)*q**24.0_dp
        elseif (A(1)==0 .and. A(2)==12 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(33.0_dp/2414425.0_dp)*q**24.0_dp
        elseif (A(1)==8 .and. A(2)==8 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/143416845.0_dp)*q**26.0_dp
        elseif (A(1)==6 .and. A(2)==10 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/15935205.0_dp)*q**26.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/13037895.0_dp)*q**26.0_dp
        elseif (A(1)==4 .and. A(2)==10 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/7243275.0_dp)*q**26.0_dp
        elseif (A(1)==2 .and. A(2)==12 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/21729825.0_dp)*q**26.0_dp
        elseif (A(1)==8 .and. A(2)==10 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/462120945.0_dp)*q**28.0_dp
        elseif (A(1)==8 .and. A(2)==8 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/378098955.0_dp)*q**28.0_dp
        elseif (A(1)==6 .and. A(2)==10 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/42010995.0_dp)*q**28.0_dp
        elseif (A(1)==4 .and. A(2)==12 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/210054975.0_dp)*q**28.0_dp
        elseif (A(1)==10 .and. A(2)==10 .and. A(3)==10) then
            P0(:,L+1,M+1,N+1)=(-21.0_dp/4775249765.0_dp)*q**30.0_dp
        elseif (A(1)==8 .and. A(2)==10 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/1302340845.0_dp)*q**30.0_dp
        elseif (A(1)==6 .and. A(2)==12 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/1302340845.0_dp)*q**30.0_dp
        elseif (A(1)==10 .and. A(2)==10 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/4775249765.0_dp)*q**32.0_dp
        elseif (A(1)==8 .and. A(2)==12 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/3907022535.0_dp)*q**32.0_dp
        elseif (A(1)==10 .and. A(2)==12 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/2170568075.0_dp)*q**34.0_dp
        elseif (A(1)==12 .and. A(2)==12 .and. A(3)==12) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/80311018775.0_dp)*q**36.0_dp

        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/15.0_dp)*q**14.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/255.0_dp)*q**16.0_dp
        elseif (A(1)==0 .and. A(2)==0 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/17.0_dp)*q**16.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/4845.0_dp)*q**18.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1615.0_dp)*q**18.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/323.0_dp)*q**18.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/33915.0_dp)*q**20.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/6783.0_dp)*q**20.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/6783.0_dp)*q**20.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/2261.0_dp)*q**20.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/260015.0_dp)*q**22.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/156009.0_dp)*q**22.0_dp
        elseif (A(1)==0 .and. A(2)==8 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/22287.0_dp)*q**22.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/52003.0_dp)*q**22.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-5.0_dp/52003.0_dp)*q**22.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/1300075.0_dp)*q**24.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/557175.0_dp)*q**24.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(3.0_dp/1300075.0_dp)*q**24.0_dp
        elseif (A(1)==0 .and. A(2)==10 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(3.0_dp/185725.0_dp)*q**24.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/260015.0_dp)*q**24.0_dp
        elseif (A(1)==0 .and. A(2)==8 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/37145.0_dp)*q**24.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/7020405.0_dp)*q**26.0_dp
        elseif (A(1)==4 .and. A(2)==8 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/5014575.0_dp)*q**26.0_dp
        elseif (A(1)==2 .and. A(2)==10 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1671525.0_dp)*q**26.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/2340135.0_dp)*q**26.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1002915.0_dp)*q**26.0_dp
        elseif (A(1)==0 .and. A(2)==12 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/1671525.0_dp)*q**26.0_dp
        elseif (A(1)==0 .and. A(2)==10 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/111435.0_dp)*q**26.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/29084535.0_dp)*q**28.0_dp
        elseif (A(1)==4 .and. A(2)==10 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/16158075.0_dp)*q**28.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/13572783.0_dp)*q**28.0_dp
        elseif (A(1)==4 .and. A(2)==8 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/9694845.0_dp)*q**28.0_dp
        elseif (A(1)==2 .and. A(2)==12 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/48474225.0_dp)*q**28.0_dp
        elseif (A(1)==2 .and. A(2)==10 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/3231615.0_dp)*q**28.0_dp
        elseif (A(1)==0 .and. A(2)==14 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(143.0_dp/48474225.0_dp)*q**28.0_dp
        elseif (A(1)==0 .and. A(2)==12 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/3231615.0_dp)*q**28.0_dp
        elseif (A(1)==8 .and. A(2)==8 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/901620585.0_dp)*q**30.0_dp
        elseif (A(1)==6 .and. A(2)==10 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/100180065.0_dp)*q**30.0_dp
        elseif (A(1)==4 .and. A(2)==12 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/500900325.0_dp)*q**30.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/60108039.0_dp)*q**30.0_dp
        elseif (A(1)==4 .and. A(2)==10 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/33393355.0_dp)*q**30.0_dp
        elseif (A(1)==2 .and. A(2)==14 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-143.0_dp/1502700975.0_dp)*q**30.0_dp
        elseif (A(1)==2 .and. A(2)==12 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/100180065.0_dp)*q**30.0_dp
        elseif (A(1)==0 .and. A(2)==14 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-143.0_dp/100180065.0_dp)*q**30.0_dp
        elseif (A(1)==8 .and. A(2)==10 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/3305942145.0_dp)*q**32.0_dp
        elseif (A(1)==6 .and. A(2)==12 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/300540195.0_dp)*q**32.0_dp
        elseif (A(1)==8 .and. A(2)==8 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/1983565287.0_dp)*q**32.0_dp
        elseif (A(1)==6 .and. A(2)==10 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/220396143.0_dp)*q**32.0_dp
        elseif (A(1)==4 .and. A(2)==14 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/1502700975.0_dp)*q**32.0_dp
        elseif (A(1)==4 .and. A(2)==12 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/100180065.0_dp)*q**32.0_dp
        elseif (A(1)==2 .and. A(2)==14 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/300540195.0_dp)*q**32.0_dp
        elseif (A(1)==0 .and. A(2)==16 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/20036013.0_dp)*q**32.0_dp
        elseif (A(1)==10 .and. A(2)==10 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1836634525.0_dp)*q**34.0_dp
        elseif (A(1)==8 .and. A(2)==12 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1502700975.0_dp)*q**34.0_dp
        elseif (A(1)==8 .and. A(2)==10 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1101980715.0_dp)*q**34.0_dp
        elseif (A(1)==6 .and. A(2)==14 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-13.0_dp/10518906825.0_dp)*q**34.0_dp
        elseif (A(1)==6 .and. A(2)==12 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/701260455.0_dp)*q**34.0_dp
        elseif (A(1)==4 .and. A(2)==14 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-13.0_dp/3506302275.0_dp)*q**34.0_dp
        elseif (A(1)==2 .and. A(2)==16 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-13.0_dp/701260455.0_dp)*q**34.0_dp
        elseif (A(1)==10 .and. A(2)==12 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/6177770675.0_dp)*q**36.0_dp
        elseif (A(1)==8 .and. A(2)==14 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/55599936075.0_dp)*q**36.0_dp
        elseif (A(1)==10 .and. A(2)==10 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(3.0_dp/13591095485.0_dp)*q**36.0_dp
        elseif (A(1)==8 .and. A(2)==12 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/3706662405.0_dp)*q**36.0_dp
        elseif (A(1)==6 .and. A(2)==14 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/25946636835.0_dp)*q**36.0_dp
        elseif (A(1)==4 .and. A(2)==16 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/8648878945.0_dp)*q**36.0_dp
        elseif (A(1)==12 .and. A(2)==12 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/240933056325.0_dp)*q**38.0_dp
        elseif (A(1)==10 .and. A(2)==14 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/18533312025.0_dp)*q**38.0_dp
        elseif (A(1)==10 .and. A(2)==12 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/16062203755.0_dp)*q**38.0_dp
        elseif (A(1)==8 .and. A(2)==14 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/11119987215.0_dp)*q**38.0_dp
        elseif (A(1)==6 .and. A(2)==16 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/5189327367.0_dp)*q**38.0_dp
        elseif (A(1)==12 .and. A(2)==14 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/759865793025.0_dp)*q**40.0_dp
        elseif (A(1)==12 .and. A(2)==12 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/658550353955.0_dp)*q**40.0_dp
        elseif (A(1)==10 .and. A(2)==14 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/50657719535.0_dp)*q**40.0_dp
        elseif (A(1)==8 .and. A(2)==16 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/30394631721.0_dp)*q**40.0_dp
        elseif (A(1)==14 .and. A(2)==14 .and. A(3)==14) then
            P0(:,L+1,M+1,N+1)=(-143.0_dp/32674229100075.0_dp)*q**42.0_dp
        elseif (A(1)==12 .and. A(2)==14 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/2178281940005.0_dp)*q**42.0_dp
        elseif (A(1)==10 .and. A(2)==16 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-3.0_dp/435656388001.0_dp)*q**42.0_dp
        elseif (A(1)==14 .and. A(2)==14 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(143.0_dp/98022687300225.0_dp)*q**44.0_dp
        elseif (A(1)==12 .and. A(2)==16 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/6534845820015.0_dp)*q**44.0_dp
        elseif (A(1)==14 .and. A(2)==16 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(-143.0_dp/307137753540705.0_dp)*q**46.0_dp
        elseif (A(1)==16 .and. A(2)==16 .and. A(3)==16) then
            P0(:,L+1,M+1,N+1)=(143.0_dp/1003316661566303.0_dp)*q**48.0_dp

        else
            print*, "Limits of derivatives of j_0 are only implemented up to g-functions.", L, M, N
            stop

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
        REAL(KIND=dp), INTENT(IN), DIMENSION(:),allocatable     :: q
        real(kind=dp):: pi

        ! For ordering
        INTEGER(KIND=ikind), DIMENSION(1:3)     :: A
        ! helper variable
        !INTEGER(KIND=ikind) :: i
        pi=dacos(-1.0_dp)
        P0=0.0_dp
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
        INTEGER(KIND=ikind), DIMENSION(1:2)     :: A
        ! helper variable
        !INTEGER(KIND=ikind) :: i
        pi=dacos(-1.0_dp)
        ! Order L, M and N
        A(1) = L
        A(2) = M
        CALL Bubble_Sort(A)

        ! These are analytical solutions to the limits of the derivatives of j_2 P_2.
        ! They were generated in Wolfram Mathematica 13.1 by Mats Simmermacher in July 2023.

        if (mod(L,2)/=0 .or. mod(M,2)/=0 .or. mod(N,2)/=0) then
            P0(:,L+1,M+1,N+1)=0.0_dp

        elseif (L+M==2*N) then
            P0(:,L+1,M+1,N+1)=0.0_dp

        elseif (A(1)==0 .and. A(2)==2 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/30.0_dp)*q**2.0_dp
        elseif (A(1)==0 .and. A(2)==0 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/15.0_dp)*q**2.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/105.0_dp)*q**4.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/210.0_dp)*q**4.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/35.0_dp)*q**4.0_dp
        elseif (A(1)==0 .and. A(2)==0 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/35.0_dp)*q**4.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/210.0_dp)*q**6.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/210.0_dp)*q**6.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/6930.0_dp)*q**8.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/3465.0_dp)*q**8.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/1155.0_dp)*q**8.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1155.0_dp)*q**8.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/15015.0_dp)*q**10.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/30030.0_dp)*q**10.0_dp

        elseif (A(1)==0 .and. A(2)==6 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/42.0_dp)*q**6.0_dp
        elseif (A(1)==0 .and. A(2)==0 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/21.0_dp)*q**6.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/693.0_dp)*q**8.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/1386.0_dp)*q**8.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-5.0_dp/1386.0_dp)*q**8.0_dp
        elseif (A(1)==0 .and. A(2)==8 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/99.0_dp)*q**8.0_dp
        elseif (A(1)==0 .and. A(2)==0 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-4.0_dp/99.0_dp)*q**8.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/9009.0_dp)*q**10.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/9009.0_dp)*q**10.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-5.0_dp/6006.0_dp)*q**10.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/6006.0_dp)*q**10.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/3003.0_dp)*q**10.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-5.0_dp/2574.0_dp)*q**10.0_dp
        elseif (A(1)==0 .and. A(2)==8 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1287.0_dp)*q**10.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/2574.0_dp)*q**10.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/30030.0_dp)*q**12.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/30030.0_dp)*q**12.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/3003.0_dp)*q**12.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/12870.0_dp)*q**12.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/6006.0_dp)*q**12.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/6435.0_dp)*q**12.0_dp
        elseif (A(1)==4 .and. A(2)==8 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/2145.0_dp)*q**12.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/2145.0_dp)*q**12.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/510510.0_dp)*q**14.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/255255.0_dp)*q**14.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/153153.0_dp)*q**14.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/153153.0_dp)*q**14.0_dp
        elseif (A(1)==4 .and. A(2)==8 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/109395.0_dp)*q**14.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/218790.0_dp)*q**14.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/43758.0_dp)*q**14.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/43758.0_dp)*q**14.0_dp
        elseif (A(1)==0 .and. A(2)==8 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/21879.0_dp)*q**14.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(5.0_dp/43758.0_dp)*q**14.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/969969.0_dp)*q**16.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1939938.0_dp)*q**16.0_dp
        elseif (A(1)==4 .and. A(2)==8 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/692835.0_dp)*q**16.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/692835.0_dp)*q**16.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(5.0_dp/831402.0_dp)*q**16.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/831402.0_dp)*q**16.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/415701.0_dp)*q**16.0_dp
        elseif (A(1)==8 .and. A(2)==8 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(28.0_dp/415701.0_dp)*q**16.0_dp
        elseif (A(1)==0 .and. A(2)==8 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-14.0_dp/415701.0_dp)*q**16.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1939938.0_dp)*q**18.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/1939938.0_dp)*q**18.0_dp
        elseif (A(1)==8 .and. A(2)==8 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/415701.0_dp)*q**18.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/831402.0_dp)*q**18.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(5.0_dp/133855722.0_dp)*q**20.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-5.0_dp/66927861.0_dp)*q**20.0_dp
        elseif (A(1)==8 .and. A(2)==8 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/9561123.0_dp)*q**20.0_dp
        elseif (A(1)==4 .and. A(2)==8 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/9561123.0_dp)*q**20.0_dp
        elseif (A(1)==8 .and. A(2)==8 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/47805615.0_dp)*q**22.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/95611230.0_dp)*q**22.0_dp

        elseif (A(1)==0 .and. A(2)==10 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-5.0_dp/286.0_dp)*q**10.0_dp
        elseif (A(1)==0 .and. A(2)==0 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(5.0_dp/143.0_dp)*q**10.0_dp
        elseif (A(1)==2 .and. A(2)==10 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/715.0_dp)*q**12.0_dp
        elseif (A(1)==0 .and. A(2)==10 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/1430.0_dp)*q**12.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-3.0_dp/1430.0_dp)*q**12.0_dp
        elseif (A(1)==0 .and. A(2)==12 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/65.0_dp)*q**12.0_dp
        elseif (A(1)==0 .and. A(2)==0 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/65.0_dp)*q**12.0_dp
        elseif (A(1)==2 .and. A(2)==10 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/36465.0_dp)*q**14.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(4.0_dp/36465.0_dp)*q**14.0_dp
        elseif (A(1)==4 .and. A(2)==10 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/24310.0_dp)*q**14.0_dp
        elseif (A(1)==0 .and. A(2)==10 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/24310.0_dp)*q**14.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(4.0_dp/12155.0_dp)*q**14.0_dp
        elseif (A(1)==2 .and. A(2)==12 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/6630.0_dp)*q**14.0_dp
        elseif (A(1)==0 .and. A(2)==12 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/3315.0_dp)*q**14.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/6630.0_dp)*q**14.0_dp
        elseif (A(1)==4 .and. A(2)==10 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/92378.0_dp)*q**16.0_dp
        elseif (A(1)==2 .and. A(2)==10 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/230945.0_dp)*q**16.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/461890.0_dp)*q**16.0_dp
        elseif (A(1)==6 .and. A(2)==10 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(4.0_dp/46189.0_dp)*q**16.0_dp
        elseif (A(1)==0 .and. A(2)==10 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/92378.0_dp)*q**16.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/92378.0_dp)*q**16.0_dp
        elseif (A(1)==2 .and. A(2)==12 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/25194.0_dp)*q**16.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/12597.0_dp)*q**16.0_dp
        elseif (A(1)==4 .and. A(2)==12 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(4.0_dp/20995.0_dp)*q**16.0_dp
        elseif (A(1)==0 .and. A(2)==12 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/20995.0_dp)*q**16.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/4199.0_dp)*q**16.0_dp
        elseif (A(1)==4 .and. A(2)==10 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-3.0_dp/3233230.0_dp)*q**18.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(3.0_dp/1616615.0_dp)*q**18.0_dp
        elseif (A(1)==6 .and. A(2)==10 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/323323.0_dp)*q**18.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/323323.0_dp)*q**18.0_dp
        elseif (A(1)==8 .and. A(2)==10 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-3.0_dp/92378.0_dp)*q**18.0_dp
        elseif (A(1)==4 .and. A(2)==12 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/146965.0_dp)*q**18.0_dp
        elseif (A(1)==2 .and. A(2)==12 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/293930.0_dp)*q**18.0_dp
        elseif (A(1)==0 .and. A(2)==10 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/92378.0_dp)*q**18.0_dp
        elseif (A(1)==0 .and. A(2)==8 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/46189.0_dp)*q**18.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(3.0_dp/293930.0_dp)*q**18.0_dp
        elseif (A(1)==6 .and. A(2)==12 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-3.0_dp/58786.0_dp)*q**18.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(3.0_dp/58786.0_dp)*q**18.0_dp
        elseif (A(1)==6 .and. A(2)==10 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/7436429.0_dp)*q**20.0_dp
        elseif (A(1)==4 .and. A(2)==10 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/14872858.0_dp)*q**20.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-5.0_dp/14872858.0_dp)*q**20.0_dp
        elseif (A(1)==8 .and. A(2)==10 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/6374082.0_dp)*q**20.0_dp
        elseif (A(1)==2 .and. A(2)==10 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/3187041.0_dp)*q**20.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-5.0_dp/6374082.0_dp)*q**20.0_dp
        elseif (A(1)==4 .and. A(2)==12 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/3380195.0_dp)*q**20.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-4.0_dp/3380195.0_dp)*q**20.0_dp
        elseif (A(1)==6 .and. A(2)==12 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/579462.0_dp)*q**20.0_dp
        elseif (A(1)==2 .and. A(2)==12 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/4056234.0_dp)*q**20.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-4.0_dp/2028117.0_dp)*q**20.0_dp
        elseif (A(1)==10 .and. A(2)==10 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(15.0_dp/1062347.0_dp)*q**20.0_dp
        elseif (A(1)==0 .and. A(2)==10 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-15.0_dp/2124694.0_dp)*q**20.0_dp
        elseif (A(1)==8 .and. A(2)==12 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(5.0_dp/289731.0_dp)*q**20.0_dp
        elseif (A(1)==0 .and. A(2)==12 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/289731.0_dp)*q**20.0_dp
        elseif (A(1)==0 .and. A(2)==8 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-4.0_dp/289731.0_dp)*q**20.0_dp
        elseif (A(1)==6 .and. A(2)==10 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/37182145.0_dp)*q**22.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/37182145.0_dp)*q**22.0_dp
        elseif (A(1)==8 .and. A(2)==10 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/10623470.0_dp)*q**22.0_dp
        elseif (A(1)==4 .and. A(2)==10 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/53117350.0_dp)*q**22.0_dp
        elseif (A(1)==4 .and. A(2)==8 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/26558675.0_dp)*q**22.0_dp
        elseif (A(1)==6 .and. A(2)==12 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/6760390.0_dp)*q**22.0_dp
        elseif (A(1)==4 .and. A(2)==12 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/16900975.0_dp)*q**22.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/4828850.0_dp)*q**22.0_dp
        elseif (A(1)==10 .and. A(2)==10 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-12.0_dp/26558675.0_dp)*q**22.0_dp
        elseif (A(1)==2 .and. A(2)==10 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(6.0_dp/26558675.0_dp)*q**22.0_dp
        elseif (A(1)==8 .and. A(2)==12 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-4.0_dp/7243275.0_dp)*q**22.0_dp
        elseif (A(1)==2 .and. A(2)==12 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/14486550.0_dp)*q**22.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/14486550.0_dp)*q**22.0_dp
        elseif (A(1)==10 .and. A(2)==12 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-33.0_dp/4828850.0_dp)*q**22.0_dp
        elseif (A(1)==0 .and. A(2)==12 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(6.0_dp/2414425.0_dp)*q**22.0_dp
        elseif (A(1)==0 .and. A(2)==10 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(21.0_dp/4828850.0_dp)*q**22.0_dp
        elseif (A(1)==8 .and. A(2)==10 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/95611230.0_dp)*q**24.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/95611230.0_dp)*q**24.0_dp
        elseif (A(1)==10 .and. A(2)==10 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/26558675.0_dp)*q**24.0_dp
        elseif (A(1)==6 .and. A(2)==12 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/60843510.0_dp)*q**24.0_dp
        elseif (A(1)==4 .and. A(2)==10 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/53117350.0_dp)*q**24.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/30421755.0_dp)*q**24.0_dp
        elseif (A(1)==8 .and. A(2)==12 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/21729825.0_dp)*q**24.0_dp
        elseif (A(1)==4 .and. A(2)==8 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/21729825.0_dp)*q**24.0_dp
        elseif (A(1)==10 .and. A(2)==12 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/4828850.0_dp)*q**24.0_dp
        elseif (A(1)==2 .and. A(2)==12 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/14486550.0_dp)*q**24.0_dp
        elseif (A(1)==2 .and. A(2)==10 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/7243275.0_dp)*q**24.0_dp
        elseif (A(1)==12 .and. A(2)==12 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(22.0_dp/7243275.0_dp)*q**24.0_dp
        elseif (A(1)==0 .and. A(2)==12 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/7243275.0_dp)*q**24.0_dp
        elseif (A(1)==8 .and. A(2)==10 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/8318177010.0_dp)*q**26.0_dp
        elseif (A(1)==8 .and. A(2)==8 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/4159088505.0_dp)*q**26.0_dp
        elseif (A(1)==10 .and. A(2)==10 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/462120945.0_dp)*q**26.0_dp
        elseif (A(1)==6 .and. A(2)==10 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/462120945.0_dp)*q**26.0_dp
        elseif (A(1)==8 .and. A(2)==12 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/378098955.0_dp)*q**26.0_dp
        elseif (A(1)==6 .and. A(2)==12 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/756197910.0_dp)*q**26.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/151239582.0_dp)*q**26.0_dp
        elseif (A(1)==10 .and. A(2)==12 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/420109950.0_dp)*q**26.0_dp
        elseif (A(1)==4 .and. A(2)==12 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/210054975.0_dp)*q**26.0_dp
        elseif (A(1)==4 .and. A(2)==10 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/84021990.0_dp)*q**26.0_dp
        elseif (A(1)==12 .and. A(2)==12 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/126032985.0_dp)*q**26.0_dp
        elseif (A(1)==2 .and. A(2)==12 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/252065970.0_dp)*q**26.0_dp
        elseif (A(1)==10 .and. A(2)==10 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/14325749295.0_dp)*q**28.0_dp
        elseif (A(1)==8 .and. A(2)==10 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/28651498590.0_dp)*q**28.0_dp
        elseif (A(1)==8 .and. A(2)==12 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/11721067605.0_dp)*q**28.0_dp
        elseif (A(1)==8 .and. A(2)==8 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-14.0_dp/11721067605.0_dp)*q**28.0_dp
        elseif (A(1)==10 .and. A(2)==12 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/520936338.0_dp)*q**28.0_dp
        elseif (A(1)==6 .and. A(2)==12 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/2604681690.0_dp)*q**28.0_dp
        elseif (A(1)==6 .and. A(2)==10 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/1302340845.0_dp)*q**28.0_dp
        elseif (A(1)==12 .and. A(2)==12 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(44.0_dp/6511704225.0_dp)*q**28.0_dp
        elseif (A(1)==4 .and. A(2)==12 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-22.0_dp/6511704225.0_dp)*q**28.0_dp
        elseif (A(1)==10 .and. A(2)==12 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/28651498590.0_dp)*q**30.0_dp
        elseif (A(1)==8 .and. A(2)==10 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/28651498590.0_dp)*q**30.0_dp
        elseif (A(1)==12 .and. A(2)==12 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1302340845.0_dp)*q**30.0_dp
        elseif (A(1)==6 .and. A(2)==12 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/2604681690.0_dp)*q**30.0_dp
        elseif (A(1)==10 .and. A(2)==12 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/47752497650.0_dp)*q**32.0_dp
        elseif (A(1)==10 .and. A(2)==10 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/23876248825.0_dp)*q**32.0_dp
        elseif (A(1)==12 .and. A(2)==12 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/19535112675.0_dp)*q**32.0_dp
        elseif (A(1)==8 .and. A(2)==12 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/19535112675.0_dp)*q**32.0_dp
        elseif (A(1)==12 .and. A(2)==12 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/80311018775.0_dp)*q**34.0_dp
        elseif (A(1)==10 .and. A(2)==12 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/160622037550.0_dp)*q**34.0_dp

        elseif (A(1)==0 .and. A(2)==14 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/510.0_dp)*q**14.0_dp
        elseif (A(1)==0 .and. A(2)==0 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/255.0_dp)*q**14.0_dp
        elseif (A(1)==2 .and. A(2)==14 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(4.0_dp/4845.0_dp)*q**16.0_dp
        elseif (A(1)==0 .and. A(2)==14 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/1938.0_dp)*q**16.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-13.0_dp/9690.0_dp)*q**16.0_dp
        elseif (A(1)==0 .and. A(2)==16 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(4.0_dp/323.0_dp)*q**16.0_dp
        elseif (A(1)==0 .and. A(2)==0 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-8.0_dp/323.0_dp)*q**16.0_dp
        elseif (A(1)==2 .and. A(2)==14 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/33915.0_dp)*q**18.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/33915.0_dp)*q**18.0_dp
        elseif (A(1)==4 .and. A(2)==14 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-3.0_dp/22610.0_dp)*q**18.0_dp
        elseif (A(1)==0 .and. A(2)==14 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/22610.0_dp)*q**18.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/11305.0_dp)*q**18.0_dp
        elseif (A(1)==2 .and. A(2)==16 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-3.0_dp/4522.0_dp)*q**18.0_dp
        elseif (A(1)==0 .and. A(2)==16 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/2261.0_dp)*q**18.0_dp
        elseif (A(1)==0 .and. A(2)==2 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(5.0_dp/4522.0_dp)*q**18.0_dp
        elseif (A(1)==4 .and. A(2)==14 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/222870.0_dp)*q**20.0_dp
        elseif (A(1)==2 .and. A(2)==14 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/780045.0_dp)*q**20.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/1560090.0_dp)*q**20.0_dp
        elseif (A(1)==6 .and. A(2)==14 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(5.0_dp/156009.0_dp)*q**20.0_dp
        elseif (A(1)==0 .and. A(2)==14 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/312018.0_dp)*q**20.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/312018.0_dp)*q**20.0_dp
        elseif (A(1)==2 .and. A(2)==16 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/44574.0_dp)*q**20.0_dp
        elseif (A(1)==2 .and. A(2)==2 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/22287.0_dp)*q**20.0_dp
        elseif (A(1)==4 .and. A(2)==16 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(5.0_dp/52003.0_dp)*q**20.0_dp
        elseif (A(1)==0 .and. A(2)==16 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/52003.0_dp)*q**20.0_dp
        elseif (A(1)==0 .and. A(2)==4 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/7429.0_dp)*q**20.0_dp
        elseif (A(1)==4 .and. A(2)==14 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/2600150.0_dp)*q**22.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/1300075.0_dp)*q**22.0_dp
        elseif (A(1)==6 .and. A(2)==14 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-4.0_dp/3900225.0_dp)*q**22.0_dp
        elseif (A(1)==2 .and. A(2)==14 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/3900225.0_dp)*q**22.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/780045.0_dp)*q**22.0_dp
        elseif (A(1)==8 .and. A(2)==14 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/1114350.0_dp)*q**22.0_dp
        elseif (A(1)==0 .and. A(2)==14 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/1114350.0_dp)*q**22.0_dp
        elseif (A(1)==0 .and. A(2)==8 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/111435.0_dp)*q**22.0_dp
        elseif (A(1)==4 .and. A(2)==16 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-4.0_dp/1300075.0_dp)*q**22.0_dp
        elseif (A(1)==2 .and. A(2)==16 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/520030.0_dp)*q**22.0_dp
        elseif (A(1)==2 .and. A(2)==4 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/2600150.0_dp)*q**22.0_dp
        elseif (A(1)==6 .and. A(2)==16 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/520030.0_dp)*q**22.0_dp
        elseif (A(1)==0 .and. A(2)==16 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/260015.0_dp)*q**22.0_dp
        elseif (A(1)==0 .and. A(2)==6 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/520030.0_dp)*q**22.0_dp
        elseif (A(1)==6 .and. A(2)==14 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/11700675.0_dp)*q**24.0_dp
        elseif (A(1)==4 .and. A(2)==14 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/23401350.0_dp)*q**24.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/7800450.0_dp)*q**24.0_dp
        elseif (A(1)==8 .and. A(2)==14 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/3343050.0_dp)*q**24.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/3343050.0_dp)*q**24.0_dp
        elseif (A(1)==4 .and. A(2)==16 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/3900225.0_dp)*q**24.0_dp
        elseif (A(1)==4 .and. A(2)==4 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/3900225.0_dp)*q**24.0_dp
        elseif (A(1)==10 .and. A(2)==14 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/557175.0_dp)*q**24.0_dp
        elseif (A(1)==6 .and. A(2)==16 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/1560090.0_dp)*q**24.0_dp
        elseif (A(1)==2 .and. A(2)==16 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/4680270.0_dp)*q**24.0_dp
        elseif (A(1)==0 .and. A(2)==14 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1114350.0_dp)*q**24.0_dp
        elseif (A(1)==0 .and. A(2)==10 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/371450.0_dp)*q**24.0_dp
        elseif (A(1)==2 .and. A(2)==6 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/2340135.0_dp)*q**24.0_dp
        elseif (A(1)==8 .and. A(2)==16 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/334305.0_dp)*q**24.0_dp
        elseif (A(1)==0 .and. A(2)==8 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/334305.0_dp)*q**24.0_dp
        elseif (A(1)==6 .and. A(2)==14 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/203591745.0_dp)*q**26.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(4.0_dp/203591745.0_dp)*q**26.0_dp
        elseif (A(1)==8 .and. A(2)==14 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/290845350.0_dp)*q**26.0_dp
        elseif (A(1)==4 .and. A(2)==14 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/290845350.0_dp)*q**26.0_dp
        elseif (A(1)==4 .and. A(2)==8 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(4.0_dp/145422675.0_dp)*q**26.0_dp
        elseif (A(1)==10 .and. A(2)==14 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/9694845.0_dp)*q**26.0_dp
        elseif (A(1)==2 .and. A(2)==14 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/48474225.0_dp)*q**26.0_dp
        elseif (A(1)==2 .and. A(2)==10 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(4.0_dp/48474225.0_dp)*q**26.0_dp
        elseif (A(1)==6 .and. A(2)==16 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/19389690.0_dp)*q**26.0_dp
        elseif (A(1)==4 .and. A(2)==16 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/67863915.0_dp)*q**26.0_dp
        elseif (A(1)==4 .and. A(2)==6 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/135727830.0_dp)*q**26.0_dp
        elseif (A(1)==8 .and. A(2)==16 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/5816907.0_dp)*q**26.0_dp
        elseif (A(1)==2 .and. A(2)==16 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/58169070.0_dp)*q**26.0_dp
        elseif (A(1)==2 .and. A(2)==8 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/58169070.0_dp)*q**26.0_dp
        elseif (A(1)==12 .and. A(2)==14 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-143.0_dp/96948450.0_dp)*q**26.0_dp
        elseif (A(1)==0 .and. A(2)==14 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/19389690.0_dp)*q**26.0_dp
        elseif (A(1)==0 .and. A(2)==12 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(44.0_dp/48474225.0_dp)*q**26.0_dp
        elseif (A(1)==10 .and. A(2)==16 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-13.0_dp/6463230.0_dp)*q**26.0_dp
        elseif (A(1)==0 .and. A(2)==16 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/3231615.0_dp)*q**26.0_dp
        elseif (A(1)==0 .and. A(2)==10 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/6463230.0_dp)*q**26.0_dp
        elseif (A(1)==8 .and. A(2)==14 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/360648234.0_dp)*q**28.0_dp
        elseif (A(1)==6 .and. A(2)==14 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/901620585.0_dp)*q**28.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/1803241170.0_dp)*q**28.0_dp
        elseif (A(1)==10 .and. A(2)==14 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(4.0_dp/500900325.0_dp)*q**28.0_dp
        elseif (A(1)==4 .and. A(2)==14 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1001800650.0_dp)*q**28.0_dp
        elseif (A(1)==4 .and. A(2)==10 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/1001800650.0_dp)*q**28.0_dp
        elseif (A(1)==6 .and. A(2)==16 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(5.0_dp/841512546.0_dp)*q**28.0_dp
        elseif (A(1)==6 .and. A(2)==6 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-5.0_dp/420756273.0_dp)*q**28.0_dp
        elseif (A(1)==8 .and. A(2)==16 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(4.0_dp/300540195.0_dp)*q**28.0_dp
        elseif (A(1)==4 .and. A(2)==16 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/300540195.0_dp)*q**28.0_dp
        elseif (A(1)==4 .and. A(2)==8 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/60108039.0_dp)*q**28.0_dp
        elseif (A(1)==12 .and. A(2)==14 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(121.0_dp/3005401950.0_dp)*q**28.0_dp
        elseif (A(1)==2 .and. A(2)==14 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-22.0_dp/1502700975.0_dp)*q**28.0_dp
        elseif (A(1)==2 .and. A(2)==12 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-77.0_dp/3005401950.0_dp)*q**28.0_dp
        elseif (A(1)==10 .and. A(2)==16 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/200360130.0_dp)*q**28.0_dp
        elseif (A(1)==2 .and. A(2)==16 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/200360130.0_dp)*q**28.0_dp
        elseif (A(1)==2 .and. A(2)==10 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/20036013.0_dp)*q**28.0_dp
        elseif (A(1)==14 .and. A(2)==14 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(1001.0_dp/1502700975.0_dp)*q**28.0_dp
        elseif (A(1)==0 .and. A(2)==14 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-1001.0_dp/3005401950.0_dp)*q**28.0_dp
        elseif (A(1)==12 .and. A(2)==16 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(77.0_dp/100180065.0_dp)*q**28.0_dp
        elseif (A(1)==0 .and. A(2)==16 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-22.0_dp/100180065.0_dp)*q**28.0_dp
        elseif (A(1)==0 .and. A(2)==12 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/20036013.0_dp)*q**28.0_dp
        elseif (A(1)==8 .and. A(2)==14 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/19835652870.0_dp)*q**30.0_dp
        elseif (A(1)==8 .and. A(2)==8 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/9917826435.0_dp)*q**30.0_dp
        elseif (A(1)==10 .and. A(2)==14 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1101980715.0_dp)*q**30.0_dp
        elseif (A(1)==6 .and. A(2)==10 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/1101980715.0_dp)*q**30.0_dp
        elseif (A(1)==12 .and. A(2)==14 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/333933550.0_dp)*q**30.0_dp
        elseif (A(1)==8 .and. A(2)==16 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/661188429.0_dp)*q**30.0_dp
        elseif (A(1)==6 .and. A(2)==16 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1322376858.0_dp)*q**30.0_dp
        elseif (A(1)==4 .and. A(2)==14 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/1001800650.0_dp)*q**30.0_dp
        elseif (A(1)==4 .and. A(2)==12 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/500900325.0_dp)*q**30.0_dp
        elseif (A(1)==6 .and. A(2)==8 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/440792286.0_dp)*q**30.0_dp
        elseif (A(1)==10 .and. A(2)==16 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-3.0_dp/734653810.0_dp)*q**30.0_dp
        elseif (A(1)==4 .and. A(2)==10 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(3.0_dp/734653810.0_dp)*q**30.0_dp
        elseif (A(1)==14 .and. A(2)==14 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-26.0_dp/1502700975.0_dp)*q**30.0_dp
        elseif (A(1)==2 .and. A(2)==14 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/1502700975.0_dp)*q**30.0_dp
        elseif (A(1)==12 .and. A(2)==16 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/100180065.0_dp)*q**30.0_dp
        elseif (A(1)==2 .and. A(2)==16 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/200360130.0_dp)*q**30.0_dp
        elseif (A(1)==2 .and. A(2)==12 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/66786710.0_dp)*q**30.0_dp
        elseif (A(1)==14 .and. A(2)==16 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(-13.0_dp/40072026.0_dp)*q**30.0_dp
        elseif (A(1)==0 .and. A(2)==16 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/100180065.0_dp)*q**30.0_dp
        elseif (A(1)==0 .and. A(2)==14 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/66786710.0_dp)*q**30.0_dp
        elseif (A(1)==10 .and. A(2)==14 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/16529710725.0_dp)*q**32.0_dp
        elseif (A(1)==8 .and. A(2)==14 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/33059421450.0_dp)*q**32.0_dp
        elseif (A(1)==8 .and. A(2)==10 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/6611884290.0_dp)*q**32.0_dp
        elseif (A(1)==12 .and. A(2)==14 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/3005401950.0_dp)*q**32.0_dp
        elseif (A(1)==6 .and. A(2)==14 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/10518906825.0_dp)*q**32.0_dp
        elseif (A(1)==6 .and. A(2)==12 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/4207562730.0_dp)*q**32.0_dp
        elseif (A(1)==8 .and. A(2)==16 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/9917826435.0_dp)*q**32.0_dp
        elseif (A(1)==8 .and. A(2)==8 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-4.0_dp/9917826435.0_dp)*q**32.0_dp
        elseif (A(1)==10 .and. A(2)==16 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/2203961430.0_dp)*q**32.0_dp
        elseif (A(1)==6 .and. A(2)==16 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/15427730010.0_dp)*q**32.0_dp
        elseif (A(1)==6 .and. A(2)==10 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-4.0_dp/7713865005.0_dp)*q**32.0_dp
        elseif (A(1)==14 .and. A(2)==14 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/10518906825.0_dp)*q**32.0_dp
        elseif (A(1)==4 .and. A(2)==14 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-13.0_dp/21037813650.0_dp)*q**32.0_dp
        elseif (A(1)==12 .and. A(2)==16 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/701260455.0_dp)*q**32.0_dp
        elseif (A(1)==4 .and. A(2)==16 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/3506302275.0_dp)*q**32.0_dp
        elseif (A(1)==4 .and. A(2)==12 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-4.0_dp/3506302275.0_dp)*q**32.0_dp
        elseif (A(1)==14 .and. A(2)==16 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(169.0_dp/21037813650.0_dp)*q**32.0_dp
        elseif (A(1)==2 .and. A(2)==16 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-13.0_dp/4207562730.0_dp)*q**32.0_dp
        elseif (A(1)==2 .and. A(2)==14 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-52.0_dp/10518906825.0_dp)*q**32.0_dp
        elseif (A(1)==16 .and. A(2)==16 .and. N==0) then
            P0(:,L+1,M+1,N+1)=(104.0_dp/701260455.0_dp)*q**32.0_dp
        elseif (A(1)==0 .and. A(2)==16 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-52.0_dp/701260455.0_dp)*q**32.0_dp
        elseif (A(1)==10 .and. A(2)==14 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/67955477425.0_dp)*q**34.0_dp
        elseif (A(1)==10 .and. A(2)==10 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/67955477425.0_dp)*q**34.0_dp
        elseif (A(1)==12 .and. A(2)==14 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/22239974430.0_dp)*q**34.0_dp
        elseif (A(1)==8 .and. A(2)==14 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/111199872150.0_dp)*q**34.0_dp
        elseif (A(1)==8 .and. A(2)==12 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/55599936075.0_dp)*q**34.0_dp
        elseif (A(1)==10 .and. A(2)==16 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/16309314582.0_dp)*q**34.0_dp
        elseif (A(1)==8 .and. A(2)==16 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/40773286455.0_dp)*q**34.0_dp
        elseif (A(1)==8 .and. A(2)==10 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(7.0_dp/81546572910.0_dp)*q**34.0_dp
        elseif (A(1)==14 .and. A(2)==14 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-52.0_dp/389199552525.0_dp)*q**34.0_dp
        elseif (A(1)==6 .and. A(2)==14 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(26.0_dp/389199552525.0_dp)*q**34.0_dp
        elseif (A(1)==12 .and. A(2)==16 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-4.0_dp/25946636835.0_dp)*q**34.0_dp
        elseif (A(1)==6 .and. A(2)==16 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/51893273670.0_dp)*q**34.0_dp
        elseif (A(1)==6 .and. A(2)==12 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/7413324810.0_dp)*q**34.0_dp
        elseif (A(1)==14 .and. A(2)==16 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(-143.0_dp/259466368350.0_dp)*q**34.0_dp
        elseif (A(1)==4 .and. A(2)==16 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(26.0_dp/129733184175.0_dp)*q**34.0_dp
        elseif (A(1)==4 .and. A(2)==14 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/37066624050.0_dp)*q**34.0_dp
        elseif (A(1)==16 .and. A(2)==16 .and. N==2) then
            P0(:,L+1,M+1,N+1)=(-13.0_dp/3706662405.0_dp)*q**34.0_dp
        elseif (A(1)==2 .and. A(2)==16 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(13.0_dp/7413324810.0_dp)*q**34.0_dp
        elseif (A(1)==12 .and. A(2)==14 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/160622037550.0_dp)*q**36.0_dp
        elseif (A(1)==10 .and. A(2)==12 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/160622037550.0_dp)*q**36.0_dp
        elseif (A(1)==14 .and. A(2)==14 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/55599936075.0_dp)*q**36.0_dp
        elseif (A(1)==10 .and. A(2)==16 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(3.0_dp/353368482610.0_dp)*q**36.0_dp
        elseif (A(1)==8 .and. A(2)==14 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/111199872150.0_dp)*q**36.0_dp
        elseif (A(1)==10 .and. A(2)==10 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-3.0_dp/176684241305.0_dp)*q**36.0_dp
        elseif (A(1)==12 .and. A(2)==16 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/48186611265.0_dp)*q**36.0_dp
        elseif (A(1)==8 .and. A(2)==12 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/48186611265.0_dp)*q**36.0_dp
        elseif (A(1)==14 .and. A(2)==16 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/17297757890.0_dp)*q**36.0_dp
        elseif (A(1)==6 .and. A(2)==16 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/51893273670.0_dp)*q**36.0_dp
        elseif (A(1)==6 .and. A(2)==14 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/25946636835.0_dp)*q**36.0_dp
        elseif (A(1)==16 .and. A(2)==16 .and. N==4) then
            P0(:,L+1,M+1,N+1)=(2.0_dp/8648878945.0_dp)*q**36.0_dp
        elseif (A(1)==4 .and. A(2)==16 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/8648878945.0_dp)*q**36.0_dp
        elseif (A(1)==12 .and. A(2)==14 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/19756510618650.0_dp)*q**38.0_dp
        elseif (A(1)==12 .and. A(2)==12 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/9878255309325.0_dp)*q**38.0_dp
        elseif (A(1)==14 .and. A(2)==14 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/759865793025.0_dp)*q**38.0_dp
        elseif (A(1)==10 .and. A(2)==14 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/759865793025.0_dp)*q**38.0_dp
        elseif (A(1)==12 .and. A(2)==16 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/658550353955.0_dp)*q**38.0_dp
        elseif (A(1)==10 .and. A(2)==16 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/1317100707910.0_dp)*q**38.0_dp
        elseif (A(1)==10 .and. A(2)==12 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/263420141582.0_dp)*q**38.0_dp
        elseif (A(1)==14 .and. A(2)==16 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(-7.0_dp/911838951630.0_dp)*q**38.0_dp
        elseif (A(1)==8 .and. A(2)==16 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/455919475815.0_dp)*q**38.0_dp
        elseif (A(1)==8 .and. A(2)==14 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/182367790326.0_dp)*q**38.0_dp
        elseif (A(1)==16 .and. A(2)==16 .and. N==6) then
            P0(:,L+1,M+1,N+1)=(-5.0_dp/212762422047.0_dp)*q**38.0_dp
        elseif (A(1)==6 .and. A(2)==16 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(5.0_dp/425524844094.0_dp)*q**38.0_dp
        elseif (A(1)==14 .and. A(2)==14 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/32674229100075.0_dp)*q**40.0_dp
        elseif (A(1)==12 .and. A(2)==14 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/65348458200150.0_dp)*q**40.0_dp
        elseif (A(1)==12 .and. A(2)==16 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/28317665220065.0_dp)*q**40.0_dp
        elseif (A(1)==12 .and. A(2)==12 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-22.0_dp/28317665220065.0_dp)*q**40.0_dp
        elseif (A(1)==14 .and. A(2)==16 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/871312776002.0_dp)*q**40.0_dp
        elseif (A(1)==10 .and. A(2)==16 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/4356563880010.0_dp)*q**40.0_dp
        elseif (A(1)==10 .and. A(2)==14 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/2178281940005.0_dp)*q**40.0_dp
        elseif (A(1)==16 .and. A(2)==16 .and. N==8) then
            P0(:,L+1,M+1,N+1)=(4.0_dp/1306969164003.0_dp)*q**40.0_dp
        elseif (A(1)==8 .and. A(2)==16 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-2.0_dp/1306969164003.0_dp)*q**40.0_dp
        elseif (A(1)==14 .and. A(2)==16 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/65348458200150.0_dp)*q**42.0_dp
        elseif (A(1)==12 .and. A(2)==14 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(11.0_dp/65348458200150.0_dp)*q**42.0_dp
        elseif (A(1)==16 .and. A(2)==16 .and. N==10) then
            P0(:,L+1,M+1,N+1)=(-1.0_dp/2178281940005.0_dp)*q**42.0_dp
        elseif (A(1)==10 .and. A(2)==16 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(1.0_dp/4356563880010.0_dp)*q**42.0_dp
        elseif (A(1)==14 .and. A(2)==16 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(143.0_dp/9214132606221150.0_dp)*q**44.0_dp
        elseif (A(1)==14 .and. A(2)==14 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-143.0_dp/4607066303110575.0_dp)*q**44.0_dp
        elseif (A(1)==16 .and. A(2)==16 .and. N==12) then
            P0(:,L+1,M+1,N+1)=(22.0_dp/307137753540705.0_dp)*q**44.0_dp
        elseif (A(1)==12 .and. A(2)==16 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(-11.0_dp/307137753540705.0_dp)*q**44.0_dp
        elseif (A(1)==16 .and. A(2)==16 .and. N==14) then
            P0(:,L+1,M+1,N+1)=(-143.0_dp/15049749923494545.0_dp)*q**46.0_dp
        elseif (A(1)==14 .and. A(2)==16 .and. N==16) then
            P0(:,L+1,M+1,N+1)=(143.0_dp/30099499846989090.0_dp)*q**46.0_dp

        else
            print*, "Limits of derivatives of j_2 P_2 are only implemented up to g-functions.", L, M, N
            stop

        endif

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