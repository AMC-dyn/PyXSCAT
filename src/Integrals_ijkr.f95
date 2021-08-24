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
        real(kind=dp),dimension(:,:,:), allocatable :: dx1red,dy1red,dz1red,dx2red,dy2red,dz2red
        REAL(kind=dp), intent(in), dimension(:,:,:,:) :: p0matrix
        REAL(kind=dp), intent(in), dimension(:,:,:) ::e12
        real(kind=dp), intent(in), dimension(:,:,:)::z1,z2
        real(kind=dp), dimension(:,:,:,:), allocatable::zcontrred,zcontrred2
        REAL(kind=dp), intent(in), dimension(:,:) :: px,py,pz
        REAL(kind=dp), intent(in), dimension(:) :: q
        REAL(kind=dp), intent(in) :: cutoffz, cutoffmd,cutoffcentre
        
        
        REAL(kind=dp), dimension(size(q)) :: f
        integer(kind=ikind), dimension(:,:), allocatable :: posits
        real(kind=dp),dimension(:,:), allocatable :: za,zb,cmat
        integer(kind=ikind),dimension(:), allocatable ::posi,posj,posk,posr
        REAL(kind=dp), intent(out), dimension(size(q)) :: tsi
        real(kind=dp) :: hx,hy,hz,h
        integer(kind=ikind) :: nq,i,j,k,r,count,napos,ii,jj
        integer(kind=ikind) :: spi, spj, spk, spr, szo

        napos=size(apos)
        allocate(posits(napos,(maxval(ll)+1)*(maxval(ll)+2)/2))
        do i=1,napos
            do j=0,(maxval(ll)+1)*(maxval(ll)+2)/2-1
            posits(i,j+1)=apos(i)+j
            end do
        end do
        print*,'posits created'


        szo = size(z1(:,1,1))        
        
        nq= size(q)
        tsi=0.0_dp

        !First big loop
        count=0
        do i=1,ncap
            do j=i+1,ncap
                do k=i+1,ncap
                    do r=k+1,ncap
                        hx = px(k, r) - px(i, j)
                        hy = py(k, r) - py(i, j)
                        hz = pz(k, r) - pz(i, j)
                        h = sqrt((hx * hx + hy * hy + hz * hz))

                        allocate(posI((ll(i)+1)*(ll(i)+2)/2),posJ((ll(j)+1)*(ll(j)+2)/2) &
                                ,posK((ll(k)+1)*(ll(k)+2)/2),posR((ll(r)+1)*(ll(r)+2)/2))


                        posI = posits(i,1:(ll(i)+1)*(ll(i)+2)/2)
                        posJ = posits(j,1:(ll(j)+1)*(ll(j)+2)/2)
                        posK = posits(k,1:(ll(k)+1)*(ll(k)+2)/2)
                        posR = posits(r,1:(ll(r)+1)*(ll(r)+2)/2)
                        
                        spi = size(posI)
                        spj = size(posJ)
                        spk = size(posK)
                        spr = size(posR)

                        allocate(zcontrred(spj, spk, spr, spi), za(szo, spj), zb(szo, spk), &
                       &         cmat(spj, spk), zcontrred2(spr, spi, spj, spk))

                          do ii = 1, spi
                            do jj = 1, spr
!                                za = transpose(z1(:,posj,posi(ii)))
                                za = z1(:,posj,posi(ii))
                                zb = z2(:,posk,posr(jj))
!                                cmat = matmul(za,zb)
                                call dgemm('t','n', spj, spk, szo, 1.0_dp/8.0_dp, za, &
                               &           szo, zb, szo, 0.0_dp, cmat, spj)
                                zcontrred(:,:,jj,ii) = cmat
                                zcontrred2(jj,ii,:,:) = cmat
                            enddo
                          enddo

!                        zcontrred=zcontrred/8.0
!                        zcontrred2=zcontrred2/8.0

                        allocate(dx1red(ll(i)+ll(j)+1,ll(j)+1,ll(i)+1),dy1red(ll(i)+ll(j)+1,ll(j)+1,ll(i)+1), &
                                dz1red(ll(i)+ll(j)+1,ll(j)+1,ll(i)+1))
                        allocate(dx2red(ll(k)+ll(r)+1,ll(r)+1,ll(k)+1),dy2red(ll(k)+ll(r)+1,ll(r)+1,ll(k)+1), &
                                dz2red(ll(k)+ll(r)+1,ll(r)+1,ll(k)+1))

                        dx1red=dx(:(ll(i)+ll(j)+1),:ll(j)+1,:ll(i)+1,j,i)
                        dy1red=dy(:(ll(i)+ll(j)+1),:ll(j)+1,:ll(i)+1,j,i)
                        dz1red=dz(:(ll(i)+ll(j)+1),:ll(j)+1,:ll(i)+1,j,i)
                        dx2red=dx(:(ll(k)+ll(r)+1),:ll(r)+1,:ll(k)+1,r,k)
                        dy2red=dy(:(ll(k)+ll(r)+1),:ll(r)+1,:ll(k)+1,r,k)
                        dz2red=dz(:(ll(k)+ll(r)+1),:ll(r)+1,:ll(k)+1,r,k)

                        if (h < cutoffcentre) then

                            call integral_ijkr_pzero(nq, ll(i), ll(j), ll(k), ll(r), p0matrix, dx1red, dy1red, &
                                    dz1red,dx2red,dy2red,dz2red, i, j, k, r, &
                                    zcontrred, zcontrred2, apos, cutoffz, cutoffmd, f)

                        else

                            call tot_integral_k_ijkr(q, ll(i), ll(j), ll(k), ll(r), hx, hy, hz, h, dx1red, dy1red, &
                                    dz1red,dx2red,dy2red,dz2red, i, j, k, r, &
                                    zcontrred, zcontrred2, apos, cutoffz, cutoffmd, f)



                        end if
                        tsi = tsi + 8.000 * f * e12(:, i, j) * e12(:, k, r)
                        count=count+1
                        deallocate(posI,posJ,posK,posR,za,zb,cmat)
                        deallocate(zcontrred, zcontrred2)
                        deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)


                    end do
                end do
            end do
        end do
        print*,'la rata '
        do i=1,ncap
            do j=i+1,ncap
                do r=i+1,ncap
                    hx = px(i, r) - px(i, j)
                    hy = py(i, r) - py(i, j)
                    hz = pz(i, r) - pz(i, j)
                    h = sqrt((hx * hx + hy * hy + hz * hz))
                    allocate(posI((ll(i)+1)*(ll(i)+2)/2),posJ((ll(j)+1)*(ll(j)+2)/2) &
                            ,posR((ll(r)+1)*(ll(r)+2)/2))
                        posI=posits(i,1:(ll(i)+1)*(ll(i)+2)/2)
                        posJ=posits(j,1:(ll(j)+1)*(ll(j)+2)/2)
                        posR=posits(r,1:(ll(r)+1)*(ll(r)+2)/2)

                        spi = size(posI)
                        spj = size(posJ)
                        spr = size(posR)

                        allocate(zcontrred(spj, spi, spr, spi), za(szo, spj), zb(szo, spi), &
                       &         cmat(spj, spi), zcontrred2(spr, spi, spj, spi))

                          do ii = 1, spi
                            do jj = 1, spr
!                                za = transpose(z1(:,posj,posi(ii)))
                                za = z1(:,posj,posi(ii))
                                zb = z2(:,posi,posr(jj))
!                                cmat = matmul(za,zb)
                                call dgemm('t','n', spj, spi, szo, 1.0_dp/8.0_dp, za, &
                               &           szo, zb, szo, 0.0_dp, cmat, spj)
                                zcontrred(:,:,jj,ii) = cmat
                                zcontrred2(jj,ii,:,:) = cmat
                            enddo
                          enddo

                    allocate(dx1red(ll(i)+ll(j)+1,ll(j)+1,ll(i)+1),dy1red(ll(i)+ll(j)+1,ll(j)+1,ll(i)+1),&
                            dz1red(ll(i)+ll(j)+1,ll(j)+1,ll(i)+1))
                    allocate(dx2red(ll(i)+ll(r)+1,ll(r)+1,ll(i)+1),dy2red(ll(i)+ll(r)+1,ll(r)+1,ll(i)+1),&
                            dz2red(ll(i)+ll(r)+1,ll(r)+1,ll(i)+1))

                    dx1red=dx(:(ll(i)+ll(j)+1),:ll(j)+1,:ll(i)+1,j,i)
                    dy1red=dy(:(ll(i)+ll(j)+1),:ll(j)+1,:ll(i)+1,j,i)
                    dz1red=dz(:(ll(i)+ll(j)+1),:ll(j)+1,:ll(i)+1,j,i)
                    dx2red=dx(:(ll(i)+ll(r)+1),:ll(r)+1,:ll(i)+1,r,i)
                    dy2red=dy(:(ll(i)+ll(r)+1),:ll(r)+1,:ll(i)+1,r,i)
                    dz2red=dz(:(ll(i)+ll(r)+1),:ll(r)+1,:ll(i)+1,r,i)
                    
!                        zcontrred=zcontrred/8.0
!                        zcontrred2=zcontrred2/8.0

                    if (h < cutoffcentre) then
                        call integral_ijkr_pzero(nq, ll(i), ll(j), ll(i), ll(r), p0matrix, dx1red, dy1red, &
                                    dz1red,dx2red,dy2red,dz2red, i, j, i, r, &
                               zcontrred,  zcontrred2, apos, cutoffz, cutoffmd, f)
                    else

                        call tot_integral_k_ijkr(q, ll(i), ll(j), ll(i), ll(r), hx, hy, hz, h, dx1red, dy1red, &
                                    dz1red,dx2red,dy2red,dz2red, i, j, i, r, &
                                zcontrred,  zcontrred2, apos, cutoffz, cutoffmd, f)



                    end if
                    tsi = tsi + 4.000 * f * e12(:, i, j) * e12(:, i, r)
                    count=count+1
                    deallocate(posI,posJ,posR,za,zb,cmat)
                    deallocate(zcontrred, zcontrred2)
                    deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

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

                    allocate(posI((ll(i)+1)*(ll(i)+2)/2),posK((ll(k)+1)*(ll(k)+2)/2) &
                            ,posR((ll(r)+1)*(ll(r)+2)/2))
                    posI=posits(i,1:(ll(i)+1)*(ll(i)+2)/2)
                    posK=posits(k,1:(ll(k)+1)*(ll(k)+2)/2)
                    posR=posits(r,1:(ll(r)+1)*(ll(r)+2)/2)

                    spi = size(posI)
                    spk = size(posK)
                    spr = size(posR)

                    allocate(zcontrred(spi, spk, spr, spi), za(szo, spi), zb(szo, spk), &
                   &         cmat(spi, spk), zcontrred2(spr, spi, spi, spk))

                    do ii = 1, spi
                        do jj = 1, spr
!                            za = transpose(z1(:,posi,posi(ii)))
                            za = z1(:,posi,posi(ii))
                            zb = z2(:,posk,posr(jj))
!                            cmat = matmul(za,zb)
                            call  dgemm('t','n', spi, spk, szo, 1.0_dp/8.0_dp, za, &
                           &           szo, zb, szo, 0.0_dp, cmat, spi)
                            zcontrred(:,:,jj,ii) = cmat
                            zcontrred2(jj,ii,:,:) = cmat
                        enddo
                    enddo
                    
                    allocate(dx1red(ll(i)+ll(i)+1,ll(i)+1,ll(i)+1),dy1red(ll(i)+ll(i)+1,ll(i)+1,ll(i)+1), &
                                dz1red(ll(i)+ll(i)+1,ll(i)+1,ll(i)+1))
                    allocate(dx2red(ll(k)+ll(r)+1,ll(r)+1,ll(k)+1),dy2red(ll(k)+ll(r)+1,ll(r)+1,ll(k)+1), &
                                dz2red(ll(k)+ll(r)+1,ll(r)+1,ll(k)+1))

                    dx1red=dx(:(ll(i)+ll(i)+1),:ll(i)+1,:ll(i)+1,i,i)
                    dy1red=dy(:(ll(i)+ll(i)+1),:ll(i)+1,:ll(i)+1,i,i)
                    dz1red=dz(:(ll(i)+ll(i)+1),:ll(i)+1,:ll(i)+1,i,i)
                    dx2red=dx(:(ll(k)+ll(r)+1),:ll(r)+1,:ll(k)+1,r,k)
                    dy2red=dy(:(ll(k)+ll(r)+1),:ll(r)+1,:ll(k)+1,r,k)
                    dz2red=dz(:(ll(k)+ll(r)+1),:ll(r)+1,:ll(k)+1,r,k)

!                    zcontrred=zcontrred/8.0
!                    zcontrred2=zcontrred2/8.0

                    if (h < cutoffcentre) then
                        call integral_ijkr_pzero(nq, ll(i), ll(i), ll(k), ll(r), p0matrix, dx1red, dy1red, &
                                    dz1red,dx2red,dy2red,dz2red, i, i, k, r, &
                                    zcontrred, zcontrred2, apos, cutoffz, cutoffmd, f)
                    else

                        call tot_integral_k_ijkr(q, ll(i), ll(i), ll(k), ll(r), hx, hy, hz, h, dx1red, dy1red, &
                                    dz1red,dx2red,dy2red,dz2red, i, i, k, r, &
                                    zcontrred, zcontrred2, apos, cutoffz, cutoffmd, f)



                    end if
                    tsi = tsi+ 4.000 * f * e12(:, i, i) * e12(:, k, r)
                    count=count+1
                    deallocate(posI,posK,posR,za,zb,cmat)
                    deallocate(zcontrred, zcontrred2)
                    deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

                end do
            end do
        end do


        do i=1,ncap
            do k=i+1,ncap

                hx = px(k, k) - px(i, i)
                hy = py(k, k) - py(i, i)
                hz = pz(k, k) - pz(i, i)
                h = sqrt((hx * hx + hy * hy + hz * hz))

                allocate(posI((ll(i)+1)*(ll(i)+2)/2),posK((ll(k)+1)*(ll(k)+2)/2))
                posI=posits(i,1:(ll(i)+1)*(ll(i)+2)/2)

                posK=posits(k,1:(ll(k)+1)*(ll(k)+2)/2)

                spi = size(posI)
                spk = size(posK)

                allocate(zcontrred(spi, spk, spk, spi), za(szo, spi), zb(szo, spk), &
               &         cmat(spi, spk), zcontrred2(spk, spi, spi, spk))

                do ii = 1, spi
                    do jj = 1, spk
!                        za = transpose(z1(:,posi,posi(ii)))
                        za = z1(:,posi,posi(ii))
                        zb = z2(:,posk,posk(jj))
!                        cmat = matmul(za,zb)
                        call dgemm('t','n', spi, spk, szo, 1.0_dp/8.0_dp, za, &
                       &           szo, zb, szo, 0.0_dp, cmat, spi)
                        zcontrred(:,:,jj,ii) = cmat
                        zcontrred2(jj,ii,:,:) = cmat
                    enddo
                enddo
                
                allocate(dx1red(ll(i)+ll(i)+1,ll(i)+1,ll(i)+1),dy1red(ll(i)+ll(i)+1,ll(i)+1,ll(i)+1), &
                                dz1red(ll(i)+ll(i)+1,ll(i)+1,ll(i)+1))
                allocate(dx2red(ll(k)+ll(k)+1,ll(k)+1,ll(k)+1),dy2red(ll(k)+ll(k)+1,ll(k)+1,ll(k)+1), &
                                dz2red(ll(k)+ll(k)+1,ll(k)+1,ll(k)+1))

                dx1red=dx(:(ll(i)+ll(i)+1),:ll(i)+1,:ll(i)+1,i,i)
                dy1red=dy(:(ll(i)+ll(i)+1),:ll(i)+1,:ll(i)+1,i,i)
                dz1red=dz(:(ll(i)+ll(i)+1),:ll(i)+1,:ll(i)+1,i,i)
                dx2red=dx(:(ll(k)+ll(k)+1),:ll(k)+1,:ll(k)+1,k,k)
                dy2red=dy(:(ll(k)+ll(k)+1),:ll(k)+1,:ll(k)+1,k,k)
                dz2red=dz(:(ll(k)+ll(k)+1),:ll(k)+1,:ll(k)+1,k,k)

!                zcontrred=zcontrred/8.0
!                zcontrred2=zcontrred2/8.0

                if (h < cutoffcentre) then
                    call integral_ijkr_pzero(nq, ll(i), ll(i), ll(k), ll(k), p0matrix, dx1red, dy1red, &
                                    dz1red,dx2red,dy2red,dz2red, i, i, k, k, &
                                    zcontrred, zcontrred2, apos, cutoffz, cutoffmd, f)
                else

                    call tot_integral_k_ijkr(q, ll(i), ll(i), ll(k), ll(k), hx, hy, hz, h, dx1red, dy1red, &
                                    dz1red,dx2red,dy2red,dz2red, i, i, k, k, &
                                    zcontrred, zcontrred2, apos, cutoffz, cutoffmd, f)
                end if
                tsi = tsi+ 2.000 * f * e12(:, i, i) * e12(:, k, k)
                deallocate(posI,posK,za,zb,cmat)
                deallocate(zcontrred, zcontrred2)
                deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

            end do
        end do




        do i=1,ncap
            allocate(posI((ll(i)+1)*(ll(i)+2)/2))
            posI=posits(i,1:(ll(i)+1)*(ll(i)+2)/2)


            spi = size(posI)

            allocate(zcontrred(spi, spi, spi, spi), za(szo, spi), zb(szo, spi), &
           &         cmat(spi, spi), zcontrred2(spi, spi, spi, spi))

            do ii = 1, spi
                do jj = 1, spi
!                    za = transpose(z1(:,posi,posi(ii)))
                    za = z1(:,posi,posi(ii))
                    zb = z2(:,posi,posi(jj))
!                    cmat = matmul(za,zb)
                    call dgemm('t','n', spi, spi, szo, 1.0_dp/8.0_dp, za, &
                   &           szo, zb, szo, 0.0_dp, cmat, spi)
                    zcontrred(:,:,jj,ii) = cmat
                    zcontrred2(jj,ii,:,:) = cmat
                enddo
            enddo
            
            allocate(dx1red(ll(i)+ll(i)+1,ll(i)+1,ll(i)+1),dy1red(ll(i)+ll(i)+1,ll(i)+1,ll(i)+1), &
                                dz1red(ll(i)+ll(i)+1,ll(i)+1,ll(i)+1))
            allocate(dx2red(ll(i)+ll(i)+1,ll(i)+1,ll(i)+1),dy2red(ll(i)+ll(i)+1,ll(i)+1,ll(i)+1), &
                                dz2red(ll(i)+ll(i)+1,ll(i)+1,ll(i)+1))

            dx1red=dx(:(ll(i)+ll(i)+1),:ll(i)+1,:ll(i)+1,i,i)
            dy1red=dy(:(ll(i)+ll(i)+1),:ll(i)+1,:ll(i)+1,i,i)
            dz1red=dz(:(ll(i)+ll(i)+1),:ll(i)+1,:ll(i)+1,i,i)
            dx2red=dx(:(ll(i)+ll(i)+1),:ll(i)+1,:ll(i)+1,i,i)
            dy2red=dy(:(ll(i)+ll(i)+1),:ll(i)+1,:ll(i)+1,i,i)
            dz2red=dz(:(ll(i)+ll(i)+1),:ll(i)+1,:ll(i)+1,i,i)

!            zcontrred=zcontrred/8.0
!            zcontrred2=zcontrred2/8.0

            call integral_ijkr_pzero(nq, ll(i), ll(i), ll(i), ll(i), p0matrix, dx1red, dy1red, &
                                    dz1red,dx2red,dy2red,dz2red, i, i, i, i, &
                                    zcontrred, zcontrred2, apos, cutoffz, cutoffmd, f)

            tsi = tsi + f * e12(:, i, i) * e12(:, i, i)
            count=count+1
            deallocate(posI,za,zb,cmat)
            deallocate(zcontrred, zcontrred2)
            deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

        end do

    print*, count

    end subroutine integration



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
        real(kind=dp),dimension(:,:,:),allocatable :: z1,z2
        real(kind=dp), intent(in), dimension(:,:,:,:) :: p0matrix
        real(kind=dp), intent(in) :: cutoffmd, cutoffz,cutoffcentre

        real(kind=dp), intent(out), dimension(nq):: result

         real(kind=dp),  dimension(maxl*2+1,maxl+1,maxl+1,napos,napos) :: ddx,ddy,ddz
        real(kind=dp), dimension(napos,napos) :: px,py,pz
        real(kind=dp),  dimension(:,:,:,:), allocatable :: zcontr
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
        allocate(zcontr(nipos,nipos,nipos,nipos))
        !allocate(z1(size(m1), nipos, nipos), z2(size(m1), nipos, nipos))

        nmat=size(m1)


        call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,e12,ll,maxl, ipos,nipos,apos,napos,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q,nq,list1,listN1,list2,listN2)

        call cpu_time(time3)
        print*,'Time variables', time3-time2

        call integration(napos,px,py,pz,ll,p0matrix,ddx,ddy,ddz,z1,z2,apos,cutoffz,cutoffmd, cutoffcentre,q,e12,result)

        call cpu_time(time4)

         print*,'Time calc', time4-time3

        end subroutine total_scattering_calculation

        end module integrals_ijkr
