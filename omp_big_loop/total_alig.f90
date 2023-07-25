module total_alig
    use bessel_calcs
    implicit none
    contains

    subroutine tot_integration_aligned(ncap,px,py,pz,l,m,n,dx,dy,dz,z1,z2,group_start,group_count,group, &
            cutoffz,cutoffmd, cutoffcentre,q_al,e12,exponent1,exponent2,tsi)



        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(kind=ikind), INTENT(IN) :: ncap
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:) :: l,m,n,group,group_start,group_count

        REAL(kind=dp), intent(in), dimension(:,:,:,:,:),target :: dx,dy,dz
        complex(kind=dp), intent(in), dimension(:,:,:,:):: exponent1,exponent2
        real(kind=dp),dimension(:,:,:), pointer :: dx1,dy1,dz1,dx2,dy2,dz2

        REAL(kind=dp), intent(in), dimension(:,:,:) ::e12
        real(kind=dp), intent(in), dimension(:,:,:)::z1,z2
        real(kind=dp), dimension(:,:,:,:), allocatable::zcontrred,zcontrred2
        REAL(kind=dp), intent(in), dimension(:,:) :: px,py,pz
        REAL(kind=dp), intent(in), dimension(:,:) :: q_al
        REAL(kind=dp), intent(in) :: cutoffz, cutoffmd,cutoffcentre


        complex(kind=dp), dimension(size(q_al(1,:))) :: f
        integer(kind=ikind), dimension(:,:), allocatable :: posits
        real(kind=dp),dimension(:,:), allocatable :: za,zb,cmat
        integer(kind=ikind),dimension(:), allocatable ::posi,posj,posk,posr
        Complex(kind=dp), intent(out), dimension(size(q_al(1,:))) :: tsi
        real(kind=dp) :: hx,hy,hz,h
        integer(kind=ikind),dimension(size(l)) :: ll
        integer(kind=ikind) :: nq,i,j,k,r,count,ng,ii,jj
        integer(kind=ikind) :: spi, spj, spk, spr, szo,nt

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



        szo = size(z1(:,1,1))

        nq= size(q_al(1,:))


        !First big loop

        tsi=0.0_dp


        do i=1,ncap
            do j=1,ncap
                do k=1,ncap
                    do r=1,ncap
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

                        allocate(zcontrred(spj, spk, spr, spi), za(szo, spj), zb(szo, spk), &
                       &         cmat(spj, spk), zcontrred2(spr, spi, spj, spk))

                          do ii = 1, spi
                            do jj = 1, spr
                                za = z1(:,posj,posi(ii))
                              !  za = transpose(z1(:,posj,posi(ii)))
                                zb = z2(:,posk,posr(jj))
                               ! cmat = matmul(za,zb)
                                call dgemm('t','n', spj, spk, szo, 1.0_dp/8.0_dp, za, &
                              &           szo, zb, szo, 0.0_dp, cmat, spj)
                                zcontrred(:,:,jj,ii) = cmat

                                za = z2(:,posj,posi(ii))
                              !  za = transpose(z1(:,posj,posi(ii)))
                                zb = z1(:,posk,posr(jj))
                               ! cmat = matmul(za,zb)
                                call dgemm('t','n', spj, spk, szo, 1.0_dp/8.0_dp, za, &
                              &           szo, zb, szo, 0.0_dp, cmat, spj)

                                zcontrred2(jj,ii,:,:) = cmat
                            enddo
                          enddo



                         deallocate(posI,posJ,posK,posR,za,zb,cmat)
                         dx1=>dx(:,:,:,j,i)
                         dy1=>dy(:,:,:,j,i)
                         dz1=>dz(:,:,:,j,i)
                         dx2=>dx(:,:,:,r,k)
                         dy2=>dy(:,:,:,r,k)
                         dz2=>dz(:,:,:,r,k)

                        call tot_integral_k_ijkr_alig(q_al, l,m,n,group_start, group_count, hx, hy, hz, h, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, k, r, &
                                zcontrred,  zcontrred2,  cutoffz, cutoffmd,f,exponent1,exponent2)





                        tsi = tsi + f * e12(:, i, j) * e12(:, k, r)


                        deallocate(zcontrred, zcontrred2)
               !        deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)


                    end do
                end do
            end do
        end do
    end subroutine tot_integration_aligned

     subroutine tot_integral_k_ijkr_alig(q_al,l,m,n,group_start,group_count,hx,hy,hz,h,dx1,dy1,dz1,dx2,dy2,dz2, gi,gj,gk, gr,&
            zcontr, zcontr2, &
            cutoff1, cutoff2,f,exp1,exp2)



        implicit none


        !INTEGER, PARAMETER :: dp = selected_real_kind(2*precision(1.0_dp))
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=selected_int_kind(8)), intent(in)  :: gi,gj,gk,gr
        integer(kind=selected_int_kind(8)), dimension(:), intent(in) :: l,m,n,group_start,group_count
        real(kind=dp), intent(in)              :: cutoff1,cutoff2, hx, hy, hz, h
        real(kind=dp), intent(in), dimension(:,:,:,:) :: zcontr,zcontr2
        complex(kind=dp),intent(in),dimension(:,:,:,:) :: exp1,exp2
        real(kind=dp), intent(in), dimension(:,:)   ::  q_al

        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1,dx2,dy2,dz2
        !real(kind=dp), pointer,dimension(:,:,:) :: dxx,dyy,dzz, dxx2,dyy2,dzz2



        integer(kind=selected_int_kind(8)) :: ka, posi, posj, posk, posr, ra,i,j,k,r
        integer(kind=selected_int_kind(8)) ::  h1
        integer(kind=selected_int_kind(8)) ::  ll, mm, nn, llp, mmp, nnp,llmax

        real(kind=dp)  ::coeff,prodd,ztot,mdn, mdl, mdm, mdlp, mdmp,mdnp,z11,z22
        INTEGER(kind=ikind), parameter :: dim = 13
        real(kind=dp) :: prod1,prod2,prod3,prod4,prod5,prod6
        REAL(kind=dp), DIMENSION(dim,dim)           :: a, b, c
        complex(kind=dp), DIMENSION(size(q_al(1,:)))               :: h_saved
        REAL(kind=dp), DIMENSION(dim, dim, dim, dim)  :: h_pre2
        REAL(kind=dp), DIMENSION(dim)               :: BD



      !  real(kind=dp), dimension(:), allocatable :: pmu, h_0, h_1, h_r, muoh,zij, zij2, zkr, zkr2

        complex(kind=dp), intent(out), dimension(size(q_al(1,:))) :: f
        real(kind=dp), allocatable,dimension(:) :: l1vec,l2vec,suml1l2


        LLmax = l(group_start(gi)) + m(group_start(gi)) + n(group_start(gi)) + &
                l(group_start(gj)) + m(group_start(gj)) + n(group_start(gj)) + &
                l(group_start(gk)) + m(group_start(gk)) + n(group_start(gk)) + &
                l(group_start(gr)) + m(group_start(gr)) + n(group_start(gr))


        if (LLmax + 1 > dim) then
            print*, "only s,p,d and f type GTOs are supported"
            stop
        end if




       ! posi=apos(i)



        h_saved=(0.0_dp,0.0_dp)

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
                                    H1=1
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
                                                    h_saved = h_saved + exp1(:,llp+1,mmp+1,nnp+1)*exp2(:,ll+1,mm+1,nn+1)*prod6

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

       f=h_saved*exp((0,1.00_dp)*q_al(1,:)*(Hx))*exp((0,1.0_dp)*q_al(2,:)*(Hy))*exp((0,1.0_dp)*q_al(3,:)*(Hz))

    end subroutine

end module total_alig