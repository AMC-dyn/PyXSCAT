Module elastic_alig
    use bessel_calcs
    implicit none
    contains
     subroutine elastic_integration_alig(ncap,px,py,pz,l,m,n,dx,dy,dz,z,group_start,group_count,group, &
            cutoffz,cutoffmd, cutoffcentre,q_al,e12,exponent1,tsi)


        use omp_lib
        implicit none

              INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(kind=ikind), INTENT(IN) :: ncap
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:) :: l,m,n,group,group_start,group_count

        REAL(kind=dp), intent(in), dimension(:,:,:,:,:),target :: dx,dy,dz
        complex(kind=dp), intent(in), dimension(:,:,:,:) :: exponent1
        real(kind=dp),dimension(:,:,:), pointer :: dx1,dy1,dz1,dx2,dy2,dz2

        REAL(kind=dp), intent(in), dimension(:,:,:) ::e12
        real(kind=dp), intent(in), dimension(:,:)::z

        REAL(kind=dp), intent(in), dimension(:,:) :: px,py,pz
        REAL(kind=dp), intent(in), dimension(:,:) :: q_al
        REAL(kind=dp), intent(in) :: cutoffz, cutoffmd,cutoffcentre


        complex(kind=dp), dimension(size(q_al(1,:))) :: f
        integer(kind=ikind), dimension(:,:), allocatable :: posits
        real(kind=dp),dimension(:,:), allocatable :: za,zb
        integer(kind=ikind),dimension(:), allocatable ::posi,posj,posk,posr
        complex(kind=dp), intent(out), dimension(size(q_al(1,:))) :: tsi
        real(kind=dp) :: hx,hy,hz,h
        integer(kind=ikind),dimension(size(l)) :: ll
        integer(kind=ikind) :: nq,i,j,k,r,count,ng,ii,jj
        integer(kind=ikind) :: spi, spj, spk, spr,nt,counter1,counter2

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
        print*,'posits created'
        print*,size(group_start),sum(group_count)




        nq= size(q_al(1,:))

      !  print*,OMP_get_num_threads()
      !  call omp_set_num_threads(16)
      !  print*,OMP_get_num_threads()
        !First big loop

        tsi=(0.0_dp,0.0_dp)
     !  !$OMP PARALLEL do private(posI,posJ,spi,spj,za), &
      !  !$OMP& private(f,ii,jj,h,hx,hy,hz,i,j,dx1,dy1,dz1) shared(q_al,l,m,n), &
     !   !$OMP& shared( cutoffz, posits,cutoffmd,group_count,group_start,z) REDUCTION(+:tsi)
        print*,Z(2,9), E12(1,16,16)
        counter1=0

        do i=1,ncap
             counter1=counter1+group_count(i)
             counter2=0
            do j=1,ncap

                        counter2=counter2+group_count(j)
                        hx =  px(i, j)
                        hy =  py(i, j)
                        hz =  pz(i, j)
                        h = (hx * hx + hy * hy + hz * hz)**0.5

                        allocate(posI(size(posits(i,:group_count(i)))),posJ(size(posits(j,:group_count(j)))))


                        posI = posits(i,:group_count(i))
                        posJ = posits(j,:group_count(j))


                        spi = size(posI)
                        spj = size(posJ)


                        allocate(za(spi, spj))

                         za=z(posI,posJ)


                         deallocate(posI,posJ)
                         dx1=>dx(:,:,:,j,i)
                         dy1=>dy(:,:,:,j,i)
                         dz1=>dz(:,:,:,j,i)



                        call integral_k_ijkr_alig(q_al,exponent1,l,m,n,group_start, group_count, hx, hy, hz, h,dx1,dy1,dz1,i, j &
                                    ,za, cutoffz, cutoffmd,f)






                        tsi = tsi + f * e12(:, i, j)

                        count=count+1
                        deallocate(za)
               !        deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)


                    end do
                end do
        print*,tsi(1),q_al(1,1), q_al(2,1), q_al(3,1)
      ! !$OMP END parallel DO



    end subroutine elastic_integration_alig
     subroutine integral_k_ijkr_alig(q_al,exponent1,l,m,n,group_start,group_count,hx,hy,hz,h,dx1,dy1,dz1, gi,gj,&
            za,&
            cutoff1, cutoff2,f)



        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        !INTEGER, PARAMETER :: dp = selected_real_kind(2*precision(1.0_dp))
        integer(kind=selected_int_kind(8)), intent(in)  :: gi,gj
        integer(kind=selected_int_kind(8)), dimension(:), intent(in) :: l,m,n,group_start,group_count
        real(kind=dp), intent(in)              :: cutoff1,cutoff2, hx, hy, hz, h
        real(kind=dp), intent(in), dimension(:,:) :: za
        real(kind=dp), intent(in), dimension(:,:)   ::  q_al

        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1
        complex(kind=dp), intent(in),  dimension(:,:,:,:) :: exponent1
        !real(kind=dp), pointer,dimension(:,:,:) :: dxx,dyy,dzz, dxx2,dyy2,dzz2



        integer(kind=selected_int_kind(8)) :: ka, posi, posj, posk, posr, ra,i,j,k,r
        integer(kind=selected_int_kind(8)) ::  h1
        integer(kind=selected_int_kind(8)) ::  ll, mm, nn, llp, mmp, nnp,llmax

        real(kind=dp)  ::coeff,prodd,ztot,mdn, mdl, mdm, mdlp, mdmp,mdnp,z1
        INTEGER(kind=ikind), parameter :: dim = 13
        real(kind=dp) :: prod1,prod2,prod3,prod4,prod5,prod6
        REAL(kind=dp), DIMENSION(dim,dim)           :: a, b, c
        complex(kind=dp), DIMENSION(size(q_al(1,:)))               :: h_saved
        REAL(kind=dp), DIMENSION(dim, dim, dim, dim)  :: h_pre2
        REAL(kind=dp), DIMENSION(dim)               :: BD



      !  real(kind=dp), dimension(:), allocatable :: pmu, h_0, h_1, h_r, muoh,zij, zij2, zkr, zkr2

        complex(kind=dp), intent(out), dimension(size(q_al(1,:))) :: f
        real(kind=dp), allocatable,dimension(:) :: l1vec,l2vec,suml1l2








        h_saved=(0.0_dp,0.0_dp)

! loop through all possible ways to get total angular momentum lmax1


        posI=1


        ! loop through all possible ways to get total angular momentum lmax1
        do i = group_start(gi), group_start(gi) + group_count(gi) - 1
            posj=1
            do j = group_start(gj), group_start(gj) + group_count(gj) - 1

                ztot=Za(posi,posj)


!
                        ! continue only if larger
                     !   if (abs(ztot) < cutoff1) then
                      !          posj=posj+1
                       !         cycle
                      !      end if
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
                                    H1=1.0
                                    MDN=Dz1(nn+1,n(i)+1,n(j)+1)
                                    if (abs(MDN)<1.0e-30) cycle ! check if MD coeff is 0
                                    prod3 = MDN  * prod2

                                    if (abs(prod3)>1E-30) then
                                    ! MD coeff 4
                                         h_saved = h_saved + exponent1(:,ll+1,mm+1,nn+1)*prod3

                                    end if

                                end do
                            end do
                           end do

                    posj=posj+1
            end do
        posi=posi+1
        end do

      F=h_saved*exp((0,1.00_dp)*q_al(1,:)*(Hx))*exp((0,1.0_dp)*q_al(2,:) &
              *(Hy))*exp((0,1.0_dp)*q_al(3,:)*(Hz))

    end subroutine

End Module elastic_alig