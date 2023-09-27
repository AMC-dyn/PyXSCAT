Module TSj0contr
    use Bessels_j0
    use p0_cases
    use twordms
    use MD
    implicit none
contains

    subroutine total_scattering_j0_ncontr(q, l, m, n, ngto, ng, nq, maxl, typec, state1, &
            state2, ncontr, group, gs, gf, gc, confs, ga, xx, yy, zz, coeffs, mmod, civecs, geom, &
            cutoffmd, cutoffz, cutoffcentre, result2)
        implicit none
        !Precision parameters
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

        !Entering variables from readers
        integer(kind = ikind), intent(in) :: ngto, ng, nq, maxl, typec, state1, state2, ncontr
        integer(kind = ikind), intent(in), dimension(:), allocatable :: l, m, n, group, gs, gf, gc
        integer(kind = ikind), dimension(:, :), allocatable, intent(in) :: confs
        real(kind = dp), intent(in), dimension(:), allocatable :: ga, xx, yy, zz, q, coeffs
        real(kind = dp), intent(in), dimension(:, :), allocatable :: mmod, civecs, geom
        real(kind = dp), intent(in) :: cutoffmd, cutoffz, cutoffcentre

        !twordm variables
        integer(kind = ikind), dimension(:), allocatable :: m1, m2, m3, m4
        integer(kind = ikind), dimension(:, :), allocatable :: mat, ep3, ndiff2
        integer(kind = ikind) :: nmat, i, j
        real(kind = dp), dimension(:), allocatable :: total

        !Variables to create total_variables
        real(kind = dp), DIMENSION(size(q), 4 * maxval(l) + 1, 4 * maxval(l) + 1, 4 * maxval(l) + 1) :: P0matrix
        real(kind = dp), dimension(maxl * 2 + 1, maxl + 1, maxl + 1, ngto, ngto) :: ddx, ddy, ddz

        real(kind = dp), dimension(:,:),allocatable :: px, py, pz
        real(kind = dp), dimension(:, :, :),allocatable :: e12


        !Result out
        real(kind = dp), allocatable, dimension(:), intent(out) :: Result2


        P0matrix = 0.0_dp
        call set_P0(P0matrix, 4 * maxval(l), q)

        call maxcoincidence(confs, ep3, ndiff2)

        call createtwordm(confs, civecs, ndiff2, ep3, mat, total, state1, state2)

        call system('python3 Zcontr.py')
        !allocate(m1(size(mat(:, 1))), m2(size(mat(:, 1))), m3(size(mat(:, 1))), m4(size(mat(:, 1))))
        !m1 = mat(:, 1)
        !m2 = mat(:, 2)
        !m3 = mat(:, 3)
        !m4 = mat(:, 4)
        !nmat = size(m1)
        call variables_total_big(px, py, pz, ddx, ddy, ddz, &
                e12, maxl, ngto, ga, l, m, n, xx, yy, zz, &
                mmod, nmat, q,coeffs, nq)

        call tot_integration_contract(ncontr,coeffs, px, py, pz, l, m, n, p0matrix, ddx, ddy, ddz, &
                e12,gs,gf,gc,cutoffz, cutoffmd, cutoffcentre, q, result2)

    end subroutine


     subroutine variables_total_big(px,py,pz,ddx,ddy,ddz,&
        e12,maxl,ngto,ga,l,m,n,xx,yy,zz, &
        mmod,nmat,q,coeffs,nq)


        implicit none


        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=ikind), intent(in)::ngto, nmat, nq, maxl
        integer(kind=ikind), intent(in),dimension(:),allocatable :: l, m, n


        real(kind=dp), intent(in),dimension(:),allocatable :: ga, xx, yy, zz,q,coeffs

        real(kind=dp), intent(in),dimension(:,:),allocatable :: mmod
        real(kind=dp), intent(out), dimension(:,:), allocatable:: px,py,pz


        !real(kind=dp), intent(out), dimension(ngto,ngto,ngto,ngto) :: zcontr

        real(kind=dp), intent(out), dimension(:,:,:),allocatable :: e12

        real(kind=dp), intent(out), dimension(maxl*2+1,maxl+1,maxl+1,ngto,ngto) :: ddx,ddy,ddz


        integer(kind=ikind) :: max4l, max2l, i, j, k, N1, N2, nn1, nn2,ii,jj,ls, ms, ns,count
        integer(kind=ikind), dimension(ngto) :: ll
        real(kind=dp) :: pi, gap,time1,time2,time3,time4

        integer(kind=ikind), dimension(:), allocatable :: iduplicates,jduplicates,m11,m22,m33,m44
        real(kind=dp),  dimension(ngto,ngto,maxl*2+1,maxl+1,maxl+1):: dx,dy,dz
        real(kind=dp),  dimension(nq):: preexp
        real(kind=dp),  dimension(:), allocatable:: temp1,temp2

        integer(kind=ikind),dimension(:,:),allocatable :: matfin
        real(kind=dp),dimension(:), allocatable :: totalfin
        real(kind=dp)   :: obwohl
        integer(kind=ikind) :: counter
        logical :: divided

        pi = acos(-1.0000)
        max4l=maxval(l)*4
        max2l=maxval(l)*2+1




        allocate(e12(nq,ngto,ngto), px(ngto,ngto), py(ngto,ngto), pz(ngto,ngto))
        gap=0.0
        px=0.0
        py=0.0
        pz=0.0
        e12=0.0
        ddx=0.0
        ddy=0.0
        ddz=0.0
        dx=0.0
        dy=0.0
        dz=0.0
        print*,'hijos de la rata '



        call cpu_time(time3)

        print*,'hijos de la rata I'

        call fill_md_table(dx,l,xx,ga)
        call fill_md_table(dy,m,yy,ga)
        call fill_md_table(dz,n,zz,ga)

         print*,'hijos de la rata II'

       do j = 1,Ngto
            ! essentially taking the values for the first GTO in the group.
            ! All gtos in the group have the same prefactors as the prefactors
            ! do not't depend on l, m and n

            do i = 1, Ngto

                gaP=ga(i)+ga(j)
                Px(i,j)=(ga(i)*xx(i) + ga(j)*xx(j))/gaP
                Py(i,j)=(ga(i)*yy(i) + ga(j)*yy(j))/gaP
                Pz(i,j)=(ga(i)*zz(i) + ga(j)*zz(j))/gaP

                E12(:,i,j) = (pi/gaP)**1.5 * exp(-q*q*0.25/gaP) &
                    * exp(-ga(i)*ga(j)/gaP*((xx(i)-xx(j))**2. + (yy(i)-yy(j))**2. + (zz(i)-zz(j))**2.))*coeffs(i)*coeffs(j)
            end do
        end do



     do j = 1, Ngto

            do i = 1, Ngto

                do ls = 1, l(i)+l(j)+1
                    Ddx(ls, l(i)+1,l(j)+1,j,i)=Dx(i,j,ls, l(i)+1,l(j)+1)
                end do

                do ms = 1, m(i)+m(j)+1
                    Ddy(ms, m(i)+1,m(j)+1,j,i)=Dy(i,j,ms, m(i)+1,m(j)+1)
                end do

                do ns = 1, n(i)+n(j)+1
                    Ddz(ns, n(i)+1,n(j)+1,j,i)=Dz(i,j,ns, n(i)+1,n(j)+1)
                end do
            end do
        end do

    print*,ngto,size(E12(:,1,1)),size(E12(1,:,1)), size(E12(1,1,:))
    print*, 'leaving variables total'

    end subroutine variables_total_big

    SUBROUTINE tot_integral_ijkr_pzero_read(nq,l,m,n,p0mat,dx1,dy1,dz1,dx2,dy2,dz2,i,j,k,r, &
             cutoff1,cutoff2,f)



        implicit none

        ! definition of input
              INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER(kind=ikind), INTENT(IN)                       :: nq,i,j,k,r
        REAL(kind=dp), INTENT(IN)                             :: cutoff1, cutoff2
        REAL(kind=dp), INTENT(IN), DIMENSION(:,:,:,:)               :: p0mat
        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1,dx2,dy2,dz2



        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:)         :: l,m,n
        ! definition of output
        REAL(kind=dp), INTENT(OUT), DIMENSION(nq)             :: f
        ! definition of loop indices
        INTEGER(kind=ikind)                                   :: h1
        INTEGER(kind=ikind)                                   :: ll, mm, nn, llp, mmp, nnp
        ! definition of internal variables


        real(dp), external :: ddot
        REAL(kind=dp)                                         :: mdl, mdm, mdn, mdlp, mdmp,mdnp
        real(kind=dp) :: prod6, prod5,prod4,prod3,prod2,prod1
!        REAL(kind=dp), DIMENSION(size(Z1(:,1,1)))                          :: zij1, zkr1
!        REAL(kind=dp), DIMENSION(size(Z2(:,1,1)))                          :: zij2, zkr2



       ! posI=apos(i)

        f=0.0
        ! loop through all possible ways to get total angular momentum lmax1

!

                        ! continue only if larger

                     do ll = 0, l(i)+l(j)
                            MDL = Dx1(ll+1,l(i)+1,l(j)+1)
                            if (abs(MDL)<1.0e-30) cycle
                            prod1 = MDL
                            ! MD coeff 2
                            do mm = 0, m(i)+m(j)
                                MDM = Dy1(mm+1,m(i)+1,m(j)+1)
                                if (abs(MDM)<1.0e-30) cycle
                                prod2 = MDM * prod1
                                ! MD coeff 3
                                do nn =0, n(i)+n(j)
                                    H1=(-1.0)**(ll+mm+nn)
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
                                                    F = F + prod6 * P0mat(:,ll+llp+1,mm+mmp+1,nn+nnp+1)
                                                      !  print*,prod6, maxval(abs(P0mat(:,ll+llp+1,mm+mmp+1,nn+nnp+1))), ll+llp+1,mm+mmp+1,nn+nnp+1
                                                    !Int=Int+prodD.*f;

                                                end if
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do

     if (maxval(abs(F))>1E20) then
            print*, 'pzero is the problem'

        end if

    END SUBROUTINE

    subroutine tot_integral_k_ijkr_read(mu,l,m,n,h,hx,hy,hz,dx1,dy1,dz1,dx2,dy2,dz2, i,j,k,r,&
            BesselQ,BesselQ2,cutoff1, cutoff2,f)



        implicit none


        !INTEGER, PARAMETER :: dp = selected_real_kind(2*precision(1.0_dp))
              INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=selected_int_kind(8)), intent(in)  :: i,j,k,r
        integer(kind=selected_int_kind(8)), dimension(:), intent(in) :: l,m,n
        real(kind=dp), intent(in)              :: cutoff1,cutoff2,hx, hy, hz, h

        real(kind=dp), intent(in), dimension(:),allocatable   ::  mu

        real(kind=dp), intent(in),  dimension(:,:,:) :: dx1,dy1,dz1,dx2,dy2,dz2
        !real(kind=dp), pointer,dimension(:,:,:) :: dxx,dyy,dzz, dxx2,dyy2,dzz2


        real(kind=dp), dimension(:,:), intent(in):: BesselQ,BesselQ2
        integer(kind=selected_int_kind(8)) :: ka, ra
        integer(kind=selected_int_kind(8)) ::  h1
        integer(kind=selected_int_kind(8)) ::  ll, mm, nn, llp, mmp, nnp,llmax,ii,jj,kk

        real(kind=dp)  ::coeff,prodd,ztot,mdn, mdl, mdm, mdlp, mdmp,mdnp,z11,z22
        INTEGER(kind=ikind), parameter :: dim = 20
        real(kind=dp) :: prod1,prod2,prod3,prod4,prod5,prod6
        REAL(kind=dp), DIMENSION(dim,dim)           :: a, b, c
        REAL(kind=dp), DIMENSION(dim)               :: h_saved
        REAL(kind=dp), DIMENSION(dim, dim, dim, dim)  :: h_pre2
        REAL(kind=dp), DIMENSION(dim)               :: BD



      !  real(kind=dp), dimension(:), allocatable :: pmu, h_0, h_1, h_r, muoh,zij, zij2, zkr, zkr2

        real(kind=dp), intent(out), dimension(size(mu)) :: f
        real(kind=dp), allocatable,dimension(:) :: l1vec,l2vec,suml1l2


        LLmax = l(i) + m(i) + n(i) + &
                l(j) + m(j) + n(j) + &
                l(k) + m(k) + n(k) + &
                l(r) + m(r) + n(r)


        if (LLmax + 1 > dim) then
            print*, "only s, p, d, f, and g type GTOs are supported"
            stop
        end if

        a=0.0_dp
        b=0.0_dp
        c=0.0_dp

        call rrdj0(LLmax, Hx,a)
        call rrdj0(LLmax, Hy,b)
        call rrdj0(LLmax, Hz,c)




        bd=0.0_dp

        do kk = 0, LLmax
            do jj = 0, LLmax - kk
                do ii = 0, LLmax - kk - jj
                    call BesselDeriv(BD, kk, jj, ii, a, b, c, LLmax)
                    h_pre2(:,kk+1,jj+1,ii+1) = BD
                end do
            end do
        end do

       ! posi=apos(i)



        h_saved=0.0_dp

! loop through all possible ways to get total angular momentum lmax1





        ! loop through all possible ways to get total angular momentum lmax1
!
                        ! continue only if larger

                         do ll = 0, l(i)+l(j)
                            MDL = Dx1(ll+1,l(i)+1,l(j)+1)
                          !  if (abs(MDL)<1.0e-30) cycle
                            prod1 = MDL
                            ! MD coeff 2
                            do mm = 0, m(i)+m(j)
                                MDM = Dy1(mm+1,m(i)+1,m(j)+1)
                              !  if (abs(MDM)<1.0e-30) cycle
                                prod2 = MDM * prod1
                                ! MD coeff 3
                                do nn =0, n(i)+n(j)
                                    H1=(-1.0)**(ll+mm+nn)
                                    MDN=Dz1(nn+1,n(i)+1,n(j)+1)
                                 !   if (abs(MDN)<1.0e-30) cycle ! check if MD coeff is 0
                                    prod3 = MDN * H1  * prod2
                                    ! MD coeff 4
                                    do llp = 0, l(k)+l(r)
                                        MDLp=Dx2(llp+1,l(k)+1,l(r)+1)
                                  !      if (abs(MDLp)<1.0e-30) cycle ! check if MD coeff is 0
                                        prod4 = MDLp * prod3
                                        ! MD coeff 5
                                        do mmp = 0, m(k)+m(r)
                                            MDMp = Dy2(mmp+1,m(k)+1,m(r)+1)
                                    !        if (abs(MDMp)<1.0e-30) cycle ! check if MD coeff is 0
                                            prod5 = MDMp * prod4
                                            ! MD coeff 6
                                            do nnp = 0, n(k)+n(r)
                                                MDNp = Dz2(nnp+1,n(k)+1,n(r)+1)
                                                prod6 = MDNp * prod5
                                                ! cutoff after MD
                                                if (abs(prod6)>cutoff2) then
                                                    ! add the contribution to the total
                                                    h_saved = h_saved + h_pre2(:,ll+llp+1,mm+mmp+1,nn+nnp+1)*prod6

                                                end if
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do


      !  CALL bessels0rr(F, llmax, mu, H,h_saved)

        CALL bessel0sum(F,llmax, BesselQ,BesselQ2,mu,H, h_saved)


        !CALL BesselSum(F, mu, H, LLmax, h_saved)


    end subroutine


     subroutine tot_integration_contract(ncontr,coeffs,px,py,pz,l,m,n,p0matrix,dx,dy,dz,&
            e12,gs,gf,gc, &
            cutoffz,cutoffmd, cutoffcentre,q,tsi)


        use omp_lib
        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

        real(kind=kind(1.d0)), external :: ddot
        INTEGER(kind=ikind), INTENT(IN) :: ncontr
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:) :: l,m,n,gs,gf,gc

        REAL(kind=dp), intent(in), dimension(:,:,:,:,:),target :: dx,dy,dz
        real(kind=dp), intent(in), dimension(:):: coeffs
        real(kind=dp),dimension(:,:,:), pointer :: dx1,dy1,dz1,dx2,dy2,dz2
        REAL(kind=dp), intent(in), dimension(:,:,:,:) :: p0matrix
        REAL(kind=dp), intent(in), dimension(:,:,:),allocatable ::e12

        real(kind=dp) :: zcontrred,zcontrred2
        REAL(kind=dp), intent(in), dimension(:,:),allocatable :: px,py,pz
        REAL(kind=dp), intent(in), dimension(:),allocatable :: q
        REAL(kind=dp), intent(in) :: cutoffz, cutoffmd,cutoffcentre


        REAL(kind=dp), dimension(size(q)) :: f
        integer(kind=ikind), dimension(:,:), allocatable :: posits
        integer(kind=ikind), dimension(:), allocatable :: indvec
        real(kind=dp),dimension(:,:), allocatable :: za,zb,cmat
        integer(kind=ikind),dimension(:), allocatable ::posi,posj,posk,posr
        REAL(kind=dp), intent(out), dimension(:),allocatable :: tsi
        real(kind=dp),dimension(size(q)):: tsi2
        real(kind=dp) :: hx,hy,hz,h
        integer(kind=ikind),dimension(:,:), allocatable:: bigvec
        integer(kind=ikind),dimension(size(l)) :: ll,veccounts
        integer(kind=ikind) :: nq,i,j,k,r,count,ng,ii,jj,kk,rr,iii
        integer(kind=ikind) :: spi, spj, spk, spr, szo,nt,ngto,dimvec,initj,initk,initr,ncontr2,newgto
        real(kind=dp), dimension(:,:,:,:), allocatable:: Zbig
        real(kind=dp), dimension(206,size(q)):: BesselQ, BesselQ2


        ll=l+m+n

        allocate(tsi(size(q)))


       print*,'looking for the intruder'

        nq= size(q)


     !
        !First big loop

        tsi=0.0_dp
        tsi2=0.0_dp

        if (any(isnan(p0matrix))) print*,'ouch'


         ngto=sum(gc)
          print*,'Main routine with ', ncontr, ' and ', ngto

         call dphrec(q,BesselQ,206,206,nq)

         do i=0,199
              BesselQ2(i+1,:)=(-1.0_dp)**float(i)/fact(dble(i))*(q/2.d0)**float(i)
         end do


        allocate(Zbig(ncontr,ncontr,ncontr,ncontr))

        open(40, file='Zcotr.dat', status='old', access='stream', form='unformatted')
        read(40) Zbig
        close(40)

        print*, Zbig(1,2,2,3)
        print*,ncontr,gs(:)
        print*,gf(:)
        do i=1,ncontr
            veccounts(gs(i):gf(i))=i

        end do
        newgto=CEILING(1.d0/24.d0*(ngto-2)*(ngto-1)*(3*ngto-1)*ngto)
        ncontr2=2
        count=1
        allocate(bigvec(4,newgto))
        do i=1,ngto
            do j=i+1,ngto
                do k=i+1,ngto
                    do r=k+1,ngto
                        bigvec(1,count)=i
                        bigvec(2,count)=j
                        bigvec(3,count)=k
                        bigvec(4,count)=r
                        count=count+1


                    enddo
                enddo
            enddo
        enddo

        !$OMP PARALLEL  do  num_threads(12) private(ii,jj,kk,rr,Zcontrred), &
        !$OMP& private(f,h,hx,hy,hz,i,j,k,r,iii,dx1,dx2,dy1,dy2,dz1,dz2) shared(bigvec,q,l,m,n,gf,gc,gs,p0matrix,Zbig), &
        !$OMP& shared( cutoffz,cutoffmd,BesselQ,BesselQ2) REDUCTION(+:tsi), &
        !$OMP& schedule(dynamic)

        do iii=1,newgto
            !print*,i,j,k,r
            !Sprint*,ii,jj,kk,rr,i,j,k,r
            i=bigvec(1,iii)
            j=bigvec(2,iii)
            k=bigvec(3,iii)
            r=bigvec(4,iii)
            ii=veccounts(i)
            jj=veccounts(j)
            kk=veccounts(k)
            rr=veccounts(r)
            Zcontrred=Zbig(ii,jj,kk,rr)
          !  if (Zcontrred>cutoffz) then
            hx = px(k, r) - px(i, j)
            hy = py(k, r) - py(i, j)
            hz = pz(k, r) - pz(i, j)
            h = (hx * hx + hy * hy + hz * hz)**0.5

            dx1=>dx(:,:,:,j,i)
            dy1=>dy(:,:,:,j,i)
            dz1=>dz(:,:,:,j,i)
            dx2=>dx(:,:,:,r,k)
            dy2=>dy(:,:,:,r,k)
            dz2=>dz(:,:,:,r,k)

            if (h < cutoffcentre) then

                call tot_integral_ijkr_pzero_read(nq,l,m,n, p0matrix, dx1,dy1,dz1,dx2,&
                        dy2,dz2,i, j, k, r, &
                        cutoffz, cutoffmd,f)


            else

                call tot_integral_k_ijkr_read(q,l,m,n,h,hx,hy,hz,dx1,dy1,dz1,dx2,&
                        dy2,dz2,i, j, k, r, &
                        BesselQ,BesselQ2,cutoffz, cutoffmd,f)





            end if
            tsi = tsi + 8.000_dp *f*zcontrred* e12(:,i,j)*e12(:,k,r)

        !    end if

        end do
        !$OMP END parallel DO

      print*,OMP_get_num_threads()
      print*,'intermediate step', tsi(1),count

        !$OMP PARALLEL do private(zcontrred,ii,jj,kk,rr), &
        !$OMP& private(f,h,hx,hy,hz,i,j,k,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n,gf,gc,gs,p0matrix,Zbig), &
        !$OMP& shared( cutoffz,cutoffmd,BesselQ,BesselQ2) REDUCTION(+:tsi), &
        !$OMP& schedule(dynamic)
        do i=1,ngto
            do j=i+1,ngto
                do r=i+1,ngto
                    hx = px(i, r) - px(i, j)
                    hy = py(i, r) - py(i, j)
                    hz = pz(i, r) - pz(i, j)
                    h = sqrt((hx * hx + hy * hy + hz * hz))
                    ii=veccounts(i)
                    jj=veccounts(j)
                    kk=veccounts(i)
                    rr=veccounts(r)
                    Zcontrred=Zbig(ii,jj,kk,rr)


                    dx1=>dx(:,:,:,j,i)
                    dy1=>dy(:,:,:,j,i)
                    dz1=>dz(:,:,:,j,i)
                    dx2=>dx(:,:,:,r,i)
                    dy2=>dy(:,:,:,r,i)
                    dz2=>dz(:,:,:,r,i)


                    if (h < cutoffcentre) then
                        call tot_integral_ijkr_pzero_read(nq, l,m,n, p0matrix, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, i, r, &
                                cutoffz, cutoffmd,f)
                    else

                        call tot_integral_k_ijkr_read(q, l,m,n,h,hx, hy, hz, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, j, i, r, &
                                besselq,BesselQ2,cutoffz, cutoffmd,f)



                    end if
                    tsi = tsi + 4.000 * f * zcontrred *e12(:, i, j) * e12(:, i, r)
                    count=count+1

                  !   deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

                end do
            end do

        end do
    !$OMP END parallel DO
         print*,'intermediate step', tsi(1)

       !$OMP PARALLEL do private(ii,kk,rr,Zcontrred), &
        !$OMP& private(f,h,hx,hy,hz,i,k,r,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
        !$OMP& shared( cutoffz,cutoffmd,BesselQ,BesselQ2) REDUCTION(+:tsi), &
        !$OMP& schedule(dynamic)
        do i=1,ngto
            do k=1,ngto
                do r=k+1,ngto
                    hx = px(k, r) - px(i, i)
                    hy = py(k, r) - py(i, i)
                    hz = pz(k, r) - pz(i, i)
                    h = sqrt((hx * hx + hy * hy + hz * hz))
                    ii=veccounts(i)

                    kk=veccounts(k)
                    rr=veccounts(r)
                    Zcontrred=Zbig(ii,ii,kk,rr)


                    dx1=>dx(:,:,:,i,i)
                    dy1=>dy(:,:,:,i,i)
                    dz1=>dz(:,:,:,i,i)
                    dx2=>dx(:,:,:,r,k)
                    dy2=>dy(:,:,:,r,k)
                    dz2=>dz(:,:,:,r,k)
!                    zcontrred=zcontrred
!                    zcontrred2=zcontrred2

                    if (h < cutoffcentre) then
                        call tot_integral_ijkr_pzero_read(nq,l,m,n, p0matrix, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, k, r, &
                                cutoffz, cutoffmd, f)
                    else

                        call tot_integral_k_ijkr_read(q, l,m,n,h, hx, hy, hz, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, k, r, &
                                besselq,BesselQ2,cutoffz, cutoffmd, f)



                    end if
                    tsi = tsi+ 4.000 * Zcontrred* f * e12(:, i, i) * e12(:, k, r)
                    count=count+1
               !     deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

                end do
            end do
        end do

        !$OMP END parallel DO
        print*,'intermediate step', tsi(1)
!
       !$OMP PARALLEL do private(ii,kk,zcontrred), &
        !$OMP& private(f,h,hx,hy,hz,i,k,dx1,dx2,dy1,dy2,dz1,dz2) shared(q,l,m,n, p0matrix), &
          !$OMP& shared( cutoffz, cutoffmd,BesselQ,BesselQ2) REDUCTION(+:tsi), &
              !$OMP& schedule(dynamic)
        do i=1,ngto
            do k=i+1,ngto

                hx = px(k, k) - px(i, i)
                hy = py(k, k) - py(i, i)
                hz = pz(k, k) - pz(i, i)
                h = sqrt((hx * hx + hy * hy + hz * hz))

                 ii=veccounts(i)

                 kk=veccounts(k)

                 Zcontrred=Zbig(ii,ii,kk,kk)


                dx1=>dx(:,:,:,i,i)
                dy1=>dy(:,:,:,i,i)
                dz1=>dz(:,:,:,i,i)
                dx2=>dx(:,:,:,k,k)
                dy2=>dy(:,:,:,k,k)
                dz2=>dz(:,:,:,k,k)


                if (h < cutoffcentre) then
                    call tot_integral_ijkr_pzero_read(nq, l,m,n, p0matrix, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, k, k, &
                                cutoffz, cutoffmd, f)
                else

                    call tot_integral_k_ijkr_read(q,l,m,n,h, hx, hy, hz, dx1,dy1,dz1,dx2,&
                                dy2,dz2,i, i, k, k, &
                                 besselq,BesselQ2,cutoffz, cutoffmd, f)
                end if
                tsi = tsi+ 2.000 * f *Zcontrred* e12(:, i, i) * e12(:, k, k)

               ! deallocate(dx1red, dy1red,dz1red,dx2red,dy2red,dz2red)

            end do
        end do


        !$OMP END parallel DO
        print*,tsi(1)
        !$OMP PARALLEL do private(ii,Zcontrred), &
        !$OMP& private(f,h,hx,hy,hz,i,dx1,dy1,dz1) shared(q,ll, p0matrix), &
          !$OMP& shared( cutoffz, cutoffmd,BesselQ,BesselQ2) REDUCTION(+:tsi), &
            !$OMP& schedule(dynamic)
        do i=1,ngto

             ii=veccounts(i)
             Zcontrred=Zbig(ii,ii,ii,ii)

            dx1=>dx(:,:,:,i,i)
            dy1=>dy(:,:,:,i,i)
            dz1=>dz(:,:,:,i,i)

!            zcontrred=zcontrred/8.0
!            zcontrred2=zcontrred2/8.0

            call tot_integral_ijkr_pzero_read(nq, l,m,n,p0matrix, dx1,dy1,dz1,dx1,&
                                dy1,dz1,i, i, i, i, &
                                 cutoffz, cutoffmd,f)

            tsi = tsi + Zcontrred*f * e12(:, i, i) * e12(:, i, i)
            count=count+1



        end do
        !$OMP END parallel DO
        print*,tsi(1)

    end subroutine tot_integration_contract

real(kind=SELECTED_REAL_KIND(15))  recursive function fact(n) result (ans)
implicit none
real(kind=SELECTED_REAL_KIND(15)), intent(in) :: n
    if(n == 0) then
        ans=1
        return
    end if

    if (n == 1) then
        ans = 1

    else

        ans = n * fact(n-1)
    end if

    return
end function fact




End Module TSj0contr