program Zcontrred
    use omp_lib
    implicit none
     INTEGER, PARAMETER :: dpp = SELECTED_REAL_KIND(kind(1.d0))
     INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
     INTEGER, PARAMETER :: ikind = SELECTED_int_KIND(8)
     integer(kind=ikind), dimension(:), allocatable:: m1,m2,m3,m4
     real(kind=dp), dimension(:), allocatable::totalfin
     real(kind=dp), allocatable,dimension(:,:,:):: z1,z2
     real(kind=dpp), allocatable, dimension(:,:,:,:)::zcontr
     real(kind=dp), allocatable, dimension(:,:):: mmod
     integer(kind=ikind)::i,j,k,r,nmat,norbs,ngtos,ncap,ic,jc,kc,rc



    open(unit=15, file='MOs.dat')
    read(15,*)norbs,ngtos
    allocate(mmod(norbs,ngtos))
    print*,norbs,ngtos
    do i=1,norbs
        read(15,*) (mmod(i,j), j=1, ngtos)
    end do
    close(15)
    nmat=1324

    allocate(m1(nmat), m2(nmat),m3(nmat), m4(nmat), totalfin(nmat))
    open(unit=15,file='twordm_fortran.dat')
            do i=1,1324
            read(15,*)m1(i),m2(i),m3(i),m4(i), totalfin(i)
            enddo
            close(15)


    allocate(z1(nmat,ngtos,ngtos), z2(nmat,ngtos,ngtos))
    ncap=ngtos
    ic=0
    jc=0
    kc=0
    rc=0




    ngtos=240
     do  i=1,ngtos
            do j=1,ngtos


                        z1(:,i, j) = totalfin * (mmod(m1, i) * mmod(m2, j) + mmod(m1, j) * mmod(m2, i))
                        z2(:, i, j)= mmod(m3, i) * mmod(m4, j) + mmod(m3, j) * mmod(m4, i)
            enddo
        enddo

     print*,'allocating the megabig matrix'

     allocate(zcontr(ngtos,ngtos,ngtos,ngtos))

     zcontr=0.0_dpp
     print*,'can we populate it?'
     do i=1,ngtos
         do j=i+1,ngtos

              call dgemm('t','n', ngtos, ngtos, nmat, 1.0_dp/8.0_dp, z1(:,:,i), &
                              &           nmat, z2(:,:,j), nmat, 0.0_dp, zcontr(:,:,j,i), ngtos)

         end do
    end do






end program Zcontrred