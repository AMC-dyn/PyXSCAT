
Module Zcontr
    use Reader
    use onerdm
    use iso_fortran_env, only:  int8, int16, int32, int64
    contains
    subroutine Contraction(prueba2)
    implicit none
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
    character(len=60):: conf1, file_out,var,var2rdm,file_in
    character(len=100)::file_bit
    real(kind=dp), dimension(:), allocatable:: atoms
    integer(kind=ikind):: ngtos,norbs,nconfs,nstates,state1,state2,natoms,ncontr,norbs2,ngtos2
    integer(kind=int64)::numberlines
    integer(kind=ikind),dimension(:,:), allocatable:: confs

    real(kind=dp),dimension(:), allocatable:: q
    real(kind=dp),dimension(:), allocatable:: coeffs,total
    real(kind=dp),dimension(:,:), allocatable:: civs
    real(kind=dp),dimension(:,:), allocatable:: mmod
    real(kind=dp),dimension(:,:), allocatable:: geom
    real(kind=dp),dimension(:), allocatable:: ga, xx,yy,zz
    real(kind=dp), dimension(:), allocatable:: result,q_abs
    integer*8, dimension(1):: nnn,start1,end1, nnn2, start_2, end_2,ordering1,ordering2
    integer*8, dimension(1,1):: newdat
    real(kind=dp):: cutoffcentre,cutoffz,cutoffmd,qmin,qmax,time1,time2,time3
    integer(kind=ikind),dimension(:), allocatable:: l,m,n,group,gs,gf,gc,contrvec
    integer(kind=ikind):: typec, i, j,k,ll, npoints,ncivs,lconfs,maxl,ng,nq,count
    logical:: jeremyR, mcci, hf,molpro,molcas,bagel,bitwise,fci
    integer(kind = int64), dimension(:, :), allocatable :: mat
    real(kind=dp), dimension(:,:,:,:), allocatable :: twordm,zcontract,prueba1
    real(kind=dp), dimension(:,:,:,:), allocatable, intent(out) :: prueba2
    real(kind=dp), dimension(:,:), allocatable:: interm,interm2,mos




    call read_files(nconfs,ngtos,norbs,ng,ncontr,state1,state2,natoms,typec,maxl,npoints, &
               cutoffcentre,cutoffz,cutoffmd,atoms,coeffs,xx,yy,zz,ga,l,m,n,group,mmod,geom, &
               jeremyR,mcci,hf,molpro,molcas,bagel,qmin,qmax,&
               file_out,file_bit,var,gs,gc,gf,confs,civs,lconfs,ncivs,contrvec,bitwise)

            open(unit=15, file='bitwise.dat')
              numberlines=0
              do while(.true.)
                  read (15, *, end=999) i
                  numberlines=numberlines+1
              end do


999 continue
        close(15)

    file_in='bitwise.dat'
    file_out='es.dat'
    call mcci_to_bit(file_in,file_out,numberlines)
    call createtwordm_bit(file_out,numberlines,mat,total)
    !norbs=maxval(mat(:,1))
    allocate(twordm(norbs,norbs,norbs,norbs), zcontract(ncontr,ncontr,ncontr,ncontr))
    open(unit=15, file='MOs2.dat')
    read(15,*)norbs2,ngtos2
    allocate(mos(norbs2,ngtos2))

    do i=1,norbs2
        read(15,*) (mos(i,j), j=1, ngtos2)
    end do
    close(15)
    !print*,shape(mos),norbs,norbs2
    !print*,mos(1,:)
    twordm=0.0_dp
    do i=1,size(total)

        twordm(mat(i,1), mat(i,2), mat(i,3),mat(i,4))=total(i)

    end do
     call cpu_time(time1)
    !print*,'twordm', twordm(1,2,3,4)

    allocate(prueba1(ncontr,norbs,norbs,norbs))

    call dgemm('n', 'n', ncontr, norbs*norbs*norbs, norbs, 1.0_dp ,mos , &
                                        &           ncontr,twordm, norbs, 0.0_dp, prueba1, ncontr)

    !print*,'first mult', prueba1(2,3,3,7)
    allocate(prueba2(norbs,ncontr,norbs,norbs))

    !print*,'reshaping '

   ! prueba2=reshape(prueba1,(/norbs,ncontr,norbs,norbs/), order=[2,1,3,4])
     forall(ll=1:norbs,k=1:norbs,j=1:ncontr,i=1:norbs) &
      prueba2(i,j,k,ll)=prueba1(j,i,k,ll)
    !print*,'reallocating'

    deallocate(prueba1)
    allocate(prueba1(ncontr,ncontr,norbs,norbs))
    !print*,'entering 2nd dgemm'
    prueba1=0.0_dp
    call dgemm('n', 'n', ncontr, norbs*norbs*ncontr, norbs, 1.0_dp ,mos , &
                                        &           ncontr,prueba2, norbs, 0.0_dp, prueba1, ncontr)

    deallocate(prueba2)

    allocate(prueba2(norbs,ncontr,ncontr,norbs))
    !print*,'second mult', prueba1(3,2,3,7)
    !print*,size(prueba1)

    do ll=1,norbs
        do k=1,ncontr
            do j=1,ncontr
                do i=1,norbs
                    prueba2(i,j,k,ll)=prueba1(k,j,i,ll)
                end do
            end do
        end do

     end do
    deallocate(prueba1)

    allocate(prueba1(ncontr,ncontr,ncontr,norbs))
    !print*,'entering 3rd dgemm'
    call dgemm('n', 'n', ncontr,norbs*ncontr*ncontr, norbs, 1.0_dp ,mos , &
                                        &           ncontr,prueba2, norbs, 0.0_dp, prueba1, ncontr)



    !print*,'third mult', prueba1(3,2,3,7)
    deallocate(prueba2)
    allocate(prueba2(norbs,ncontr,ncontr,ncontr))

     do ll=1,ncontr
        do k=1,ncontr
            do j=1,ncontr
                do i=1,norbs
                    prueba2(i,j,k,ll)=prueba1(ll,j,k,i)

                end do
            end do
        end do

     end do
    !print*,'trying shape', shape(prueba2)
    deallocate(prueba1)
    allocate(prueba1(ncontr,ncontr,ncontr,ncontr))
    call dgemm('n', 'n', ncontr,ncontr*ncontr*ncontr, norbs, 1.0_dp ,mos , &
                                        &           ncontr,prueba2, norbs, 0.0_dp, prueba1, ncontr)

    !print*,prueba1(1,2,2,3)
    deallocate(prueba2)
    !print*,'before last reshape'
      !print*,shape(prueba1)
    !print*,'before last reshape 2'

 !print*,'fourth mult', prueba1(3,3,7,2)
 allocate(prueba2,mold=prueba1)
    ! !dir$ ivdep
    do ll=1,ncontr
      !   !dir$ ivdep
        do k=1,ncontr
     !        !dir$ ivdep
            do j=1,ncontr
      !           !dir$ ivdep
                do i=1,ncontr
                    prueba2(i,j,k,ll)=prueba1(j,k,ll,i)

                end do
            end do
        end do
    end do
   !print*,'fourth mult', prueba2(2,3,3,7)
   prueba1=prueba2
   prueba2=0.0
    !print*,'here we are'

    call cpu_time(time2)
    !print*,time2-time1

    !$OMP PARALLEL DO private(ll,k,j,i) shared(prueba1) schedule(dynamic)
      do ll=1,ncontr
          !dir$ ivdep
        do k=1,ncontr
            !dir$ ivdep
            do j=1,ncontr
                !dir$ ivdep
                do i=1,ncontr
              prueba2(i,j,k,ll)=prueba1(i,j,k,ll)+prueba1(j,i,k,ll)+prueba1(i,j,ll,k)+prueba1(j,i,ll,k)+ &
              prueba1(k,ll,i,j)+prueba1(k,ll,j,i)+prueba1(ll,k,i,j)+prueba1(ll,k,j,i)
                    end do
                    end do
                    end do
              end do
      !$OMP END parallel DO
      prueba2=prueba2/8.00_dp
      call cpu_time(time3)
    !print*,time3-time2
end subroutine

End Module Zcontr