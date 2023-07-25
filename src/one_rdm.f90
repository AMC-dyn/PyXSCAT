
   module onerdm
       use twordms
       contains

        function excitations(det1,det2,Nint) result(n_excitations)
           implicit none
            INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
           integer(kind=ikind), intent(in) :: det1(Nint,2), det2(Nint,2)
           integer(kind=ikind) , intent(in) :: Nint
           integer(kind=ikind) :: n_excitations

           integer(kind=ikind) :: l

           n_excitations = &
                   popcnt(xor( det1(1,1), det2(1,1)) ) + &
                           popcnt(xor( det1(1,2), det2(1,2)) )

           do l=2,Nint
               n_excitations = n_excitations + &
                       popcnt(xor( det1(l,1), det2(l,1)) ) + &
                       popcnt(xor( det1(l,2), det2(l,2)) )
           end do
           n_excitations = ishft(n_excitations,-1)

       end function

     subroutine onerdm_creat(confs,civs,onerdm,maxnmo,state1,state2)
            INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
         integer(kind=ikind), intent(in), dimension(:,:) :: confs
         real(kind=dp), intent(in), dimension(:,:) :: civs
         real (kind=dp), intent(out), dimension(:,:), allocatable :: onerdm
         integer(kind=ikind),dimension(:,:), allocatable :: ep2, ndiff
         integer(kind=ikind), intent(out) :: maxnmo
            integer(kind=ikind),intent(in):: state1,state2



         call maxcoincidence(confs, ep2,ndiff)

         call createonerdm(confs,civs,ndiff,ep2,onerdm,maxnmo,state1,state2)



         end subroutine


!     subroutine maxcoincidence(confs, ep2,ndiff)
!        use types
!        implicit none
!        integer(kind=ikind), intent(in), dimension(:,:) :: confs
!        integer(kind=ikind), intent(out), dimension(:,:), allocatable :: ep2, ndiff
!        integer(kind=ikind), dimension(size(confs(:,1)), size(confs(1,:))):: matdum
!        integer(kind=ikind), dimension(:,:), allocatable :: mat1
!        integer(kind=ikind) :: i,j, c1,c2, count
!
!
!        allocate(ep2(size(confs(:,1)),size(confs(:,1))))
!        ep2=1
!        count=0
!        matdum=0
!        print*,'holaa'
!            do i=1,size(confs(:,1))
!                count=0
!                do j=1,size(confs(1,:))
!
!
!                    if (confs(i,j)/=0) then
!                        count=count+1
!                        matdum(i,count)=j
!                        end if
!                end do
!            end do
!           print*, 'matdum constructed'
!            allocate(mat1(size(confs(:,1)),count), ndiff(size(confs(:,1)),size(confs(:,1))))
!            ndiff=0
!            mat1=matdum(:,1:count)
!            print*, 'mat1 constructed'
!            do c1=1, size(confs(:,1))
!                do c2=c1+1,size(confs(:,1))
!                    do i=1,size(mat1(1,:))
!                        if  (mat1(c1,i) /= mat1(c2,i)) then
!
!
!                            do j=1,size(mat1(1,:))
!                                if (mat1(c1,i) /= mat1(c2,j)) then
!                                    ep2(c1,c2)=-ep2(c1,c2)
!                                    ep2(c2,c1)=ep2(c1,c2)
!
!                                end if
!                            end do
!                        end if
!                    end do
!
!                    do j=1,size(confs(1,:))
!                        if (confs(c1,j)/=confs(c2,j)) then
!                           ndiff(c1,c2)= ndiff(c1,c2) +1
!                           ndiff(c2,c1)= ndiff(c2,c1) +1
!                        end if
!                    end do
!
!                end do
!            end do
!
!    print*,ndiff(1,2)
!    end subroutine maxcoincidence


    subroutine createonerdm(confs,civs,ndiff,ep2,onerdm,maxnmo,state1,state2)

    implicit none

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

    integer(kind=ikind), intent(in),dimension(:,:) :: confs
    integer(kind=ikind), intent(in), dimension(:,:) :: ep2, ndiff
    real(kind=dp), intent(in), dimension(:,:) :: civs
    integer(kind=ikind), intent(in) :: state1,state2


    real(kind=dp), dimension(:,:), intent(out),allocatable :: onerdm
    integer(kind=ikind), intent(out) :: maxnmo

    integer(kind=ikind) :: ep,nc1,lconfs,norbs,nc2,sorb,rorb,qorb,porb,p,q,r,s,c1,c2,count1,count

    integer(kind=ikind) :: i,i1,n,count2,eg, ndiff1


    integer(kind=ikind):: diffs1,diffs2
    integer(kind=ikind) :: spin1,spin2

    lconfs=size(confs(1,:))
    maxnmo=lconfs/2
    nc1=1
    nc2=1
    ep=1
    allocate(onerdm(maxnmo,maxnmo))
    onerdm=0.0_dp
    do c1=1,size(confs(:,1))
        do c2=1,size(confs(:,1))

            ndiff1=ndiff(c1,c2)
            ep=ep2(c1,c2)


            if (ndiff1==2) then
                if (c1<=20) print*,'conguration of interest',c1,c2,ep

                do n=1,size(confs(c1,:))
                    if (confs(c1,n) /= confs(c2,n)) then
                        if (confs(c1,n) /= 0) THEN
                            diffs1=nint((n) / 2.0_dp + 0.1)
                            spin1=confs(c1,n)



                        elseif (confs(c2,n) /= 0) THEN
                            diffs2=nint((n) / 2.0_dp + 0.1)
                            spin2=confs(c2,n)
                        end if
                    end if
                enddo


                    sorb = diffs1
                    qorb = diffs2


                    eg = 1.00


                        onerdm(sorb , qorb ) = onerdm(sorb , qorb ) + &
                                civs(c1,state1) * civs(c2,state2) * ep

!                        print*,'conf1', pack(confs(c1,:),confs(c1,:)/=0)
!                        print*,'conf2',pack(confs(c2,:),confs(c2,:)/=0)
!                        print*,'phase', c1,c2,ep



            elseif (ndiff1== 0) then

                ep = 1

                do i1=1,size(confs(c1,:))
                    if (confs(c1,i1)/=0) then

                        sorb=nint((i1) / 2.0_dp + 0.1)

                        onerdm( sorb , sorb)  = onerdm(sorb , sorb) + civs(c1,state1) * civs(c1,state2) * ep


                    end if
                    enddo

            end if


        end do

    end do
    count=0.0
    open(file='onerdm.dat', unit=15)
    do i=1,size(onerdm(:,1))

         write(15,'(1000F14.7)')( onerdm(i,c1) ,c1=1,size(onerdm(i,:)))
    end do
    close(15)

     open(file='configurations.dat', unit=15)
    do i=1,size(confs(:,1))

         write(15,*)( confs(i,c1) ,c1=1,size(confs(i,:)))
    end do
    close(15)

    print*, 'onerdm calculated'



    end subroutine


     subroutine one_rdm_slow(file_read,onerdm,maxnmo)

          Use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

    character(len=20),intent(in) :: file_read



    real(kind=dp), dimension(:,:), intent(out),allocatable :: onerdm
    integer(kind=ikind), intent(inout) :: maxnmo

    integer(kind=ikind) :: ep,nc1,lconfs,norbs,nc2,sorb,rorb,qorb,porb,n1a,n1b,error_1,error_2
          integer(kind=ikind),allocatable, dimension(:) ::Nalpha,Nbeta
          real(kind=dp), dimension(:), allocatable :: civs
          real(kind=dp):: civ1,civ2

    integer(kind=ikind) :: i,i1,n,count2,eg, ndiff1,n2a,n2b,j,count,c1,c2
    integer(kind=ikind), dimension(:), allocatable :: c1a,c1b,c2a,c2b





    nc1=1
    nc2=1
    ep=1
    allocate(onerdm(maxnmo,maxnmo))
    onerdm=0.0_dp
    open(15,file=file_read)

    count=0
    do while (error_1 /= -1)
        read(15,*, iostat = error_1) j
        count=count+1

    enddo
    print*,count
    close(15)
      allocate(civs(count-1), Nalpha(count-1), Nbeta(count-1))
          civs=0
          Nalpha=0
          Nbeta=0
        open(15,file=file_read)
        do i=1,count-1
            read(15,*)j, civs(i), Nalpha(i), Nbeta(i)
        enddo
        close(15)
count=10
         do c1=1,count-1

             do c2=1,count-1

                 N1a=Nalpha(c1)
                 N2a=Nalpha(c2)

                 N1b=Nbeta(c1)
                 N2b=Nbeta(c2)
                 civ1=civs(c1)
                 civ2=civs(c2)


            ep=1

            if (N1a==N2a .and. N1b==N2b) then

            !case where everything is the same
                ep=1
                c1a=integer2binary_orbs(N1a)

                c1b=integer2binary_orbs(N1b)

                do i=1,size(c1a)
                   sorb=c1a(i)

                   onerdm( sorb , sorb)  = onerdm(sorb , sorb) + civ1 * civ2
                end do

                 do i=1,size(c1b)
                   sorb=c1b(i)

                   onerdm( sorb , sorb)  = onerdm(sorb , sorb) + civ1 * civ2

                end do




            elseif (N1a==N2a) then
                ndiff1=0
                c1b=integer2binary_orbs(N1b)
                c2b=integer2binary_orbs(N2b)


                do i=1,size(c1b)
                    if (c1b(i)/=c2b(i)) then
                        sorb = c1b(i)
                        qorb = c2b(i)

                        do j=1,size(c2b)
                            if (c1b(i)/=c2b(j)) then
                                ep=-ep

                            end if
                        end do
                         print*,c1,c2,ep
                        onerdm(sorb , qorb ) = onerdm(sorb , qorb ) + &
                                civ1*civ2*ep



                    end if
                end do


            elseif (N1b==N2b) then
                     ndiff1=0
                c1a=integer2binary_orbs(N1a)
                c2a=integer2binary_orbs(N2a)


                do i=1,size(c1a)
                    if (c1a(i)/=c2a(i)) then
                        sorb = c1a(i)
                        qorb = c2a(i)

                        do j=1,size(c2a)
                            if (c1a(i)/=c2a(j)) then
                                ep=-ep

                            end if
                        end do
                        print*,c1,c2,ep
                        onerdm(sorb , qorb ) = onerdm(sorb , qorb ) + &
                                civ1*civ2*ep



                    end if
                end do
                end if

            end do

enddo



    open(file='onerdm.dat', unit=15)
    do i=1,size(onerdm(:,1))

         write(15,'(1000F14.7)')( onerdm(i,c1) ,c1=1,size(onerdm(i,:)))
    end do
    close(15)


    print*, 'onerdm calculated'



    end subroutine

       subroutine createtwordm_bit_fci(file_read,numberlines,newdat,irep,start,end,mat,total)


     Use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
     character(len=20), intent(in) :: file_read
     integer(kind=ikind) :: norbs
     integer(KIND=IKIND), intent(in) :: numberlines
     integer(kind=ikind), intent(in), dimension(:,:):: newdat
     integer(kind=ikind),intent(in),dimension(:):: irep,start,end
     double precision, parameter :: phase_dbl(0:1)=(/1.d0,-1.d0/)
    real(kind=dp), intent(out), dimension(:), allocatable :: total
    integer(kind=ikind), intent(out), dimension(:,:), allocatable :: mat

    real(kind=dp), dimension(:,:,:,:), allocatable :: twordm
    real(kind=dp), dimension(:,:), allocatable :: onerdm
    integer(kind=ikind) :: ep, i, j, k, l, compnum, c1
     integer(kind=ikind):: X_vec(numberlines)
     integer(kind=ikind), dimension(:,:,:), allocatable :: nalphbet
     real(kind=dp), dimension(:), allocatable :: civs
    real(kind=dp) :: time1, time2, cutoff
     double precision:: civ1,civ2,c,c_term,phase
    integer(kind=ikind) :: buffer, buffer_2, hole_int, particle_int, buffer_p, index, eg,index1
    integer(kind=ikind) :: porb,rorb,qorb,sorb, count2, p, q, r, s, dummy, n_single

    logical(4), dimension(:), allocatable :: logicaltwordms
    integer(kind=ikind),  dimension(:,:), allocatable :: matdum
    real(kind=dp), dimension(:), allocatable :: totaldum
    integer(kind=ikind), dimension(:,:), allocatable :: single_exc,hole_exc, particle_exc,phases,counter_exc,starts,ends
    integer(kind=ikind), dimension(:,:), allocatable :: double_exc,hole_exc_1,hole_exc_2
integer(kind=ikind), dimension(:,:),allocatable :: particle_exc_1,particle_exc_2, phases_2, counter_exc_2,starts_2,ends_2
    integer(kind=ikind), dimension(:,:,:), allocatable :: beta_excs
    integer(kind=ikind) :: buffer_prime, buffer_prime_2, tz,tz2, count_p,n,ss,ee,spin11,spin22,buffer2,buffer_3,buffer_4,newdat_1,len_rm

    integer(kind=ikind) :: n_double,count,nel,index2,myhigh,mylow,myhigh1,myhigh2,mylow1,mylow2,count_total
    integer(kind=ikind),dimension(:), allocatable :: find_ind,find_ind_2,all_indexes,occ,unocc
      double precision, external      :: ddot

    allocate(nalphbet(1,2,numberlines), civs(numberlines), all_indexes(numberlines))
    open(15,file=file_read)
!
     count_total=0
     civs=0
     Nalphbet=0
     all_indexes=0
     do i=1,numberlines
         all_indexes(i)=i
     end do

!
     do i=1,numberlines
            read(15,*)dummy, civs(i), Nalphbet(1,1,i), Nalphbet(1,2,i)

     enddo
        close(15)
!
     norbs=maxval(integer2binary_orbs(maxval(Nalphbet(1,1,:))))



     print*,'Number of orbitals',norbs,numberlines
     print*, popcnt(xor(3187,3159))
!    print*, 'read file'
    allocate(twordm(norbs,norbs,norbs,norbs))
!     print*, 'allocataed twordm'

    do i=1,norbs
        do j=1,norbs
            do k=1,norbs
                do l=1,norbs
                  twordm(i,j,k,l)=0.0_dp
                end do
            end do
        end do
     end do

!    nc1=1
!    nc2=1
    ep=1

      compnum=0
     do i=1,norbs
         compnum=compnum+2**(i-1)
     end do


     n_single=popcnt(xor(newdat(i,1), compnum))*popcnt(newdat(1,1))


     allocate(single_exc(size(newdat(:,1)), n_single))
     allocate(hole_exc(size(newdat(:,1)), n_single), particle_exc(size(newdat(:,1)), n_single),phases(size(newdat(:,1)), n_single))
     allocate(counter_exc(size(newdat(:,1)), n_single))
     allocate(starts(size(newdat(:,1)), n_single), ends(size(newdat(:,1)), n_single))

    ! Calculate all single and double excitations for the unique determinants
     call cpu_time(time1)
     single_exc=0
     particle_exc=0
     hole_exc=0
     counter_exc=0
     phases=0
     starts=0
     ends=0
     print*,size(irep),size(newdat(:,1))

     do i=1,size(newdat(:,1))

               if (i==1000) then
                    call cpu_time(time2)
                    print*,time2-time1
               end if
               count2=1
               buffer=xor(newdat(i,1), compnum)
               do while (buffer /= 0_8)
                   j = trailz(buffer)
                   hole_int=ibset(newdat(i,1),j)
                   buffer_2=newdat(i,1)
                   do while(buffer_2/=0_8)
                      k=trailz(buffer_2)
                      particle_int=ibclr(hole_int,k)
                      if (particle_int==newdat(i,1)) cycle
                      buffer_2 = iand(buffer_2,buffer_2-1_8)

                      index=0
                      do l=1,size(newdat(:,1))
                          if (particle_int==newdat(l,1)) then
                              index=l
                              exit
                          end if
                      end do
                      if (index==0) cycle

                      count_p=irep(index)
                      if (count_p==0) cycle

                      single_exc(i,count2)=particle_int
                      hole_exc(i,count2)=k+1
                      particle_exc(i,count2)=j+1
                      counter_exc(i,count2)=count_p
                      starts(i,count2)=start(index)
                      ends(i,count2)=end(index)

                      ep=1
                      mylow=min(j+1,k+1)
                      myhigh=max(j+1,k+1)
                      ep=POPCNT(IAND(newdat(i,1),IAND(ibset(0_8,myhigh-1_8)-1_8,ibclr(-1,mylow)+1_8)))
!                      buffer_prime=newdat(i,1)
!                      buffer_prime_2=particle_int
!
!                        do while (buffer_prime/=0_8)
!                                tz=trailz(buffer_prime)
!                                tz2=trailz(buffer_prime_2)
!
!                                if (tz/=tz2 .and. btest(Nalphbet(1,1,index),tz)) then
!
!                                    ep=-ep
!                                end if
!                                buffer_prime=iand(buffer_prime, buffer_prime-1_8)
!                                buffer_prime_2=iand(buffer_prime_2, buffer_2-1_8)
!                        end do
                        phases(i,count2)= ep

                      count2=count2+1
                   end do
                    buffer = iand(buffer,buffer-1_8)
               end do
     end do
     call cpu_time(time2)
     print*,'created single excitations', time2-time1, count2, n_single,compnum
     !print*, newdat(2,1), single_exc(2,:)
     call cpu_time(time1)
     nel=popcnt(newdat(1,1))
     n_double=combo(nel,2)*combo(norbs-nel,2)

     allocate(occ(nel), unocc(norbs-nel))
     allocate(double_exc(size(newdat(:,1)), n_double), phases_2(size(newdat(:,1)), n_double))
     allocate(hole_exc_1(size(newdat(:,1)), n_double), particle_exc_1(size(newdat(:,1)), n_double))
     allocate(hole_exc_2(size(newdat(:,1)), n_double), particle_exc_2(size(newdat(:,1)), n_double))
     allocate(counter_exc_2(size(newdat(:,1)), n_double),starts_2(size(newdat(:,1)), n_double),ends_2(size(newdat(:,1)), n_double) )

     double_exc=0
     phases_2=0
     hole_exc_1=0
     hole_exc_2=0
     particle_exc_1=0
     particle_exc_2=0
     counter_exc_2=0
     starts_2=0
     ends_2=0
     occ=0
     unocc=0


      do n=1,size(newdat(:,1))
            count2=1
            count=1
            buffer=newdat(n,1)

            do while(buffer/=0_8)
                 i=trailz(buffer)
                 occ(count)=i
                 count=count+1
                 buffer=iand(buffer,buffer-1_8)
            end do
           buffer=xor(newdat(n,1), compnum)
           count=1
           do while(buffer/=0_8)
                 i=trailz(buffer)
                 unocc(count)=i
                 count=count+1
                 buffer=iand(buffer,buffer-1_8)
           end do

           do i=1,size(occ)-1
                do j=i+1,size(occ)
                    do k=1,size(unocc)-1
                        do l=k+1,size(unocc)
                             buffer=ibclr(newdat(n,1),occ(i))

                             buffer=ibclr(buffer,occ(j))

                             buffer=ibset(buffer,unocc(k))

                             buffer=ibset(buffer,unocc(l))
                             if (buffer==newdat(n,1)) cycle
                             index=findloc(newdat(:,1),buffer,1)
                             if (popcnt(xor(buffer, newdat(n,1)))/=4) then
                                 print*,popcnt(xor(buffer, newdat(n,1)))
                             end if
                             if (index==0) cycle

                             double_exc(n,count2)=buffer
                             hole_exc_1(n,count2)=occ(i)+1
                             hole_exc_2(n,count2)=occ(j)+1
                             particle_exc_1(n,count2)=unocc(k)+1
                             particle_exc_2(n,count2)=unocc(l)+1

                             count_p=irep(index)
                             if (count_p==0) cycle

                             counter_exc_2(n,count2)=count_p
                             starts_2(n,count2)=start(index)
                             ends_2(n,count2)=end(index)

                             ep=0
                             mylow1=min(occ(i)+1,unocc(k)+1)
                             myhigh1=max(occ(i)+1,unocc(k)+1)
                             ep=POPCNT(IAND(newdat(n,1),IAND(ibset(0,myhigh1-1)-1,ibclr(-1,mylow1)+1)))
                             mylow2=min(occ(j)+1,unocc(l)+1)
                             myhigh2=max(occ(j)+1,unocc(l)+1)
                             ep=ep+POPCNT(IAND(newdat(n,1),IAND(ibset(0,myhigh2-1)-1,ibclr(-1,mylow2)+1)))

                             if((mylow2>mylow1).AND.(mylow2<myhigh1).AND.(myhigh2>myhigh1)) ep=ep+1

                             phases_2(n,count2)=ep

                             count2=count2+1

                        enddo
                    end do
                end do
           end do


          enddo

    call cpu_time(time2)
     print*,'created double excitations', time2-time1, count2, n

     call cpu_time(time1)
     print*,sum(newdat(:,1)-newdat(:,2))




do c1=1,numberlines



         if (c1==100000) then

                call cpu_time(time2)
                print*,time2-time1

         end if
        if (c1==700000) then

                call cpu_time(time2)
                print*,time2-time1

         end if

         civ1=civs(c1)


        c=civ1*civ1


        do spin11=1,2
            buffer=Nalphbet(1,spin11,c1)
            do while(buffer/=0_8)
                sorb=trailz(buffer)+1
                buffer=IAND(buffer,buffer-1_8)
                do spin22=1,2
                    buffer2=Nalphbet(1,spin22,c1)
                    do while(buffer2/=0_8)
                        qorb=trailz(buffer2)+1
                        buffer2=IAND(buffer2,buffer2-1_8)
                        if((qorb==sorb).AND.(spin11==spin22)) cycle
                        if(spin11==spin22) THEN
                            twordm(sorb,qorb,sorb,qorb)=twordm(sorb,qorb,sorb,qorb)-c !porb=sorb rorb=qorb

                        end if
                        twordm(qorb,sorb,sorb,qorb)=twordm(qorb,sorb,sorb,qorb)+c !porb=sorb rorb=qorb


                    end do
                end do


            enddo
        end do




         index1=findloc(newdat(:,1), Nalphbet(1,1,c1),1)

         index2=findloc(newdat(:,1), Nalphbet(1,2,c1),1)

         do i=1,size(single_exc(index1,:))
                index=0
                civ2=0.0_dp
                particle_int=single_exc(index1,i)
                if (particle_int==0) cycle
                count2=counter_exc(index1,i)
                if (count2==0) cycle
                ss = starts(index1,i)
                ee = ends(index1,i)


                index=findloc(Nalphbet(1,2,ss:ee),nalphbet(1,2,c1), 1)

                if (index>0) then

                index=index+ss-1
                civ2=civs(index)
                porb=hole_exc(index1,i)
                qorb=particle_exc(index1,i)
                ep=phases(index1,i)

                buffer=IBCLR(Nalphbet(1,1,c1),porb-1)
                !if (particle_int/=nalphbet(1,2,c1) ) then
                c_term=phase_dbl(iand(ep,1))*2.00*civ1*civ2
              !  else
                !     c_term=phase_dbl(iand(ep,1))*2.00*civ1*civ2
              !  end if

                do while(buffer /= 0_8)
                    rorb=trailz(buffer)+1
                    buffer=IAND(buffer,buffer-1_8)

                    twordm(porb,rorb,rorb,qorb)=twordm(porb,rorb,rorb,qorb)+c_term

                    twordm(rorb,porb,qorb,rorb)=twordm(rorb,porb,qorb,rorb)+c_term !sorb=rorb
                    !case 2 spins and spinr are the same by construction as are spinp and spinq
                    twordm(rorb,porb,rorb,qorb)=twordm(rorb,porb,rorb,qorb)-c_term !sorb=rorb
                    !case 3 s and p are the different ones - sign change so both operators to come before occupied when acting on their det (left det for p and r, right det for s and q)
                    twordm(porb,rorb,qorb,rorb)=twordm(porb,rorb,qorb,rorb)-c_term !sorb=rorb




                end do

                buffer=nalphbet(1,2,c1)
                do while(buffer /= 0_8)
                    rorb=trailz(buffer)+1
                    buffer=IAND(buffer,buffer-1_8)

                    twordm(porb,rorb,rorb,qorb)=twordm(porb,rorb,rorb,qorb)+c_term

                    twordm(rorb,porb,qorb,rorb)=twordm(rorb,porb,qorb,rorb)+c_term!sorb=rorb
                    !case 2 spins and spinr are the same by construction as are spinp and spinq

                end do
                end if


              !double exc with alpha and beta
              do j=1,size(single_exc(index2,:))
                  index=0
                  civ2=0.0_dp
                  particle_int=single_exc(index2,j)

                if (particle_int==0) cycle
                  ss = starts(index1,i)
                  ee = ends(index1,i)
                index=findloc(Nalphbet(1,2,ss:ee),particle_int, 1)
                if (index==0) cycle
                  index=index+ss-1
                  civ2=civs(index)
                  rorb=hole_exc(index1,i)
                  sorb=particle_exc(index1,i)
                  porb=hole_exc(index2,j)
                  qorb=particle_exc(index2,j)

                  ep=phases(index1,i)+phases(index2,j)

                  c_term=1.00*phase_dbl(iand(ep,1))*civ1*civ2


                   twordm(porb,rorb,sorb,qorb)=twordm(porb,rorb,sorb,qorb)&
                        +c_term

                   twordm(rorb,porb,qorb,sorb)=twordm(rorb,porb,qorb,sorb)&
                       + c_term


                end do




         end do





         do i=1,size(double_exc(index1,:))
             index=0
             civ2=0.0_dp

             particle_int=double_exc(index1,i)
             if (particle_int==0) cycle
             count2=counter_exc_2(index1,i)
             if (count2==0) cycle
             ss = starts_2(index1,i)
             ee = ends_2(index1,i)
             index=findloc(Nalphbet(1,2,ss:ee),nalphbet(1,2,c1), 1)
             if (index==0) cycle
             index=index+ss-1
             civ2=civs(index)

              !samespin.eq.1 or   2
  !Case 1
            ep=phases_2(index1,i)
            sorb=particle_exc_1(index1,i) !korb
            qorb=particle_exc_2(index1,i)!lorb
            porb=hole_exc_2(index1,i) !jorb
            rorb=hole_exc_1(index1,i) !iorb

            ! if (particle_int/=nalphbet(1,2,c1) .and. (nalphbet(1,2,c1)/=nalphbet(1,1,c1))) then
                c_term=phase_dbl(iand(ep,1))*2.00*civ1*civ2

            !c_term=4.00*phase_dbl(iand(ep,1))*civ1*civ2

            twordm(porb,rorb,sorb,qorb)=twordm(porb,rorb,sorb,qorb)&
                +c_term
            !all same spin so all swaps are allowed
            !Case 2 Swap p and r introduces negative sign

            twordm(rorb,porb,sorb,qorb)=twordm(rorb,porb,sorb,qorb)&
                -c_term
   !Case 3 from Case 1 swap s and q to give negative sign
            twordm(porb,rorb,qorb,sorb)=twordm(porb,rorb,qorb,sorb)&
            -c_term
   !Case 4 from Case 1 swap s and q then swap p and r so no sign change

            twordm(rorb,porb,qorb,sorb)=twordm(rorb,porb,qorb,sorb)&
                + c_term



         end do






end do







     print*,'Value not appearing', twordm(1,3,1,4)
    print*,'everything finishes'
    cutoff = 1E-30
    count_p=0
    count2=1
    print*,sum(twordm)
    allocate(logicaltwordms(norbs**4))
    allocate(totaldum(norbs**4), matdum(norbs**4,4))
    logicaltwordms(:)=.False.
    open (unit = 15, file = 'twordm_fortran_bit.dat')
    do p=1,norbs
        do q=1,norbs
            do r=1,norbs
                do s=1,norbs
                    totaldum(count2)=twordm(p,q,r,s)
                    matdum(count2,:)=(/p,s,q,r/)
                    if (abs(twordm(p,q,r,s))>=cutoff) then
                        count_p=count_p+1
                        logicaltwordms(count2)=.True.
                        write(15,"(4I3, E30.16)") p,s,q,r,twordm(p,q,r,s)

                    end if
                    count2=count2+1
                end do
            end do
        end do
    end do
    close(15)

    allocate(mat(count_p, 4), total(count_p))
    count_p=1
    do i=1,count2-1
        if (logicaltwordms(i)) then
            mat(count_p,:)=matdum(i,:)
            total(count_p)=totaldum(i)

            count_p=count_p+1
        end if
    end do

        print*, 'twordm calculated'


!
    end subroutine createtwordm_bit_fci



     subroutine createtwordm_bit(file_read,length, mat,total)


      !length of wavefunction, energy from 2RDM (without nuclear contribution)
 ! basis functions, electrons, coefficients, list of orbitals in SD
implicit none
 INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
     character(len=30), intent(in) :: file_read

    real(kind=dp), intent(out), dimension(:), allocatable :: total
    integer(kind=ikind), intent(out), dimension(:,:), allocatable :: mat

integer(kind=ikind) ::  length
double precision, allocatable :: SpinFree2RDM(:,:,:,:)
integer(kind=ikind) porb,rorb,sorb,qorb
integer(kind=ikind) spins,spinq,spinp,spinr
double precision dtemp,ep,eg
integer(kind=ikind) ici,jci,ndiff,idiff1,idiff2,countcivs
integer(kind=ikind) i,i1,i2,l,l2,k,k2
double precision Eone_e,Etwo_e,TwoRDM_e,cutoff
double precision civ1, civ2
integer(kind=ikind) newdiff,mytemp,n,count2,count_p,p,q,r,s
  integer(kind=ikind) hole,part,tz,tz2,buffer,nword
  integer(kind=ikind) myhigh,mylow,nperm
  integer(kind=ikind) idx_part,idx_hole
  integer(kind=ikind) exc(0:2,2,2),mya,myb,myc,myd
  integer(kind=ikind) tz3
  integer(kind=ikind) myspin,samespin,ispin,nbft
  integer(kind=ikind) buffer2,ispin2, dummy,length2,buffer3
  integer(kind=ikind) buffer_esp1, buffer_esp2, allc1,allc2
  double precision, parameter :: phase_dbl(0:1)=(/1.d0,-1.d0/)
  double precision c_term,temp
  double precision :: c1(length), c2(length)
  integer(kind=ikind) icij(2,1,length)
 logical(4), dimension(:), allocatable :: logicaltwordms
    integer(kind=ikind),  dimension(:,:), allocatable :: matdum
    real(kind=dp), dimension(:), allocatable :: totaldum



open(15,file=file_read)
nword=1
print*,'reading ', length, ' lines'

 do i=1,length


            read(15,*)dummy, c1(i),c2(i),icij(1,1,i), icij(2,1,i)

     enddo
close(15)


 nbft=maxval(integer2binary_orbs(maxval(icij)))
print*,'max orbs', nbft
ALLOCATE(SpinFree2RDM(nbft,nbft,nbft,nbft))

print*,'allocated matrix'
  !set to zero (a while ago Array=0.0D0 did not work for some compilers)
 do porb=1,nbft
 do rorb=1,nbft
 do sorb=1,nbft
 do qorb=1,nbft
 SpinFree2RDM(porb,rorb,sorb,qorb)=0.0D0
 end do
 end do
 end do
 end do

print*, 'spinfreee filled'
dtemp=0.0D0
!ici=jci first


length2=length


do ici=1,length2


!only zero differences by construction  ! could try to make further improvements but there are only length terms for no differences not O(length**2) as for 1 and 2
c_term=c1(ici)*c2(ici) !calculate once as  common to all zero differences for this ici
do ispin=1,2
buffer=icij(ispin,1,ici)

do while(buffer.ne.0)
sorb=trailz(buffer)+1
buffer=IAND(buffer,buffer-1)

do ispin2=1,2
buffer2=icij(ispin2,1,ici)
do while(buffer2.ne.0)
qorb=trailz(buffer2)+1
buffer2=IAND(buffer2,buffer2-1)

if((qorb.eq.sorb).AND.(ispin.eq.ispin2)) cycle ! all possible choices of two so can't pick the same twice

  !Case 1 p is same as s and r is q
  if(ispin.eq.ispin2) THEN

SpinFree2RDM(sorb,qorb,sorb,qorb)=SpinFree2RDM(sorb,qorb,sorb,qorb)-c_term !porb=sorb rorb=qorb
  END IF

!Case 2 p is same as q and r is same as s
!All s and q spins valid for this contribution

SpinFree2RDM(qorb,sorb,sorb,qorb)=SpinFree2RDM(qorb,sorb,sorb,qorb)+c_term !p=q r=s

end do
end do
end do
end do ! end of zero differences


end do !end of ici loop

print*,'starting big loop'




!now jci>ici and double values so we don't need to do jci<ici
do ici=1,length2
do jci=1,length2

temp=SpinFree2RDM(11,11,12,14)



!!!!!!!!!!!!!
newdiff=0
do n=1,nword
mytemp=IEOR(icij(1,n,ici),icij(1,n,jci))

newdiff=newdiff+POPCNT(mytemp) ! calcs number of bits set to 1 as we used xor (IEOR) bitwise that must be - note one difference adds two to newdiff as there are two places where the bits will be set to 1

mytemp=IEOR(icij(2,n,ici),icij(2,n,jci))

newdiff=newdiff+POPCNT(mytemp)

end do
 if (newdiff.gt.4) cycle !more than two differences so matrix element is zero






if(newdiff.eq.4) THEN ! get differences and phase

exc(0,1,1)=0
exc(0,2,1) =0
exc(0,1,2)=0
exc(0,2,2)=0
nperm=0
samespin=0
do myspin=1,2
idx_part=0
idx_hole=0
mytemp=(IEOR(icij(myspin,1,ici),icij(myspin,1,jci)))
hole= IAND(mytemp,icij(myspin,1,ici))
part= IAND(mytemp,icij(myspin,1,jci))

do while(part.ne.0)
tz=trailz(part)
idx_part=idx_part+1
exc(0,2,myspin)=exc(0,2,myspin)+1
exc(idx_part,2,myspin)=tz+1
part=iand(part,part-1)
end do

do while(hole.ne.0)
tz=trailz(hole)
idx_hole=idx_hole+1
exc(0,1,myspin)=exc(0,1,myspin)+1
exc(idx_hole,1,myspin)=tz+1
hole=iand(hole,hole-1)
end do

do i=1,exc(0,1,myspin)
mylow=min(exc(i,1,myspin),exc(i,2,myspin))
myhigh=max(exc(i,1,myspin),exc(i,2,myspin))
nperm=nperm+POPCNT(IAND(icij(myspin,1,ici),&
IAND(ibset(0,myhigh-1)-1,ibclr(-1,mylow)+1)))


end do
!both holes have same spin
if(exc(0,1,myspin).eq.2) THEN !this seems less efficient...
samespin=myspin
mya=min(exc(1,1,myspin),exc(1,2,myspin))
myb=max(exc(1,1,myspin),exc(1,2,myspin))
myc=min(exc(2,1,myspin),exc(2,2,myspin))
myd=max(exc(2,1,myspin),exc(2,2,myspin))

if((myc>mya).AND.(myc<myb).AND.(myd>myb)) nperm=nperm+1
exit
END IF ! end of check if both are same spin

end do !loop over myspin

  c_term=c1(ici)*c2(jci)*phase_dbl(iand(nperm,1)) ! calculate once as common to all in this loop
  ! c_term=c(ici)*c(jci)*ep
   if(samespin.eq.0) THEN

      !Case 1
     sorb=exc(1,2,1) !korb
     qorb=exc(1,2,2)!lorb
     porb=exc(1,1,2) !jorb
     rorb=exc(1,1,1) !iorb
 ! print*,ici,jci,c_term
 ! print*, porb,rorb,sorb,qorb
      !spins.eq.spinr and spinq.eq.spinp by construction
SpinFree2RDM(porb,rorb,sorb,qorb)=SpinFree2RDM(porb,rorb,sorb,qorb)&
+c_term

SpinFree2RDM(rorb,porb,qorb,sorb)=SpinFree2RDM(rorb,porb,qorb,sorb)&
+c_term
   !Case 2 Swap p and r introduces negative sign but means spins.ne.spinr so no contribution
   !Case 3 from Case 1 swap s and q to give negative sign but spins.ne.spinr so no contribution
   !Case 4 from Case 1 swap s and q then swap p and r so no sign change

         !spins.eq.spinr and spinq.eq.spinp by construction
!SpinFree2RDM(rorb,porb,qorb,sorb)=SpinFree2RDM(rorb,porb,qorb,sorb)&
!+c_term


!SpinFree2RDM(porb,rorb,sorb,qorb)=SpinFree2RDM(porb,rorb,sorb,qorb)&
!+c_term
    ELSE
  !samespin.eq.1 or   2
  !Case 1
     sorb=exc(1,2,samespin) !korb
     qorb=exc(2,2,samespin)!lorb
     porb=exc(2,1,samespin) !jorb
     rorb=exc(1,1,samespin) !iorb
SpinFree2RDM(porb,rorb,sorb,qorb)=SpinFree2RDM(porb,rorb,sorb,qorb)&
+c_term

!SpinFree2RDM(qorb,sorb,rorb,porb)=SpinFree2RDM(qorb,sorb,rorb,porb)&
!+c_term
!all same spin so all swaps are allowed
   !Case 2 Swap p and r introduces negative sign

SpinFree2RDM(rorb,porb,sorb,qorb)=SpinFree2RDM(rorb,porb,sorb,qorb)&
-c_term

!SpinFree2RDM(sorb,qorb,rorb,porb)=SpinFree2RDM(sorb,qorb,rorb,porb)&
!-c_term
   !Case 3 from Case 1 swap s and q to give negative sign
SpinFree2RDM(porb,rorb,qorb,sorb)=SpinFree2RDM(porb,rorb,qorb,sorb)&
-c_term

!SpinFree2RDM(qorb,sorb,porb,rorb)=SpinFree2RDM(qorb,sorb,porb,rorb)&
!-c_term
   !Case 4 from Case 1 swap s and q then swap p and r so no sign change

SpinFree2RDM(rorb,porb,qorb,sorb)=SpinFree2RDM(rorb,porb,qorb,sorb)&
+c_term

   !    SpinFree2RDM(sorb,qorb,porb,rorb)=SpinFree2RDM(sorb,qorb,porb,rorb)&
!+c_term


   end if ! end of samespin check



cycle


end if ! end of two differences (newdiff.eq.4)

!one difference
if(newdiff.eq.2) THEN !  newdiff.eq.2 is one difference

do myspin=1,2
mytemp=(IEOR(icij(myspin,1,ici),icij(myspin,1,jci)))
hole= IAND(mytemp,icij(myspin,1,ici))
part= IAND(mytemp,icij(myspin,1,jci))
if(hole.ne.0) THEN
  tz=1+trailz(hole)
  tz2=1+trailz(part)
  EXIT
  END IF
  end do


!get sign for single
!myspin shows if alpha or beta at this point
mylow=min(tz,tz2)
myhigh=max(tz,tz2)
!PRINT *,mylow,myhigh
nperm=POPCNT(IAND(icij(myspin,1,ici),IAND(ibset(0,myhigh-1)-1,ibclr(-1,mylow)+1)))



  c_term=c1(ici)*c2(jci)*phase_dbl(iand(nperm,1)) ! calculate once as common to all in this loop
    !c_term=c(ici)*c(jci)*ep
!case 1 and 4 and 2 and 3
porb=tz
qorb=tz2

!ispin=myspin
buffer=IBCLR(icij(myspin,1,ici),porb-1) ! remove porb of this spin so rorb cannot be porb


do while(buffer.ne.0)
rorb=trailz(buffer)+1
buffer=IAND(buffer,buffer-1)


!case 1
SpinFree2RDM(porb,rorb,rorb,qorb)=SpinFree2RDM(porb,rorb,rorb,qorb)+c_term  !sorb=rorb
!case 4 s and r are the differences so signs of moving through occupied will cancel
SpinFree2RDM(rorb,porb,qorb,rorb)=SpinFree2RDM(rorb,porb,qorb,rorb)+c_term !sorb=rorb
!case 2 spins and spinr are the same by construction as are spinp and spinq
SpinFree2RDM(rorb,porb,rorb,qorb)=SpinFree2RDM(rorb,porb,rorb,qorb)-c_term !sorb=rorb
! case 3 s and p are the different ones - sign change so both operators to come before occupied when acting on their det (left det for p and r, right det for s and q)
SpinFree2RDM(porb,rorb,qorb,rorb)=SpinFree2RDM(porb,rorb,qorb,rorb)-c_term !sorb=rorb


!SpinFree2RDM(qorb,rorb,rorb,porb)=SpinFree2RDM(qorb,rorb,rorb,porb)+c_term  !sorb=rorb
!!case 4 s and r are the differences so signs of moving through occupied will cancel
!SpinFree2RDM(rorb,qorb,porb,rorb)=SpinFree2RDM(rorb,qorb,porb,rorb)+c_term !sorb=rorb
!!case 2 spins and spinr are the same by construction as are spinp and spinq
!SpinFree2RDM(rorb,qorb,rorb,porb)=SpinFree2RDM(rorb,qorb,rorb,porb)-c_term !sorb=rorb
!! case 3 s and p are the different ones - sign change so both operators to come before occupied when acting on their det (left det for p and r, right det for s and q)
!SpinFree2RDM(qorb,rorb,porb,rorb)=SpinFree2RDM(qorb,rorb,porb,rorb)-c_term !sorb=rorb


!  if (rorb==1 .and. porb==4 .and. qorb==3) then
!      print*, SpinFree2RDM(rorb,porb,qorb,rorb),ici,jci
!  end if

end do

buffer=icij(IEOR(myspin,3),1,ici)  !map spin 1 to spin 2 and 2 to 1

do while(buffer.ne.0)
rorb=trailz(buffer)+1
buffer=IAND(buffer,buffer-1)
!only have case 1 and 2 as spins are different
!case 1
SpinFree2RDM(porb,rorb,rorb,qorb)=SpinFree2RDM(porb,rorb,rorb,qorb)+c_term  !sorb=rorb
!case 4 s and r are the differences so signs of moving through occupied will cancel
SpinFree2RDM(rorb,porb,qorb,rorb)=SpinFree2RDM(rorb,porb,qorb,rorb)+c_term !sorb=rorb
!
!SpinFree2RDM(qorb,rorb,rorb,porb)=SpinFree2RDM(qorb,rorb,rorb,porb)+c_term  !sorb=rorb
!!case 4 s and r are the differences so signs of moving through occupied will cancel
!SpinFree2RDM(rorb,qorb,porb,rorb)=SpinFree2RDM(rorb,qorb,porb,rorb)+c_term !sorb=rorb

end do

!end of case 1 and case 4 and 2 and 3



!if (SpinFree2RDM(4,5,12,12) /= temp) print*,ici,jci,ep, SpinFree2RDM(4,5,12,12)

cycle
end if  ! end of  single difference, newdiff.eq.2





end do !end of ici loop
end do !end of jci loop

print*,'here we are'

print*, maxval(SpinFree2RDM)

print*,'final guy', SpinFree2RDM(11,11,12,14)

!Write out 2RDM


 print*,'everything finishes'
    cutoff = 1E-30
    count_p=0
    count2=1

    allocate(logicaltwordms(nbft**4))
    allocate(totaldum(nbft**4), matdum(nbft**4,4))
    logicaltwordms(:)=.False.
    open (unit = 15, file = 'twordm_fortran_bit_2.dat')
    do p=1,nbft
        do q=1,nbft
            do r=1,nbft
                do s=1,nbft
                    totaldum(count2)=SpinFree2RDM(p,q,r,s)
                    matdum(count2,:)=(/p,s,q,r/)
                    if (abs(SpinFree2RDM(p,q,r,s))>=cutoff) then
                        count_p=count_p+1
                        logicaltwordms(count2)=.True.
                        write(15,"(4I3, E30.16)") p,s,q,r,SpinFree2RDM(p,q,r,s)

                    end if
                    count2=count2+1
                end do
            end do
        end do
    end do
    close(15)

    allocate(mat(count_p, 4), total(count_p))
    count_p=1
    do i=1,count2-1
        if (logicaltwordms(i)) then
            mat(count_p,:)=matdum(i,:)
            total(count_p)=totaldum(i)

            count_p=count_p+1
        end if
    end do

        print*, 'twordm calculated'


end subroutine

subroutine one_rdm_two_rdm(mat,twordm,onerdm,nel)
    implicit none

    integer(kind=SELECTED_INT_KIND(8)),intent(in), dimension(:,:):: mat
    double precision, intent(in), dimension(:) :: twordm
    double precision, intent(out), dimension(:,:), allocatable::onerdm
    integer(kind=ikind), intent(in) :: nel
    integer(kind=ikind) :: maxnmo,i,j,k,l

    maxnmo=maxval(mat)
    print*,'maximum number or orbitals', maxnmo
    allocate(onerdm(maxnmo,maxnmo))
    onerdm=0.d0
    do i=1,size(mat(:,1))
        j=mat(i,1)
        k=mat(i,2)
        if (mat(i,3)==mat(i,4)) then
        onerdm(j,k)=onerdm(j,k)+twordm(i)/(nel-1)

        end if
    end do







    end subroutine




       subroutine one_rdm_bit(file_read,onerdm, maxnmo,numberlines,newdat,irep)

           Use, intrinsic :: iso_fortran_env, Only : iostat_end
           implicit none

           INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
           INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)


           integer(kind=ikind):: numberlines
           character(len=30),intent(in) :: file_read
           integer(kind=ikind),intent(in),dimension(:):: irep
           real(kind=dp), dimension(:), allocatable :: civs
           integer(kind=ikind),dimension(:,:,:), allocatable :: Nalphbet
           integer(kind=ikind):: j, i, error_1,count,temp,sp,diff1,diff2,popcnt1,popcnt2,compnum,buffer_prime, buffer_prime_2
            integer(kind=ikind):: particle_int,hole_int,k
           real(kind=dp),dimension(:,:), allocatable :: onerdm
           integer(kind=ikind),intent(inout):: maxnmo
           integer(kind=ikind):: exc(0:2,2,2),deg,buffer,c1,c2,ep,n1a,n2a,n1b,n2b,low,high,number,red_vec(maxnmo*2)
           double precision :: phase, c,civ1,civ2
           integer(kind=ikind):: tmp,tz,final,max,prueba(maxnmo*2),buffer_2,tz2,index,nn,nel,n,porb,qorb,newdiff,diff
           integer(kind=ikind):: mylow, myhigh
           real(kind=dp)::time1,time2
           integer(kind=ikind), intent(in), dimension(:,:):: newdat
           integer(kind=ikind), allocatable, dimension(:,:,:) :: orb_mat
           integer(kind=ikind),allocatable, dimension(:) :: orb_indexed
           integer(kind=ikind), dimension(numberlines):: sum_mat
              double precision, parameter :: phase_dbl(0:1)=(/1.d0,-1.d0/)


           onerdm=0.0_dp
           count=numberlines


           allocate(civs(count), Nalphbet(1,2,count))
           civs=0
           Nalphbet=0
           open(15,file=file_read)
           do i=1,count

               read(15,*)j, civs(i), Nalphbet(1,1,i), Nalphbet(1,2,i)
           enddo
           close(15)

           maxnmo=maxval(integer2binary_orbs(maxval(Nalphbet)))

           print*,'max orbitals',maxnmo
           nel=popcnt(Nalphbet(1,1,1))
            allocate(onerdm(maxnmo,maxnmo))

           onerdm=0.0_dp

          do c1=1,numberlines

              do i=1,2

                  buffer=Nalphbet(1,i,c1)
                  do while(buffer/=0_8)
                      j=trailz(buffer)

                      onerdm(j+1,j+1)=onerdm(j+1,j+1)+civs(c1)*civs(c1)
                      buffer=iand(buffer,buffer-1_8)
                  end do

              end do

           end do
          print*,'we have finished the elastic part'
            call cpu_time(time1)
           count=numberlines
           do c1=1,count-1

               if (c1==100000) then


                  call cpu_time(time2)
                   print*,'almost there', time2-time1
               end if
               do c2=c1+1,count

                    diff1=popcnt(IEOR(Nalphbet(1,1,c1),Nalphbet(1,1,c2)))
                    diff2=popcnt(IEOR(Nalphbet(1,2,c1),Nalphbet(1,2,c2)))
                    newdiff=diff1+diff2
                    if (newdiff/=2) cycle

                    if (diff1==2) then
                            diff=IEOR(Nalphbet(1,1,c1),Nalphbet(1,1,c2))
                            porb=trailz(iand(Nalphbet(1,1,c1),diff))+1
                            qorb=trailz(iand(Nalphbet(1,1,c2),diff))+1

                            mylow=min(porb,qorb)
                            myhigh=max(porb,qorb)
                            ep=POPCNT(IAND(Nalphbet(1,1,c1),&
                                IAND(ibset(0_8,myhigh-1_8)-1_8,ibclr(-1_8,mylow)+1_8)))
!                            ep=ep+POPCNT(IAND(Nalphbet(1,2,c1),&
!                               IAND(ibset(0,myhigh-1)-1,ibclr(-1,mylow))))



                            onerdm(porb,qorb)=onerdm(porb,qorb)+civs(c1)*civs(c2)*phase_dbl(iand(ep,1_8))

                            onerdm(qorb,porb)=onerdm(qorb,porb)+civs(c1)*civs(c2)*phase_dbl(iand(ep,1_8))

                    end if

                   if (diff2==2) then
                            diff=IEOR(Nalphbet(1,2,c1),Nalphbet(1,2,c2))
                            porb=trailz(iand(Nalphbet(1,2,c1),diff))+1
                            qorb=trailz(iand(Nalphbet(1,2,c2),diff))+1

                            mylow=min(porb,qorb)
                            myhigh=max(porb,qorb)
                            ep=POPCNT(IAND(Nalphbet(1,2,c1),&
                                IAND(ibset(0_8,myhigh-1_8)-1_8,ibclr(-1,mylow)+1_8)))
                             !ep=ep+POPCNT(IAND(Nalphbet(1,1,c1),&
                             !   IAND(ibset(0,myhigh-1)-1,ibclr(-1,mylow)+1)))


                          !  if (c1<=20) print*,'conguration of interest',c1,c2,ep

                            onerdm(porb,qorb)=onerdm(porb,qorb)+civs(c1)*civs(c2)*phase_dbl(iand(ep,1_8))

                            onerdm(qorb,porb)=onerdm(qorb,porb)+civs(c1)*civs(c2)*phase_dbl(iand(ep,1_8))

                    end if


                   end do! calcs number of bits set to 1 as we used xor (IEOR) bitwise that must be - note one difference adds two to newdiff as there are two places where the bits will be set to 1
enddo


             call cpu_time(time2)
           print*,'time constructing 1rdm', time2-time1
           open(file='onerdm_mine.dat', unit=15)
             do i=1,size(onerdm(:,1))

                 write(15,'(1000F14.7)')( onerdm(i,c1) ,c1=1,size(onerdm(i,:)))
            end do
!           call compute_density_matrix(Nalphbet,count-1,civs,maxnmo, &
!               1,onerdm)
!
!
!               open(file='onerdm_2.dat', unit=15)
!             do i=1,size(onerdm(:,1))
!
!                 write(15,'(1000F14.7)')( onerdm(i,c1) ,c1=1,size(onerdm(i,:)))
!            end do
!    close(15)

       end subroutine



       subroutine compute_density_matrix(det,Ndet,coef,mo_num, &
               Nint,density_matrix)
           implicit none
           INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
           integer(kind=ikind), intent(in) :: det(Nint,2,Ndet)
           integer, intent(in) :: Ndet, Nint, mo_num
           double precision, intent(in) :: coef(Ndet)
           double precision, intent(out) :: density_matrix(mo_num,mo_num)
           integer(kind=ikind) :: i,j,k,l,ispin,ishift
           integer(kind=ikind):: buffer
           integer(kind=ikind) :: deg
           integer(kind=ikind) :: exc(0:2,2,2)
           double precision :: phase, c
           integer(kind=ikind) :: n_excitations
           density_matrix = 0.d0



           do k=1,Ndet
               do ispin=1,2
                   ishift = 1
                   do i=1,Nint
                       buffer = det(i,ispin,k)
                       do while (buffer /= 0_8)
                           j = trailz(buffer) + ishift
                           density_matrix(j,j) = density_matrix(j,j) &
                                   + coef(k)*coef(k)
                           buffer = iand(buffer,buffer-1_8)
                       end do
                       ishift = ishift+64
                   end do
               end do

               do l=1,k-1
                   if (excitations(det(1,1,k),det(1,1,l),Nint) /= 1) then
                       cycle
                   end if
                   call get_excitation(det(1,1,k),det(1,1,l),exc,deg,phase,Nint)
                   if (exc(0,1,1) == 1) then
                       i = exc(1,1,1)
                       j = exc(1,2,1)


                   else
                       i = exc(1,1,2)
                       j = exc(1,2,2)

                   end if
                   c = phase*coef(k)*coef(l)

                   c = c+c
                   density_matrix(j,i) = density_matrix(j,i) + c
                   density_matrix(i,j) = density_matrix(i,j) + c
               end do
           end do
         !   print*,'diag1',density_matrix(26,26)
       end subroutine

       subroutine get_double_excitation(det1,det2,exc,phase,Nint)
           implicit none
           integer(kind=ikind), intent(in) :: Nint
           integer(kind=ikind), intent(in) :: det1(Nint,2), det2(Nint,2)
           integer(kind=ikind), intent(out) :: exc(0:2,2,2)
           double precision, intent(out) :: phase
           integer(kind=ikind) :: l, ispin, idx_hole, idx_particle, ishift
           integer(kind=ikind) :: i,j,k,m,n,high, low,a,b,c,d,nperm,tz,nexc
           integer(kind=ikind) :: hole, particle, tmp
           double precision, parameter :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)
           exc(0,1,1) = 0
           exc(0,2,1) = 0
           exc(0,1,2) = 0
           exc(0,2,2) = 0
           nexc=0
           nperm=0
           do ispin = 1,2
               idx_particle = 0
               idx_hole = 0
               ishift = -63
               do l=1,Nint
                   ishift = ishift + 64
                   if (det1(l,ispin) == det2(l,ispin)) then
                       cycle
                   end if
                   tmp = xor( det1(l,ispin), det2(l,ispin) )
                   particle = iand(tmp, det2(l,ispin))
                   hole = iand(tmp, det1(l,ispin))
                   do while (particle /= 0_8)
                       tz = trailz(particle)
                       nexc = nexc+1
                       idx_particle = idx_particle + 1
                       exc(0,2,ispin) = exc(0,2,ispin) + 1
                       exc(idx_particle,2,ispin) = tz+ishift
                       particle = iand(particle,particle-1_8)
                   end do
                   do while (hole /= 0_8)
                       tz = trailz(hole)
                       nexc = nexc+1
                       idx_hole = idx_hole + 1
                       exc(0,1,ispin) = exc(0,1,ispin) + 1
                       exc(idx_hole,1,ispin) = tz+ishift
                       hole = iand(hole,hole-1_8)
                   end do
                   if (nexc == 4) exit
               end do

               do i=1,exc(0,1,ispin)
                   low = min(exc(i,1,ispin),exc(i,2,ispin))
                   high = max(exc(i,1,ispin),exc(i,2,ispin))
                   j = ishft(low-1,-6)+1
                   n = iand(low,63)
                   k = ishft(high-1,-6)+1
                   m = iand(high,63)
                   if (j==k) then
                       nperm = nperm + popcnt(iand(det1(j,ispin), &
                               iand( ibset(0_8,m-1)-1_8, ibclr(-1_8,n)+1_8 ) ))
                   else
                       nperm = nperm + popcnt(iand(det1(k,ispin), &
                               ibset(0_8,m-1)-1_8)) &
                               + popcnt(iand(det1(j,ispin), &
                                       ibclr(-1_8,n) +1_8))
                       do l=j+1,k-1
                           nperm = nperm + popcnt(det1(l,ispin))
                       end do
                   end if
               end do
               if (exc(0,1,ispin) == 2) then
                   a = min(exc(1,1,ispin), exc(1,2,ispin))
                   b = max(exc(1,1,ispin), exc(1,2,ispin))
                   c = min(exc(2,1,ispin), exc(2,2,ispin))
                   d = max(exc(2,1,ispin), exc(2,2,ispin))
                   if (c>a .and. c<b .and. d>b) nperm = nperm + 1
                   exit
               end if
           end do
           phase = phase_dble(iand(nperm,1))

       end

       subroutine get_single_excitation(det1,det2,exc,phase,Nint)
           implicit none
           integer(kind=ikind), intent(in) :: Nint
           integer(kind=ikind), intent(in) :: det1(Nint,2)
           integer(kind=ikind), intent(in) :: det2(Nint,2)
           integer(kind=ikind), intent(out) :: exc(0:2,2,2)
           double precision, intent(out) :: phase
           integer(kind=ikind) :: tz, l, ispin, ishift, nperm, i, j, k, m, n, high, low
           integer(kind=ikind) :: hole, particle, tmp
           double precision, parameter :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)

           exc(0,1,1) = 0
           exc(0,2,1) = 0
           exc(0,1,2) = 0
           exc(0,2,2) = 0
           do ispin = 1,2
               ishift = -63
               do l=1,Nint
                   ishift = ishift + 64
                   if (det1(l,ispin) == det2(l,ispin)) cycle
                   tmp = xor( det1(l,ispin), det2(l,ispin) )
                   particle = iand(tmp, det2(l,ispin))
                   hole = iand(tmp, det1(l,ispin))
                   if (particle /= 0_8) then
                       tz = trailz(particle)
                       exc(0,2,ispin) = 1
                       exc(1,2,ispin) = tz+ishift
                   end if
                   if (hole /= 0_8) then
                       tz = trailz(hole)
                       exc(0,1,ispin) = 1
                       exc(1,1,ispin) = tz+ishift
                   end if
                    phase=1
!                   if ( iand(exc(0,1,ispin),exc(0,2,ispin)) == 1 ) then
!                       low = min(exc(1,1,ispin),exc(1,2,ispin))
!
!                       high = max(exc(1,1,ispin),exc(1,2,ispin))
!                       j = ishft(low-1,-6)+1
!                       n = iand(low,63)
!                       k = ishft(high-1,-6)+1
!                       m = iand(high,63)
!                       if (j==k) then
!                           nperm = popcnt(iand(det1(j,ispin), &
!                                   iand( ibset(0_8,m-1)-1_8, ibclr(-1_8,n)+1_8 ) ))
!                       else
!                           nperm = popcnt(iand(det1(k,ispin), ibset(0_8,m-1)-1_8)) + &
!                                   popcnt(iand(det1(j,ispin), ibclr(-1_8,n) +1_8))
!                           do i=j+1,k-1
!                               nperm = nperm + popcnt(det1(i,ispin))
!                           end do
!                       end if
!                       phase = phase_dble(iand(nperm,1))
!                       return
!                   end if
               end do
           end do
       end


       subroutine get_excitation(det1,det2,exc,degree,phase,Nint)
           implicit none
           integer(kind=ikind), intent(in) :: Nint
           integer(kind=ikind), intent(in) :: det1(Nint,2), det2(Nint,2)
           integer(kind=ikind), intent(out) :: exc(0:2,2,2)
           integer(kind=ikind), intent(out) :: degree
           double precision, intent(out) :: phase



           degree = excitations(det1,det2,Nint)

           select case (degree)

           case (3:)
               degree = -1
               return

           case (2)
               call get_double_excitation(det1,det2,exc,phase,Nint)
               return

           case (1)
               call get_single_excitation(det1,det2,exc,phase,Nint)
               return

           case(0)
               return

           end select
       end

       subroutine mcci_to_bit(file_read,file_write,numberlines)
           implicit none
           character(len=60), intent(in):: file_read,file_write
           integer(kind=ikind),intent(in):: numberlines
           integer(kind=ikind):: icij(2,1,numberlines), j, k,i,nbft,itemp,ici,ntotal,phase,iperm
           integer(kind=ikind), allocatable, dimension(:,:):: list
           real(kind=dp):: c1(numberlines), c2(numberlines)
           double precision, parameter :: phase_dbl(0:1) = (/ 1.d0, -1.d0 /)


           open(file=file_read,unit=15)

           do i=1,numberlines
               read(15,*)ici, c1(i), c2(i), icij(1,1,i), icij(2,1,i)
           enddo
           close(15)


           nbft=maxval(integer2binary_orbs(maxval(icij)))

           allocate(list(2,popcnt(icij(1,1,1))*2))
           open(file=file_write,unit=17)
           do ici=1,numberlines
                k=0
               do j=0,nbft-1
                   if(btest(icij(1,1,ici), j)) THEN
                       k=k+1
                       list(1,k) = j+1
                   END IF
                   if(btest(icij(2,1,ici), j)) THEN
                       k=k+1
                       list(1,k) = j+1 + nbft
                   END IF

               end do

               ntotal=2*popcnt(icij(1,1,1))
               !alpha then beta format  'bitwise rep'
               k=0
               do j=0,nbft-1
                   if(btest(icij(1,1,ici), j)) THEN
                       k=k+1
                       list(2,k) = j+1
                   END IF
               end do
               do j=0,nbft-1
                   if(btest(icij(2,1,ici), j)) THEN
                       k=k+1
                       list(2,k) = j+1 + nbft
                   END IF
               end do


               iperm=0  !moves to transform list(2  to list(1

               do j=1,ntotal-1
                   if(list(1,j).ne.list(2,j)) THEN
                       do k=j+1,ntotal
                           if(list(1,j).eq.list(2,k)) THEN
                               iperm=iperm+1
                               itemp=list(2,j)
                               list(2,j)=list(2,k)
                               list(2,k)=itemp
                               exit
                           end if
                       end do
                   END IF
               end do


               phase=phase_dbl(iand(iperm,1_8))

               c1(ici)=phase*c1(ici)
               c2(ici)=phase*c2(ici)

               write(17,"(I15, 2E20.10, 2I20)")ici,c1(ici),c2(ici), icij(1,1,ici),icij(2,1,ici)


           end do !end of loop over ici

           close(17)


           end subroutine mcci_to_bit
       recursive function combo(n,k) result(cmb)
           implicit none
           integer(kind=ikind) :: cmb
           integer(kind=ikind), intent(in) :: n,k
           integer(kind=ikind) :: mm(100,100)
    if (k == n) then
    cmb = real(1,16)
 else if (k == 1) then
    cmb = real(n,16)
 else if (mm(n,k) /=0)  then
    cmb = mm(n,k)
 else if ((k /= 1) .and. (k /= n)) then
    cmb = combo(n-1,k-1) + combo(n-1,k)
    mm(n,k) = cmb
 end if
 end function

   end module
