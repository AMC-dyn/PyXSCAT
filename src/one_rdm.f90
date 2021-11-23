
   module onerdm
       use twordms
       contains

        function excitations(det1,det2,Nint) result(n_excitations)
           implicit none
           integer*8, intent(in) :: det1(Nint,2), det2(Nint,2)
           integer , intent(in) :: Nint
           integer :: n_excitations

           integer :: l

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

       subroutine createtwordm_bit(file_read,numberlines,newdat,irep,start,end,mat,total)


     Use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
     character(len=20), intent(in) :: file_read
     integer(kind=ikind) :: norbs
     integer(KIND=IKIND), intent(in) :: numberlines
     integer*8, intent(in), dimension(:,:):: newdat
     integer*8,intent(in),dimension(:):: irep,start,end
     double precision, parameter :: phase_dbl(0:1)=(/1.d0,-1.d0/)
    real(kind=dp), intent(out), dimension(:), allocatable :: total
    integer(kind=ikind), intent(out), dimension(:,:), allocatable :: mat

    real(kind=dp), dimension(:,:,:,:), allocatable :: twordm
    real(kind=dp), dimension(:,:), allocatable :: onerdm
    integer*8 :: ep, i, j, k, l, compnum, c1
     integer*8:: X_vec(numberlines)
     integer*8, dimension(:,:,:), allocatable :: nalphbet
     real(kind=dp), dimension(:), allocatable :: civs
    real(kind=dp) :: time1, time2, cutoff
     double precision:: civ1,civ2,c,c_term,phase
    integer*8 :: buffer, buffer_2, hole_int, particle_int, buffer_p, index, eg,index1
    integer*8 :: porb,rorb,qorb,sorb, count2, p, q, r, s, dummy, n_single

    logical(4), dimension(:), allocatable :: logicaltwordms
    integer(kind=ikind),  dimension(:,:), allocatable :: matdum
    real(kind=dp), dimension(:), allocatable :: totaldum
    integer*8, dimension(:,:), allocatable :: single_exc,hole_exc, particle_exc,phases,counter_exc,starts,ends
    integer*8, dimension(:,:), allocatable :: double_exc,hole_exc_1,hole_exc_2, particle_exc_1,particle_exc_2, phases_2, counter_exc_2,starts_2,ends_2
    integer*8, dimension(:,:,:), allocatable :: beta_excs
    integer*8 :: buffer_prime, buffer_prime_2, tz,tz2, count_p,n,ss,ee,spin11,spin22,buffer2,buffer_3,buffer_4,newdat_1,len_rm

    integer*8 :: n_double,count,nel,index2,myhigh,mylow,myhigh1,myhigh2,mylow1,mylow2,count_total
    integer*8,dimension(:), allocatable :: find_ind,find_ind_2,all_indexes,occ,unocc
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
     print*,'position of jci',

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
                      ep=POPCNT(IAND(newdat(i,1),IAND(ibset(0,myhigh-1)-1,ibclr(-1,mylow)+1)))
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
     print*,'created double excitations', time2-time1, count2, n_double
























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

!               index=0
!                civ2=0.0_dp
!                particle_int=single_exc(index1,i)
!                if (particle_int==0) cycle
!                count2=counter_exc(index1,i)
!                if (count2==0) cycle
!                ss = starts(index1,i)
!                ee = ends(index1,i)
!
!
!                index=findloc(Nalphbet(1,2,ss:ee),nalphbet(1,2,c1), 1)
!
!                if (index==0) cycle
!
!                index=index+ss-1
!                civ2=civs(index)
!                porb=hole_exc(index1,i)
!                qorb=particle_exc(index1,i)
!                ep=phases(index1,i)
!
!                buffer=IBCLR(Nalphbet(1,1,c1),porb-1)
!                !if (particle_int/=nalphbet(1,2,c1) ) then
!                c_term=phase_dbl(iand(ep,1))*2.00*civ1*civ2
!              !  else
!                !     c_term=phase_dbl(iand(ep,1))*2.00*civ1*civ2
!              !  end if
!
!                do while(buffer /= 0_8)
!                    rorb=trailz(buffer)+1
!                    buffer=IAND(buffer,buffer-1_8)
!
!                    twordm(porb,rorb,rorb,qorb)=twordm(porb,rorb,rorb,qorb)+c_term
!
!                    twordm(rorb,porb,qorb,rorb)=twordm(rorb,porb,qorb,rorb)+c_term !sorb=rorb
!                    !case 2 spins and spinr are the same by construction as are spinp and spinq
!                    twordm(rorb,porb,rorb,qorb)=twordm(rorb,porb,rorb,qorb)-c_term !sorb=rorb
!                    !case 3 s and p are the different ones - sign change so both operators to come before occupied when acting on their det (left det for p and r, right det for s and q)
!                    twordm(porb,rorb,qorb,rorb)=twordm(porb,rorb,qorb,rorb)-c_term !sorb=rorb
!
!
!
!
!                end do
!
!                buffer=nalphbet(1,2,c1)
!                do while(buffer /= 0_8)
!                    rorb=trailz(buffer)+1
!                    buffer=IAND(buffer,buffer-1_8)
!
!                    twordm(porb,rorb,rorb,qorb)=twordm(porb,rorb,rorb,qorb)+c_term
!
!                    twordm(rorb,porb,qorb,rorb)=twordm(rorb,porb,qorb,rorb)+c_term!sorb=rorb
!                    !case 2 spins and spinr are the same by construction as are spinp and spinq
!
!                end do




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






!
!!!!
!!
!!!!                    !case 1
                        !print*,porb,rorb,qorb
                                        !    !sorb=rorb
!!!!                    !case 4 s and r are the differences so signs of moving through occupied will cancel
!!!!                    twordm(rorb,porb,qorb,rorb)=twordm(rorb,porb,qorb,rorb)+c_term !sorb=rorb
!!!!                    !case 2 spins and spinr are the same by construction as are spinp and spinq
!!!!                    twordm(rorb,porb,rorb,qorb)=twordm(rorb,porb,rorb,qorb)-c_term !sorb=rorb
!!!!!                    case 3 s and p are the different ones - sign change so both operators to come before occupied when acting on their det (left det for p and r, right det for s and q)
!!!!                    twordm(porb,rorb,qorb,rorb)=twordm(porb,rorb,qorb,rorb)-c_term !sorb=rorb
!!!!
!                                        end do
!
!                                    exit
!                                endif
!                            end do

!                            if (index==0) print*,'wtf'


!!
!                buffer_p=nalphbet(2,c1) !map spin 1 to spin 2 and 2 to 1
!
!                do while(buffer_p/=0_8)
!                    rorb=trailz(buffer_p)+1
!                    buffer_p=IAND(buffer_p,buffer_p-1_8)
!!                   only have case 1 and 2 as spins are different
!                    !case 1
!                    twordm(porb,rorb,rorb,qorb)=twordm(porb,rorb,rorb,qorb)+c_term  !sorb=rorb
!                    !case 4 s and r are the differences so signs of moving through occupied will cancel
!                    twordm(rorb,porb,qorb,rorb)=twordm(rorb,porb,qorb,rorb)+c_term !sorb=rorb
!                end do










!        do c2=c1+1,numberlines
!            civ2=civs(c2)
!!            call combine_alpha_beta(Nalphbet(1,c2), Nalphbet(2,c2), alphbet_2)
!
!
!            diff1=xor( Nalphbet(1,c1), Nalphbet(1,c2))
!            diff2=xor(Nalphbet(2,c1), Nalphbet(2,c2))
!            if (popcnt(diff1)==4 .and.popcnt(diff2)==0 ) then
!
!                ep=1
!                hole= iand(diff1,Nalphbet(1,c1))
!                particle = iand(diff1,Nalphbet(1,c2))
!                if (particle/= 0_8) then
!                   i=trailz(particle)+1
!                   particle=iand(particle,particle-1_8)
!                   j=trailz(particle)+1
!                   high1=i
!                   high2=j
!               end if
!               if (hole/= 0_8) then
!                   k=trailz(hole)+1
!                   hole=iand(hole,hole-1_8)
!                   l=trailz(hole)+1
!                   low1=k
!                   low2=l
!               end if
!
!
!
!                 buffer_prime=Nalphbet(1,c1)
!                 buffer_prime_2=Nalphbet(1,c2)
!
!                  do while (buffer_prime/=0_8)
!                      tz=trailz(buffer_prime)
!                      tz2=trailz(buffer_prime_2)
!
!                      if (tz/=tz2 .and. btest(Nalphbet(1,c2),tz)) then
!
!                          ep=-ep
!                      end if
!                      buffer_prime=iand(buffer_prime, buffer_prime-1_8)
!                      buffer_prime_2=iand(buffer_prime_2, buffer_prime_2-1_8)
!                  end do
!                c_term=2*civ1*civ2*ep
!                sorb=i !korb
!                qorb=j!lorb
!                porb=l !jorb
!                rorb=k !iorb
!                twordm(porb,rorb,sorb,qorb)=twordm(porb,rorb,sorb,qorb)&
!                    +c_term
!!all same spin so all swaps are allowed
!   !Case 2 Swap p and r introduces negative sign
!
!                twordm(rorb,porb,sorb,qorb)=twordm(rorb,porb,sorb,qorb)&
!                    -c_term
!   !Case 3 from Case 1 swap s and q to give negative sign
!                twordm(porb,rorb,qorb,sorb)=twordm(porb,rorb,qorb,sorb)&
!                    -c_term
!   !Case 4 from Case 1 swap s and q then swap p and r so no sign change
!
!                twordm(rorb,porb,qorb,sorb)=twordm(rorb,porb,qorb,sorb)&
!                    +c_term
!
!
!            else if (popcnt(diff1)==0 .and.popcnt(diff2)==4 ) then
!
!                ep=1
!                hole= iand(diff2,Nalphbet(2,c1))
!                particle = iand(diff2,Nalphbet(2,c2))
!                if (particle/= 0_8) then
!                   i=trailz(particle)+1
!                   particle=iand(particle,particle-1_8)
!                   j=trailz(particle)+1
!                   high1=i
!                   high2=j
!               end if
!               if (hole/= 0_8) then
!                   k=trailz(hole)+1
!                   hole=iand(hole,hole-1_8)
!                   l=trailz(hole)+1
!                   low1=k
!                   low2=l
!               end if
!
!
!
!                 buffer_prime=Nalphbet(2,c1)
!                 buffer_prime_2=Nalphbet(2,c2)
!
!                  do while (buffer_prime/=0_8)
!                      tz=trailz(buffer_prime)
!                      tz2=trailz(buffer_prime_2)
!
!                      if (tz/=tz2 .and. btest(nalphbet(2,c2),tz)) then
!
!                          ep=-ep
!                      end if
!                      buffer_prime=iand(buffer_prime, buffer_prime-1_8)
!                      buffer_prime_2=iand(buffer_prime_2, buffer_prime_2-1_8)
!                  end do
!                c_term=2*civ1*civ2*ep
!                sorb=i !korb
!                qorb=j!lorb
!                porb=l !jorb
!                rorb=k !iorb
!                twordm(porb,rorb,sorb,qorb)=twordm(porb,rorb,sorb,qorb)&
!                    +c_term
!!all same spin so all swaps are allowed
!   !Case 2 Swap p and r introduces negative sign
!
!                twordm(rorb,porb,sorb,qorb)=twordm(rorb,porb,sorb,qorb)&
!                    -c_term
!   !Case 3 from Case 1 swap s and q to give negative sign
!                twordm(porb,rorb,qorb,sorb)=twordm(porb,rorb,qorb,sorb)&
!                    -c_term
!   !Case 4 from Case 1 swap s and q then swap p and r so no sign change
!
!                twordm(rorb,porb,qorb,sorb)=twordm(rorb,porb,qorb,sorb)&
!                    +c_term
!
!   else if (popcnt(diff1)==2 .and.popcnt(diff2)==2 ) then
!
!                ep=1
!                hole= iand(diff1,Nalphbet(1,c1))
!                particle = iand(diff1,Nalphbet(1,c2))
!                if (particle/= 0_8) then
!                   i=trailz(particle)+1
!
!                   high1=i
!
!               end if
!               if (hole/= 0_8) then
!                   k=trailz(hole)+1
!
!                   low1=k
!
!               end if
!
!                hole= iand(diff2,Nalphbet(2,c1))
!                particle = iand(diff2,Nalphbet(2,c2))
!                if (particle/= 0_8) then
!                  j=trailz(particle)+1
!
!                   high2=j
!
!                end if
!
!               if (hole/= 0_8) then
!                   l=trailz(hole)+1
!
!                   low1=l
!
!               end if
!
!                do spin11=1,2
!                 buffer_prime=Nalphbet(spin11,c1)
!                 buffer_prime_2=Nalphbet(spin11,c2)
!
!                  do while (buffer_prime/=0_8)
!                      tz=trailz(buffer_prime)
!                      tz2=trailz(buffer_prime_2)
!
!                      if (tz/=tz2 .and. btest(Nalphbet(spin11,c2),tz)) then
!
!                          ep=-ep
!                      end if
!                      buffer_prime=iand(buffer_prime, buffer_prime-1_8)
!                      buffer_prime_2=iand(buffer_prime_2, buffer_prime_2-1_8)
!                  end do
!
!                end do
!                c_term=2*civ1*civ2*ep
!                sorb=i !korb
!                qorb=j!lorb
!                porb=l !jorb
!                rorb=k !iorb
!                twordm(porb,rorb,sorb,qorb)=twordm(porb,rorb,sorb,qorb)&
!                    +c_term
!!all same spin so all swaps are allowed
!   !Case 2 Swap p and r introduces negative sign
!
!
!
!                twordm(rorb,porb,qorb,sorb)=twordm(rorb,porb,qorb,sorb)&
!                    +c_term
!
!
!
!
!
!
!            end if
!
!
!
!
!            end do

!            ndiff1=popcnt(tmp)
!            ndiff1=ishft(ndiff1,-1)
!
!            if (ndiff1==2) then
!
!
!               particle= iand(tmp,alphbet_2)
!               hole= iand(tmp,alphbet)
!
!               if (particle/= 0_8) then
!                   i=trailz(particle)+1
!                   particle=iand(particle,particle-1_8)
!                   j=trailz(particle)+1
!                   high1=i
!                   high2=j
!               end if
!               if (hole/= 0_8) then
!                   k=trailz(hole)+1
!                   hole=iand(hole,hole-1_8)
!                   l=trailz(hole)+1
!                   low1=k
!                   low2=l
!               end if
!
!
!                buffer=alphbet
!                buffer_2=alphbet_2
!
!                do while (buffer/=0_8)
!                  tz=trailz(buffer)
!                  tz2=trailz(buffer_2)
!
!                  if (tz/=tz2 .and. btest(alphbet_2,tz)) then
!
!                      ep=-ep
!                  end if
!                  buffer=iand(buffer, buffer-1_8)
!                  buffer_2=iand(buffer_2, buffer_2-1_8)
!                 end do
!
!
!
!                spin11=mod(low1,2)
!                spin12=mod(low2,2)
!                spin21=mod(high1,2)
!                spin22=mod(high2,2)
!
!                i = nint((i) / 2.0_dp + 0.1)
!                j = nint((j) / 2.0_dp + 0.1)
!                k = nint((k) / 2.0_dp + 0.1)
!                l = nint((l) / 2.0_dp + 0.1)
!
!
!
!
!                sorb = i
!                qorb = j
!                porb = l
!                rorb = k
!
!
!                eg = 1.00
!
!                if (spin21 == spin11 .and. spin22 == spin12) then
!
!
!                    twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
!                                civ1 * civ2 * ep * eg
!
!                end if
!
!                sorb = i
!                qorb = j
!                rorb = l
!                porb = k
!
!                eg=-1.00
!
!
!                if (spin21 == spin12 .and. spin22 == spin11) then
!
!                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
!                                civ1 * civ2 * ep * eg
!
!
!                end if
!
!                qorb = i
!                sorb = j
!                porb = l
!                rorb = k
!
!                eg=-1.00
!
!                if (spin11 == spin22 .and. spin21 == spin12) then
!
!                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
!                                civ1 * civ2 * ep * eg
!
!
!
!                end if
!
!
!                qorb = i
!                sorb = j
!                rorb = l
!                porb = k
!
!                eg=1.00
!
!                if (spin12 == spin22 .and. spin21 == spin11) then
!
!                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
!                                civ1 * civ2 * ep * eg
!
!
!                end if
!
!
!            else if (ndiff1==1) then
!
!
!                if (Nalphbet(1,c1)==Nalphbet(1,c2)) then
!                            tmp = xor( Nalphbet(2,c1), Nalphbet(2,c2))
!                            particle = iand(tmp, Nalphbet(2,c2))
!                            hole = iand(tmp, Nalphbet(2,c1))
!                            if (particle /= 0_8) then
!                                tz = trailz(particle)
!                                j=tz+1
!                            end if
!
!                            if (hole /= 0_8) then
!                                tz = trailz(hole)
!                                i=tz+1
!                            end if
!                            low=i*2
!                            high=j*2
!                else
!                            tmp = xor( Nalphbet(1,c1), Nalphbet(1,c2))
!                            particle = iand(tmp, Nalphbet(1,c2))
!                            hole = iand(tmp, Nalphbet(1,c1))
!                            if (particle /= 0_8) then
!                                tz = trailz(particle)
!                                j=tz+1
!                            end if
!
!                            if (hole /= 0_8) then
!                                tz = trailz(hole)
!                                i=tz+1
!                            end if
!                             low=i*2-1
!                             high=j*2-1
!                end if
!
!               ! number=popcnt(ibits(alphbet,low,(high-low-1)))
!                      ! print*,'newnumber',number
!!                       red_vec=integer2binary_orbs_bit(Nalphbet(1,:,c2),maxnmo*2)
!!                       number=size(pack(red_vec(low+1:high-1),red_vec(low+1:high-1)/=0))
!!                        print*,'oldnumber',number
!
!
!                !ep=(-1)**number
!                buffer=alphbet
!                buffer_2=alphbet_2
!                do while (buffer/=0_8)
!                  tz=trailz(buffer)
!                  tz2=trailz(buffer_2)
!
!                  if (tz/=tz2 .and. btest(alphbet_2,tz)) then
!
!                      ep=-ep
!                  end if
!                  buffer=iand(buffer, buffer-1_8)
!                  buffer_2=iand(buffer_2, buffer_2-1_8)
!                 end do
!
!                c = ep*civ1*civ2
!                count=1
!                do tz=1,norbs*2
!
!                    if (btest(alphbet_2, tz-1) .and. btest(alphbet,tz-1)) then
!                        occs(count)=nint((tz) / 2.0_dp + 0.1)
!                        spin(count)=mod(tz,2)
!                        count=count+1
!                    end if
!                end do
!
!                if (allocated(occs_r)) deallocate(occs_r)
!                allocate(occs_r(count-1))
!                occs_r=occs(1:count-1)
!
!                spin11=mod(low,2)
!                spin22=mod(low,2)
!
!
!                do tz=1,size(occs_r)
!                        qorb = j
!                        porb = i
!                        eg = 1.00
!
!
!!twordm(18,18,18,15)
!                        sorb=occs_r(tz)
!                        rorb=occs_r(tz)
!
!                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
!                                c* eg
!
!
!
!                        qorb = j
!                        rorb = i
!                        eg = -1.00
!                        sorb=occs_r(tz)
!                        porb=occs_r(tz)
!                        if (spin(tz)==spin11 .and. spin(tz)==spin22) then
!                               twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
!                                c * eg
!
!                        end if
!
!
!                        sorb = j
!                        porb = i
!                        eg = -1.00
!                        qorb=occs_r(tz)
!                        rorb=occs_r(tz)
!                        if (spin11==spin(tz)) then
!                                twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
!                                c * eg
!                        endif
!
!
!                        sorb = j
!                        rorb = i
!                        eg = 1.00
!                        porb=occs_r(tz)
!                        qorb=occs_r(tz)
!
!                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
!                                c * eg
!
!
!
!                end do
!
!
!            end if
!
!
!
!
!
!
!
!        end do
!
!
!
!    end do


!
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
    end subroutine createtwordm_bit











       subroutine one_rdm_bit(file_read,onerdm, maxnmo,numberlines,newdat,irep)

           Use, intrinsic :: iso_fortran_env, Only : iostat_end
           implicit none

           INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
           INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)


           integer(kind=ikind):: numberlines
           character(len=20),intent(in) :: file_read
           integer*8,intent(in),dimension(:):: irep
           real(kind=dp), dimension(:), allocatable :: civs
           integer*8,dimension(:,:,:), allocatable :: Nalphbet
           integer*8:: j, i, error_1,count,temp,sp,diff1,diff2,popcnt1,popcnt2,compnum,buffer_prime, buffer_prime_2
            integer*8:: particle_int,hole_int,k
           real(kind=dp),dimension(:,:), allocatable :: onerdm
           integer(kind=ikind),intent(inout):: maxnmo
           integer*8:: exc(0:2,2,2),deg,buffer,c1,c2,ep,n1a,n2a,n1b,n2b,low,high,number,red_vec(maxnmo*2)
           double precision :: phase, c,civ1,civ2
           integer*8:: hole, particle, tmp,tz,final,max,prueba(maxnmo*2),buffer_2,tz2,index,nn,nel,n
           real(kind=dp)::time1,time2
           integer*8, intent(in), dimension(:,:):: newdat
           integer*8, allocatable, dimension(:,:,:) :: orb_mat
           integer*8,allocatable, dimension(:) :: orb_indexed
           integer*8, dimension(numberlines):: sum_mat



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

           maxnmo=maxval(integer2binary_orbs(maxval(Nalphbet(1,1,:))))
           nel=popcnt(Nalphbet(1,1,1))
            allocate(onerdm(maxnmo,maxnmo))

           onerdm=0.0_dp

           allocate(orb_mat(2,nel,size(newdat(:,2))), orb_indexed(nel))
           print*, 'allocated orb_mat, orb_indexed'
           do i=1,size(newdat(:,1))
               count=1
               buffer=newdat(i,1)
                do while (buffer /= 0_8)
                                j = trailz(buffer)
                                orb_mat(1,count,i)=j+1
                                count=count+1
                                buffer = iand(buffer,buffer-1_8)
                end do
           end do
           do i=1,size(newdat(:,2))
               count=1
               buffer=newdat(i,2)
                do while (buffer /= 0_8)
                                j = trailz(buffer)
                                orb_mat(2,count,i)=j+1
                                count=count+1
                                buffer = iand(buffer,buffer-1_8)
                end do
           end do

            print*, 'orbmat calculated',numberlines,maxnmo

            compnum=0

            do i=1,maxnmo
                compnum=compnum+2**(i-1)
           end do
            print*,compnum
            call cpu_time(time1)

            do c1=1,numberlines


                   if (c1==1000) then
                     print*,'Almost there 100000'
                       call cpu_time(time2)
                       print*,time2-time1
                   end if

                    civ1=civs(c1)



                        do i=1,2
                            buffer = Nalphbet(1,i,c1)
!                            do n=1,size(newdat(:,i))
!                                if (newdat(n,i)==buffer) then
!                                    index=n
!                                    exit
!                                end if
!                            end do
                            index=FINDLOC(newdat(:,i),buffer,1)

                            orb_indexed=orb_mat(i,:,index)
                            do nn=1,size(orb_indexed)
                                j=orb_indexed(nn)
                                onerdm(j,j) = onerdm(j,j) + &
                                civ1*civ1
                            end do
!                            do while (buffer /= 0_8)
!                                j = trailz(buffer)
!
!                                onerdm(j+1,j+1) = onerdm(j+1,j+1) + &
!                                civ1*civ1
!
!                                buffer = iand(buffer,buffer-1_8)
!                            end do


                        end do


                   ! print*,'the fun starts'
                    buffer=xor(Nalphbet(1,1,c1), compnum)
                   ! print*,'buffer 1 '

                  !  print*,'buffer 2 '
                    do while (buffer /= 0_8)
                        j=trailz(buffer)
                        hole_int=ibset(Nalphbet(1,1,c1),j)
                        buffer_2=Nalphbet(1,1,c1)
                        do while(buffer_2 /=0_8)
                            k=trailz(buffer_2)
                            particle_int = ibclr(hole_int,k)
                            buffer_2 = iand(buffer_2,buffer_2-1_8)

                            !print*,'before index'

                           ! print*,particle_int,Nalphbet(1,2,c1)
                           ! print*,size(abs(Nalphbet(1,:,:)- spread( (/particle_int, Nalphbet(1,2,c1)/), 2, size(Nalphbet(1,1,:)))))

!                            n=1
!                            count=0
!                            do while(n/=0)
!                                count=count+1
!                                if (count>size(Nalphbet(1,1,:))) exit
!                                n=sum(abs(Nalphbet(1,:,count)-(/particle_int, Nalphbet(1,2,c1)/)))
!
!                            end do
!                            if (n/=0) cycle
                            !index=count
                             !call sgemm('t','n', size(Nalphbet(1,1,:)),1 , 2, 1.0, Nalphbet(1,:,:), &
                              !&           2, (/-Nalphbet(1,2,c1), particle_int/), 2, 0.0_dp, sum_mat,size(Nalphbet(1,1,:)))
                            index=0
                            civ2=0.0_dp
                            do i = 1,size(Nalphbet(1,1,:))
                                if (particle_int==Nalphbet(1,1,i) .and. Nalphbet(1,2,c1)==Nalphbet(1,2,i)) then
                                    index = i
                                     civ2=civs(index)
                                    exit
                                endif
                            end do
                            if (index==0) cycle
                           !     print*,particle_int
                           ! end if
                           ! index=minloc(abs(sum_mat),1)

                            !index = minloc(sum(abs(Nalphbet(1,:,:)- spread( (/particle_int, Nalphbet(1,2,c1)/), 2, size(Nalphbet(1,1,:)))), 1),1)
                           ! print*,index
                          !  if (0==index .or. index>size(Nalphbet(1,1,:))) then
                           !     cycle
                          !  end if
                          !  print*,index,particle_int, Nalphbet(1,2,c1),Nalphbet(1,1,index), Nalphbet(1,2,index)


                            buffer_prime=Nalphbet(1,1,c1)
                            buffer_prime_2=Nalphbet(1,1,index)
                            ep=1
                            do while (buffer_prime/=0_8)
                                tz=trailz(buffer_prime)
                                tz2=trailz(buffer_prime_2)

                                if (tz/=tz2 .and. btest(Nalphbet(1,1,index),tz)) then

                                    ep=-ep
                                end if
                                buffer_prime=iand(buffer_prime, buffer_prime-1_8)
                                buffer_prime_2=iand(buffer_prime_2, buffer_2-1_8)
                            end do

                            phase=ep

                           ! print*,j,k,index,phase
                            c = phase*civ1*civ2
                            onerdm(j+1,k+1) = onerdm(j+1,k+1) + c


                        end do
                        buffer = iand(buffer,buffer-1_8)
                    end do

                    buffer=xor(Nalphbet(1,2,c1), compnum)

                    do while (buffer /= 0_8)
                        j=trailz(buffer)
                        hole_int=ibset(Nalphbet(1,2,c1),j)
                        buffer_2=Nalphbet(1,2,c1)
                        do while(buffer_2 /=0_8)
                            k=trailz(buffer)
                            particle_int=ibclr(hole_int,k)
                            buffer_2 = iand(buffer_2,buffer_2-1_8)
                             do i = 1,size(Nalphbet(1,1,:))
                                if (particle_int==Nalphbet(1,2,i) .and. Nalphbet(1,1,c1)==Nalphbet(1,1,i)) then
                                    index = i
                                    exit
                                endif
                            end do
                           ! index = minloc(sum(abs(Nalphbet(1,:,:)- spread( (/Nalphbet(1,1,c1), particle_int/), 2, size(Nalphbet(1,1,:)))), dim=1),1)
                            civ2=civs(index)

                            buffer_prime=Nalphbet(1,2,c1)
                            buffer_prime_2=Nalphbet(1,2,index)
                            ep=1
                            do while (buffer_prime/=0_8)
                                tz=trailz(buffer_prime)
                                tz2=trailz(buffer_prime_2)

                                if (tz/=tz2 .and. btest(Nalphbet(1,2,index),tz)) then

                                    ep=-ep
                                end if
                                buffer_prime=iand(buffer_prime, buffer_prime-1_8)
                                buffer_prime_2=iand(buffer_prime_2, buffer_2-1_8)
                            end do

                            phase=ep
                            c = phase*civ1*civ2
                            onerdm(j+1,k+1) = onerdm(j+1,k+1) + c


                        end do
                        buffer = iand(buffer,buffer-1_8)
                    end do

            end do
                 !   print*, 'number completed'


!                    do c2=c1+1,numberlines
!                        diff1=xor( Nalphbet(1,1,c1), Nalphbet(1,1,c2))
!                        diff2=xor(Nalphbet(1,2,c1), Nalphbet(1,2,c2))
!
!                        popcnt1=popcnt(diff1)
!                        popcnt2= popcnt(diff2)
!                        !number=excitations(,Nalphbet(1,:,c2),1)
!                        number= popcnt1+popcnt2
!
!                        number=number-1
!
!                    if  (number== 1) then
!
!                        if (popcnt1==0) then
!                            tmp =diff2
!                            particle = iand(tmp, Nalphbet(1,2,c2))
!                            hole = iand(tmp, Nalphbet(1,2,c1))
!
!                            tz = trailz(particle)
!                            j=tz+1
!                            tz = trailz(hole)
!                            i=tz+1
!
!                            low=i*2
!                            high=j*2
!                        else
!                            tmp = xor( Nalphbet(1,1,c1), Nalphbet(1,1,c2))
!                            particle = iand(tmp, Nalphbet(1,1,c2))
!                            hole = iand(tmp, Nalphbet(1,1,c1))
!
!                            tz = trailz(particle)
!                            j=tz+1
!
!
!
!                            tz = trailz(hole)
!                            i=tz+1
!
!                            low=i*2-1
!                            high=j*2-1
!
!
!
!                        end if
!                        if (low>high) then
!                          temp=low
!                          low=high
!                          high=temp
!                      end if
!                        phase=1
!                        civ2=civs(c2)
!
!
!
!
!                      ! call combine_alpha_beta(Nalphbet(1,1,c2), Nalphbet(1,2,c2), final)
!
!                      ! number=popcnt(ibits(final,low,(high-low-1)))
!                      ! print*,'newnumber',number
!!                       red_vec=integer2binary_orbs_bit(Nalphbet(1,:,c2),maxnmo*2)
!!                       number=size(pack(red_vec(low+1:high-1),red_vec(low+1:high-1)/=0))
!!                        print*,'oldnumber',number
!
!
!
!                        do i=1,2
!                            buffer=Nalphbet(1,i,c1)
!                            buffer_2=Nalphbet(1,i,c2)
!
!                            do while (buffer/=0_8)
!                                tz=trailz(buffer)
!                                tz2=trailz(buffer_2)
!
!                                if (tz/=tz2 .and. btest(Nalphbet(1,i,c2),tz)) then
!
!                                    ep=-ep
!                                end if
!                                buffer=iand(buffer, buffer-1_8)
!                                buffer_2=iand(buffer_2, buffer_2-1_8)
!                            end do
!
!                        end do
!                       phase=ep
!                       c = phase*civ1*civ2
!
!                       onerdm(i,j) = onerdm(i,j) + c
!                       onerdm(j,i) = onerdm(j,i) + c
!
!
!                    end if
!                    end do

              !  end do
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
           integer*8, intent(in) :: det(Nint,2,Ndet)
           integer, intent(in) :: Ndet, Nint, mo_num
           double precision, intent(in) :: coef(Ndet)
           double precision, intent(out) :: density_matrix(mo_num,mo_num)
           integer :: i,j,k,l,ispin,ishift
           integer*8:: buffer
           integer :: deg
           integer :: exc(0:2,2,2)
           double precision :: phase, c
           integer :: n_excitations
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
           integer, intent(in) :: Nint
           integer*8, intent(in) :: det1(Nint,2), det2(Nint,2)
           integer, intent(out) :: exc(0:2,2,2)
           double precision, intent(out) :: phase
           integer :: l, ispin, idx_hole, idx_particle, ishift
           integer :: i,j,k,m,n,high, low,a,b,c,d,nperm,tz,nexc
           integer*8 :: hole, particle, tmp
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
           integer, intent(in) :: Nint
           integer*8, intent(in) :: det1(Nint,2)
           integer*8, intent(in) :: det2(Nint,2)
           integer, intent(out) :: exc(0:2,2,2)
           double precision, intent(out) :: phase
           integer :: tz, l, ispin, ishift, nperm, i, j, k, m, n, high, low
           integer*8 :: hole, particle, tmp
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
           integer, intent(in) :: Nint
           integer*8, intent(in) :: det1(Nint,2), det2(Nint,2)
           integer, intent(out) :: exc(0:2,2,2)
           integer, intent(out) :: degree
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

       recursive function combo(n,k) result(cmb)
           implicit none
           integer*8 :: cmb
           integer*8, intent(in) :: n,k
           integer*8 :: mm(100,100)
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