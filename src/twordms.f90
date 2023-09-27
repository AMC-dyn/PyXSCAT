module twordms
    implicit none

    contains



subroutine maxcoincidence(confs, ep2,ndiff)

        implicit none
                INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=ikind), intent(in), dimension(:,:) :: confs
        integer(kind=ikind), intent(out), dimension(:,:), allocatable :: ep2, ndiff
        integer(kind=ikind), dimension(size(confs(:,1)), size(confs(1,:))):: matdum
        integer(kind=ikind), dimension(:,:), allocatable :: mat1
        integer(kind=ikind) :: i,j, c1,c2, count,itemp


        allocate(ep2(size(confs(:,1)),size(confs(:,1))))
        ep2=1
        count=0
        matdum=0
        print*,'holaa'
            do i=1,size(confs(:,1))
                count=0
                do j=1,size(confs(1,:))


                    if (confs(i,j)/=0) then
                        count=count+1
                        matdum(i,count)=j
                        end if
                end do
            end do
           print*, 'matdum constructed'
            allocate(mat1(size(confs(:,1)),count), ndiff(size(confs(:,1)),size(confs(:,1))))
            ndiff=0
            mat1=matdum(:,1:count)
            print*, 'mat1 constructed 2'

            do c1=1, size(confs(:,1))
                do c2=c1+1,size(confs(:,1))

                    ep2(c1,c2)=1
                    do i=1,size(mat1(c1,:))
                        if  (mat1(c1,i) /= mat1(c2,i)) then


                            do j=1,size(mat1(c2,:))
                                if (mat1(c1,i) == mat1(c2,j)) then

                                    ep2(c1,c2)=-ep2(c1,c2)
                                    ep2(c2,c1)=ep2(c1,c2)
!                                    if (c1==1 .and. c2==85) print*,'changed',mat1(c2,i),mat1(c2,j)
!                                    itemp=mat1(c2,j)
!                                    mat1(c2,j)=mat1(c2,i)
!                                    mat1(c2,i)=itemp

                                    exit




                                end if
                            end do
                        end if
                    end do

                    do j=1,size(confs(1,:))
                        if (confs(c1,j)/=confs(c2,j)) then
                           ndiff(c1,c2)= ndiff(c1,c2) +1
                           ndiff(c2,c1)= ndiff(c2,c1) +1
                        end if
                    end do

                end do
            end do


    end subroutine maxcoincidence

    subroutine maxc_individual(c1,c2, ep,diffs1,diffs2,spin1,spin2,ndiff)

        implicit none
                INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=ikind), intent(in), dimension(:) :: c1,c2
        integer(kind=ikind), intent(out):: ep, ndiff
        integer(kind=ikind), intent(out), dimension(:), allocatable :: diffs1,diffs2,spin1,spin2
        integer(kind=ikind), dimension(:,:), allocatable :: mat1
        integer(kind=ikind) :: i,j, count1,itemp


        allocate(mat1(2,count(c1/=0)))
        ep=1
        count1=0


            do i=1,size(c1)



                    if (c1(i)/=0) then
                        count1=count1+1
                        mat1(1,count1)=i
                        end if
                end do

         count1=0
        do i=1,size(c2)



                    if (c2(i)/=0) then
                        count1=count1+1
                        mat1(2,count1)=i
                        end if
        end do


            ndiff=0


            ep=1
            do i=1,size(mat1(1,:))

                if  (mat1(1,i) /= mat1(2,i)) then


                    do j=1,size(mat1(1,:))
                        if (mat1(1,i) == mat1(2,j)) then

                                    !ndiff=ndiff-1
                                    ep=-ep

                                    itemp=mat1(2,j)
                                    mat1(2,j)=mat1(2,i)
                                    mat1(2,i)=itemp

                                    exit




                                end if
                            end do
                        end if
                    end do
                 do i=1,size(mat1(1,:))
                     if  (mat1(1,i) /= mat1(2,i)) then
                         ndiff=ndiff+1

                         end if
                 end do

                allocate(diffs1(ndiff), diffs2(ndiff), spin1(ndiff), spin2(ndiff))
                count1=0
               do i=1,size(mat1(1,:))
                    if  (mat1(1,i) /= mat1(2,i)) then
                        count1=count1+1
                        diffs1(count1)=nint(mat1(1,i) / 2.0_dp + 0.1)
                        diffs2(count1)=nint(mat1(2,i) / 2.0_dp + 0.1)
                        spin1(count1)=mod(mat1(1,i),2)
                        spin2(count1)=mod(mat1(2,i),2)
                   end if
               end do




    end subroutine maxc_individual


    subroutine createtwordm(confs,civs,ndiff,ep2,mat,total,state1,state2)

    implicit none

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

    integer(kind=ikind), intent(in),dimension(:,:) :: confs
    integer(kind=ikind), intent(in), dimension(:,:) :: ep2, ndiff
    integer(kind=ikind),intent(in):: state1,state2
    real(kind=dp), intent(in), dimension(:,:) :: civs
    real(kind=dp), intent(out), dimension(:), allocatable :: total
    integer(kind=ikind), intent(out), dimension(:,:), allocatable :: mat

    real(kind=dp), dimension(:,:,:,:), allocatable :: twordm
    integer(kind=ikind) :: ep,nc1,lconfs,norbs,nc2,sorb,rorb,qorb,porb,p,q,r,s,c1,c2,count1,count

    integer(kind=ikind) :: i,i1,i2,n,count2,eg, ndiff1,ndiff2

    integer(kind=ikind), dimension(:), allocatable :: mat1,mat2
    logical(4) :: sdef, rdef, pdef, qdef
    logical(4), dimension(:), allocatable :: logicaltwordms
    real(kind=dp) :: cutoff,temp
    integer(kind=ikind), dimension(:), allocatable :: diffs1,diffs2
    integer(kind=ikind) :: spins,spinr,spinq, spinp
    integer(kind=ikind), dimension(:), allocatable :: spin1,spin2
    integer(kind=ikind),  dimension(:,:), allocatable :: matdum
    real(kind=dp), dimension(:), allocatable :: totaldum

    lconfs=size(confs(1,:))
    norbs=lconfs/2
    allocate(twordm(norbs,norbs,norbs,norbs))
    twordm=0.0_dp
    nc1=1
    nc2=1
  !  ep=1
    print*, 'starting the twordm'
    do c1=1,size(confs(:,1))
        do c2=1,size(confs(:,1))
            temp=twordm(11,11,12,14)

            ndiff1=ndiff(c1,c2)


            if (ndiff(c1,c2)/=0 .and. ndiff1 <=4) then

                sdef = .False.
                rdef = .False.
                pdef = .False.
                qdef = .False.
!                if (allocated(diffs1)) deallocate(diffs1)
!                if (allocated(diffs2)) deallocate(diffs2)
!                if (allocated(spin1)) deallocate(spin1)
!                if (allocated(spin2)) deallocate(spin2)
!                allocate(diffs1(ndiff(c1,c2)/2), diffs2(ndiff(c1,c2)/2), spin1(ndiff(c1,c2)/2))
!                allocate(spin2(ndiff(c1,c2)/2))
!
!                count1=1
!                count2=1
!
!                do n=1,size(confs(c1,:))
!                    if (confs(c1,n) /= confs(c2,n)) then
!                        if (confs(c1,n) /= 0) THEN
!                            diffs1(count1)=nint((n) / 2.0_dp + 0.1)
!                            spin1(count1)=confs(c1,n)
!                            count1=count1+1
!
!
!                        elseif (confs(c2,n) /= 0) THEN
!                            diffs2(count2)=nint((n) / 2.0_dp + 0.1)
!                            spin2(count2)=confs(c2,n)
!                            count2=count2+1
!                        end if
!                    end if
!                enddo
!!
!                if (c1==1 .and. c2==85) print*, spin1
!                 if (c1==1 .and. c2==85) print*, spin2

                if (ndiff(c1,c2) == 4) then
                      call maxc_individual(confs(c1,:), confs(c2,:), ep,diffs1, diffs2, spin1,spin2,ndiff2)

                      !if (c1==1 .and. c2==85) print*,spin1,spin2,ep,ep2(c1,c2),count1

                    sorb = diffs2(1)
                    qorb = diffs2(2)
                    porb = diffs1(2)
                    rorb = diffs1(1)

                    eg = 1.00

                    if (spin2(1) == spin1(1) .and. spin2(2) == spin1(2)) THEN
                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
                                civs(c1,state1) * civs(c2,state2) * ep * eg



                    end if

                    sorb = diffs2(1)
                    qorb = diffs2(2)
                    rorb = diffs1(2)
                    porb = diffs1(1)
                    eg = -1.00

                    if (spin2(1) == spin1(2) .and. spin2(2) == spin1(1)) then
                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )&
+ civs(c1,state1) * civs(c2,state2) * ep * eg

                    endif
                    qorb = diffs2(1)
                    sorb = diffs2(2)
                    porb = diffs1(2)
                    rorb = diffs1(1)
                    eg = -1.00

                    if (spin1(1) == spin2(2) .and. spin2(1) == spin1(2) ) THEN
                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
                                civs(c1,state1) * civs(c2,state2) * ep * eg



                    end if

                    qorb = diffs2(1)
                    sorb = diffs2(2)
                    rorb = diffs1(2)
                    porb = diffs1(1)
                    eg = 1.00

                    if (spin1(2) == spin2(2) .and. spin2(1) == spin1(1)) then
                        twordm(porb , rorb , sorb , qorb )= twordm(porb , rorb , sorb , qorb ) &
                                + civs(c1,state1) * civs(c2,state2) * ep * eg

                    end if


                elseif (ndiff(c1,c2) == 2) THEN
                    call maxc_individual(confs(c1,:), confs(c2,:), ep,diffs1, diffs2, spin1,spin2,ndiff1)

                    qorb = diffs2(1)
                    porb = diffs1(1)
                    eg = 1.00

                    do i=1,size(confs(c2,:))

                        if (confs(c2,i) /=0) then
                            sdef = .True.

                            sorb =nint(i / 2.0_dp + 0.1)

                            spins = mod(i,2)
                        end if

                        if (confs(c1,i) /= 0) then

                            rdef = .True.
                            rorb = nint((i) / 2.0_dp + 0.1)
                            spinr = mod(i,2)
                        end if

                        if (sdef .and. rdef .and. spins == spinr .and. spin2(1) == spin1(1)) then
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.

                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) &
                                    + civs(c1,state1) * civs(c2,state2) * ep * eg







                        else
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.

                        end if
                    enddo
                    qorb = diffs2(1)
                    rorb = diffs1(1)
                    eg = -1.00

                    do  i= 1, size(confs(c2,:))

                        if (confs(c2,i) /= 0) then
                            sdef = .True.
                            sorb = nint((i) / 2.0_dp + 0.1)
                            spins =mod(i,2)
                        endif
                        if (confs(c1,i) /= 0) then
                            pdef = .True.
                            porb = nint((i) / 2.0_dp + 0.1)
                            spinp = mod(i,2)
                        end if

                        if (sdef .and. pdef .and. spin1(1) == spins .and. spin2(1) == spinp) then
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.
                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) + &
                                    civs(c1,state1) * civs(c2,state2) * ep * eg


                        else
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.

                        end if
                    end do



                    sorb = diffs2(1)
                    porb = diffs1(1)
                    eg = -1.00

                    do i=1,size(confs(c2,:))

                        if (confs(c2,i) /= 0) then
                            qdef = .True.
                            qorb = nint((i) / 2.0_dp + 0.1)
                            spinq = mod(i,2)
                        endif
                        if (confs(c1,i) /= 0) then
                            rdef = .True.
                            rorb = nint((i) / 2.0_dp + 0.1)
                            spinr = mod(i,2)
                        endif
                        if (rdef .and. qdef .and. spin1(1) == spinr .and. spin2(1) == spinq) then
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.
                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) +&
                                    civs(c1,state1) * civs(c2,state2) * ep * eg



                        else
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.

                        end if
                    enddo

                    sorb = diffs2(1)
                    rorb = diffs1(1)
                    eg = 1.00

                    do i=1,size(confs(c2,:))
                        if (confs(c2,i) /= 0) then
                            qdef = .True.
                            qorb = nint((i) / 2.0_dp + 0.1)
                            spinq = mod(i,2)

                        end if

                        if (confs(c1,i) /= 0) then
                            pdef = .True.
                            porb = nint((i) / 2.0_dp + 0.1)
                            spinp = mod(i,2)

                        end if

                        if (qdef .and. pdef .and. spinq == spinp .and. spin2(1) == spin1(1)) then
                            qdef = .False.
                            pdef = .False.

                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) &
                                    + civs(c1,state1) * civs(c2,state2) * ep * eg



                        else

                            pdef = .False.
                            qdef = .False.

                        end if
                    enddo
                endif

            elseif (ndiff(c1,c2) == 0) then

                ep = 1

                do i1=1,size(confs(c1,:))

                    do i2=1,size(confs(c2,:))

                        if (i1/=i2) then

                            if (confs(c1,i1) /= 0 .and. confs(c2,i2) /= 0) then

                                sorb = nint(i1 / 2.0_dp + 0.1)
                                qorb = nint(i2 / 2.0_dp + 0.1)

                                if (confs(c1,i1) == confs(c1,i2)) then
                                    porb = sorb
                                    rorb = qorb

                                    eg = -1.00

                                    twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) &
                                            + civs(c1,state1) * civs(c2,state2) * ep * eg


                                end if
                                porb = qorb
                                rorb = sorb
                                eg = 1.00

                                twordm(porb , rorb , sorb , qorb )  = twordm(porb , rorb , sorb , qorb ) &
                                        + civs(c1,state1) * civs(c2,state2) * ep * eg
                            endif
                        end if

                    enddo
                end do
            end if

            if (twordm(11,11,12,14) /= temp) print*,porb,rorb,sorb,qorb,twordm(11,11,12,14)
        end do

    end do



    cutoff = 1E-30
    count=0
    count2=1
    print*,'final value', twordm(11,11,12,14)
    allocate(logicaltwordms(norbs**4))
    allocate(totaldum(norbs**4), matdum(norbs**4,4))
    logicaltwordms(:)=.False.
    open (unit = 15, file = 'twordm_fortran.dat')
    do p=1,norbs
        do q=1,norbs
            do r=1,norbs
                do s=1,norbs
                    totaldum(count2)=twordm(p,q,r,s)
                    matdum(count2,:)=(/p,s,q,r/)
                    if (abs(twordm(p,q,r,s))>=cutoff) then
                        count=count+1
                        logicaltwordms(count2)=.True.
                        write(15,*) p, s, q, r, twordm(p,q,r,s)

                    end if
                    count2=count2+1
                end do
            end do
        end do
    end do
    close(15)

    allocate(mat(count, 4), total(count))
    count=1
    do i=1,count2-1
        if (logicaltwordms(i)) then
            mat(count,:)=matdum(i,:)
            total(count)=totaldum(i)

            count=count+1
        end if
    end do

        print*, 'twordm calculated'



    end subroutine createtwordm






     subroutine createtwordm_slow(file_read,mat,total)


     Use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
     character(len=30), intent(in) :: file_read
     integer :: norbs


    real(kind=dp), intent(out), dimension(:), allocatable :: total
    integer(kind=ikind), intent(out), dimension(:,:), allocatable :: mat

    real(kind=dp), dimension(:,:,:,:), allocatable :: twordm
    integer(kind=ikind) :: j,ep,nc1,lconfs,nc2,sorb,rorb,qorb,porb,p,q,r,s,c1,c2,count1,count

    integer(kind=ikind) :: i,i1,i2,n,count2,eg, ndiff1,N1a,N1b,N2a,N2b

    integer*8, dimension(:), allocatable :: mat1,mat2,Nalpha,Nbeta
    logical(4) :: sdef, rdef, pdef, qdef
    logical(4), dimension(:), allocatable :: logicaltwordms
    real(kind=dp) :: cutoff,civ1,civ2
    integer(kind=ikind), dimension(:), allocatable :: diffs1,diffs2,conf1,conf2
    integer(kind=ikind) :: spins,spinr,spinq, spinp,error_1,count_zeros_1,count_zeros_2
    integer(kind=ikind), dimension(:), allocatable :: spin1,spin2,c1a,c2a,c1b,c2b
    integer(kind=ikind),  dimension(:,:), allocatable :: matdum
    real(kind=dp), dimension(:), allocatable :: totaldum,civs




    open(15,file=file_read)

    count=0
    do while (error_1 /= -1)
        read(15,*, iostat = error_1) j
        count=count+1
        print*,error_1
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

     norbs=maxval(integer2binary_orbs(maxval(Nalpha)))
     print*,'Number of orbitals',norbs
    print*, 'read file'
    allocate(twordm(norbs,norbs,norbs,norbs))
     print*, 'allocataed twordm'
    twordm=0.0_dp
    nc1=1
    nc2=1
    ep=1

    do c1=1,count-1
        do c2=1,count-1
            ndiff1=0
            ep=1
            eg=1

            N1a=Nalpha(c1)
            N2a=Nalpha(c2)
            N1b=Nbeta(c1)
            N2b=Nbeta(c2)
            civ1=civs(c1)
            civ2=civs(c2)


            ep=1




                    c1a=integer2binary(N1a,norbs)
                    c2a=integer2binary(N2a,norbs)
                    c1b=integer2binary(N1b,norbs)
                    c2b=integer2binary(N2b,norbs)

                    if (allocated(conf1)) deallocate(conf1)
                    if (allocated(conf2)) deallocate(conf2)
                    allocate(conf1(size(c1a)+size(c1b)),conf2(size(c2a)+size(c2b)))

                    do i=0,size(c1a)-1

                             conf1(2 * i + 1) = c1a(i+1)
                             conf1( 2 * i + 2) = 2*c1b(i+1)
                             conf2( 2 * i + 1) = c2a(i+1)
                             conf2( 2 * i + 2)= 2*c2b(i+1)
                    end do
                    count_zeros_1=0
                    count_zeros_2=0
                    do j=1,size(conf1)
                        if (conf1(j) /= 0) then
                                count_zeros_1=count_zeros_1+1
                        end if
                        if (conf2(j) /= 0) then
                                count_zeros_2=count_zeros_2+1
                        end if
                        if (conf1(j)/=conf2(j)) then
                           ndiff1= ndiff1 +1
                           do i=1,size(conf2)
                            if (conf1(j)/=conf2(i)) then
                                ep=-ep
                            end if
                        end do
                    end if
                    enddo
            if (count_zeros_1/=count_zeros_2) then
                    print*,'ERROR',count_zeros_1, count_zeros_2

            end if


            if (ndiff1/=0 .and. ndiff1 <=4) then

                sdef = .False.
                rdef = .False.
                pdef = .False.
                qdef = .False.
                if (allocated(diffs1)) deallocate(diffs1)
                if (allocated(diffs2)) deallocate(diffs2)
                if (allocated(spin1)) deallocate(spin1)
                if (allocated(spin2)) deallocate(spin2)
                allocate(diffs1(ndiff1/2), diffs2(ndiff1/2), spin1(ndiff1/2))
                allocate(spin2(ndiff1/2))

                count1=1
                count2=1

                do n=1,size(conf1)
                    if (conf1(n) /= conf2(n)) then
                        if (conf1(n) /= 0) THEN
                            diffs1(count1)=nint((n) / 2.0_dp + 0.1)
                            spin1(count1)=conf1(n)
                            count1=count1+1


                        elseif (conf2(n) /= 0) THEN
                            diffs2(count2)=nint((n) / 2.0_dp + 0.1)
                            spin2(count2)=conf2(n)
                            count2=count2+1
                        end if
                    end if
                enddo

                if (ndiff1 == 4) then

                    sorb = diffs2(1)
                    qorb = diffs2(2)
                    porb = diffs1(2)
                    rorb = diffs1(1)

                    eg = 1.00

                    if (spin2(1) == spin1(1) .and. spin2(2) == spin1(2)) THEN
                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
                                civ1 * civ2 * ep * eg

                    end if

                    sorb = diffs2(1)
                    qorb = diffs2(2)
                    rorb = diffs1(2)
                    porb = diffs1(1)
                    eg = -1.00

                    if (spin2(1) == spin1(2) .and. spin2(2) == spin1(1)) then
                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ civ1 * civ2 * ep * eg
                    endif
                    qorb = diffs2(1)
                    sorb = diffs2(2)
                    porb = diffs1(2)
                    rorb = diffs1(1)
                    eg = -1.00

                    if (spin1(1) == spin2(2) .and. spin2(1) == spin1(2) ) THEN
                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
                                civ1 * civ2 * ep * eg


                    end if

                    qorb = diffs2(1)
                    sorb = diffs2(2)
                    rorb = diffs1(2)
                    porb = diffs1(1)
                    eg = 1.00

                    if (spin1(2) == spin2(2) .and. spin2(1) == spin1(1)) then
                        twordm(porb , rorb , sorb , qorb )= twordm(porb , rorb , sorb , qorb ) &
                                + civ1 * civ2 * ep * eg

                    end if


                elseif (ndiff1 == 2) THEN

                    qorb = diffs2(1)
                    porb = diffs1(1)
                    eg = 1.00

                    do i=1,size(conf2)

                        if (conf2(i) /=0) then
                            sdef = .True.

                            sorb = nint((i) / 2.0_dp + 0.1)

                            spins = conf2(i)
                        end if

                        if (conf1(i) /= 0) then

                            rdef = .True.
                            rorb = nint((i) / 2.0_dp + 0.1)
                            spinr = conf1(i)
                        end if

                        if (sdef .and. rdef .and. spins == spinr .and. spin2(1) == spin1(1)) then
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.

                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) &
                                    + civ1 * civ2 * ep * eg



                        else
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.

                        end if
                    enddo
                    qorb = diffs2(1)
                    rorb = diffs1(1)
                    eg = -1.00

                    do  i= 1, size(conf2)

                        if (conf2(i) /= 0) then
                            sdef = .True.
                            sorb = nint((i) / 2.0_dp + 0.1)
                            spins = conf2(i)
                        endif
                        if (conf1(i) /= 0) then
                            pdef = .True.
                            porb = nint((i) / 2.0_dp + 0.1)
                            spinp = conf1(i)
                        end if

                        if (sdef .and. pdef .and. spin1(1) == spins .and. spin2(1) == spinp) then
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.
                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) + &
                                    civ1 * civ2 * ep * eg


                        else
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.

                        end if
                    end do



                    sorb = diffs2(1)
                    porb = diffs1(1)
                    eg = -1.00

                    do i=1,size(conf2)

                        if (conf2(i) /= 0) then
                            qdef = .True.
                            qorb = nint((i) / 2.0_dp + 0.1)
                            spinq = conf2(i)
                        endif
                        if (conf1(i) /= 0) then
                            rdef = .True.
                            rorb = nint((i) / 2.0_dp + 0.1)
                            spinr = conf1(i)
                        endif
                        if (rdef .and. qdef .and. spins == spinr .and. spin1(1) == spin2(1)) then
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.
                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) +&
                                    civ1 * civ2 * ep * eg


                        else
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.

                        end if
                    enddo

                    sorb = diffs2(1)
                    rorb = diffs1(1)
                    eg = 1.00

                    do i=1,size(conf2)
                        if (conf2(i) /= 0) then
                            qdef = .True.
                            qorb = nint((i) / 2.0_dp + 0.1)
                            spinq = conf2(i)

                        end if

                        if (conf1(i) /= 0) then
                            pdef = .True.
                            porb = nint((i) / 2.0_dp + 0.1)
                            spinp = conf1(i)

                        end if

                        if (qdef .and. pdef .and. spinq == spinp .and. spin2(1) == spin1(1)) then
                            qdef = .False.
                            pdef = .False.

                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) &
                                    + civ1 * civ2 * ep * eg

                        else

                            pdef = .False.
                            qdef = .False.

                        end if
                    enddo
                endif
              elseif (ndiff1== 0) then

                ep = 1

                do i1=1,size(conf1)

                    do i2=1,size(conf2)

                        if (i1/=i2) then

                            if (conf1(i1) /= 0 .and. conf2(i2) /= 0) then

                                sorb = nint(i1 / 2.0_dp + 0.1)
                                qorb = nint(i2 / 2.0_dp + 0.1)

                                if (conf1(i1) == conf2(i2)) then
                                    porb = sorb
                                    rorb = qorb

                                    eg = -1.00

                                    twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) &
                                            + civ1 * civ2 * ep * eg

                                    print*,'0 case', porb,rorb,sorb,qorb

                                end if
                                porb = qorb
                                rorb = sorb
                                eg = 1.00

                                twordm(porb , rorb , sorb , qorb )  = twordm(porb , rorb , sorb , qorb ) &
                                        + civ1 * civ2 * ep * eg
                            endif
                        end if

                    enddo
                end do
            end if

        end do
        print*,c1
    end do


    print*,'everything finishes', twordm(3,3,4,3)
    cutoff = 1E-30
    count=0
    count2=1
    print*,sum(twordm)
    allocate(logicaltwordms(norbs**4))
    allocate(totaldum(norbs**4), matdum(norbs**4,4))
    logicaltwordms(:)=.False.
    open (unit = 15, file = 'twordm_fortran.dat')
    do p=1,norbs
        do q=1,norbs
            do r=1,norbs
                do s=1,norbs
                    totaldum(count2)=twordm(p,q,r,s)
                    matdum(count2,:)=(/p,s,q,r/)
                    if (abs(twordm(p,q,r,s))>=cutoff) then
                        count=count+1
                        logicaltwordms(count2)=.True.
                        write(15,*) p, s, q, r, twordm(p,q,r,s)

                    end if
                    count2=count2+1
                end do
            end do
        end do
    end do
    close(15)

    allocate(mat(count, 4), total(count))
    count=1
    do i=1,count2-1
        if (logicaltwordms(i)) then
            mat(count,:)=matdum(i,:)
            total(count)=totaldum(i)

            count=count+1
        end if
    end do

        print*, 'twordm calculated'



    end subroutine createtwordm_slow






 function integer2binary(i,n) result(b)
    integer,intent(in) :: i,n
    integer :: b(n),count,count_orbs,b_orbs(n)
    integer, allocatable,dimension(:):: b_real,b_orbs_r
    integer k,j
    b=0
    j=i
    count=0
    b_orbs=0

    count_orbs=1
    do k=1,n
      b(k)=mod(j,2)
      if (mod(j,2)==1) then
          b_orbs(count_orbs)=count+1

          count_orbs=count_orbs+1
      end if
      j=j/2
      count=count+1

    enddo
    allocate(b_real(count),b_orbs_r(count_orbs-1))
    b_real=b(1:count)
    b_orbs_r=b_orbs(1:count_orbs-1)
  end function

        function integer2binary_orbs(i) result(b_orbs_r)
    integer,intent(in) :: i
    integer :: b(32),count,count_orbs,b_orbs(32)
    integer, allocatable,dimension(:):: b_real,b_orbs_r
    integer k,j
    b=0
    j=i
    count=0
    b_orbs=0

    count_orbs=1
    do while (j>=1.00)
      b(count+1)=mod(j,2)
      if (mod(j,2)==1) then
          b_orbs(count_orbs)=count+1

          count_orbs=count_orbs+1
      end if
      j=j/2
      count=count+1

    enddo
    allocate(b_real(count),b_orbs_r(count_orbs-1))
    b_real=b(1:count)
    b_orbs_r=b_orbs(1:count_orbs-1)
  end function

function integer2binary_orbs_bit(i,maxnmo) result(b_orbs)
    integer*8,intent(in) :: i(2),maxnmo
    integer*8 :: b_orbs(maxnmo),buffer,count
    integer, allocatable,dimension(:):: b_real,b_orbs_r
    integer k,j,n

    b_orbs=0


    do n=1,2
         count=0
        buffer = i(n)
        do k=1,maxnmo
          j = trailz(buffer)

          b_orbs(2*j+n)=1



          buffer = iand(buffer,buffer-1_8)

        end do

        end do



  end function

subroutine combine_alpha_beta(alpha,beta, result)
    implicit none
    integer, intent(in):: alpha, beta
    integer*8 :: buffer,z
    integer*8,intent(out) :: result

    result=0
    buffer=alpha
    do while (buffer/=0_8)
        z=trailz(buffer)

        result=ibset(result,2*z)

        buffer=iand(buffer,buffer-1_8)
    end do
    buffer=beta
    do while (buffer/=0_8)
        z=trailz(buffer)
        result=ibset(result,2*z+1)

        buffer=iand(buffer,buffer-1_8)
    end do




end subroutine combine_alpha_beta





end module twordms
