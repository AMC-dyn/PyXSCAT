
   module onerdm
       use twordms
       contains


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

    print*, 'onerdm calculated'



    end subroutine

   end module