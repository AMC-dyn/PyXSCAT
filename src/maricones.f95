!-----------------------------------------------------------------------
! 2RDM Fortran Modules for PyXSCAT code
! Andres Moreno Carrascosa and Mats Simmermacher, (2021)
!-----------------------------------------------------------------------

MODULE types

        implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

END MODULE types



MODULE uniquemodule

    implicit none 

    contains


    SUBROUTINE unique_integer(matrix,rmat,iuni,irec)
    
        use types
        
        integer(kind=ikind), dimension(:,:), intent(in)                  :: matrix
        integer(kind=ikind), dimension(:,:), allocatable                 :: sprmat
        integer(kind=ikind), dimension(:,:), intent(out), allocatable    :: rmat
        integer(kind=ikind), dimension(:), allocatable                   :: same, iunifull
        integer(kind=ikind), dimension(:), intent(out), allocatable      :: iuni
        integer(kind=ikind), dimension(:), intent(out), allocatable      :: irec
        integer(kind=ikind)                                              :: i, cnt

        allocate(sprmat(size(matrix(:,1)),size(matrix(1,:))))
        allocate(same(size(matrix(:,1))), iunifull(size(matrix(:,1))))
        allocate(irec(size(matrix(:,1))))
        
        irec = 0
        cnt = 0
            
        do i = 1, size(matrix(:,1))   
            if (irec(i) == 0) then
                cnt = cnt + 1
                iunifull(cnt) = i
                sprmat = spread(matrix(i,:), 1, size(matrix(:,i)))
                same = sum(abs(matrix - sprmat), dim=2)
                same = (-1)*same + 1
                same = (abs(same) + same)/2
                irec = irec + same*cnt
            endif
        enddo

        allocate(iuni(cnt))
        allocate(rmat(cnt,size(matrix(1,:))))

        iuni = iunifull(1:cnt)
        rmat = matrix(iuni,:)

    END SUBROUTINE


    SUBROUTINE unique_real(matrix,rmat,iuni,irec)
    
        use types
        
        real(kind=dp), dimension(:,:), intent(in)                      :: matrix
        real(kind=dp), dimension(:,:), allocatable                     :: sprmat
        real(kind=dp), dimension(:,:), intent(out), allocatable        :: rmat
        integer(kind=ikind), dimension(:), allocatable                 :: same, iunifull
        integer(kind=ikind), dimension(:), intent(out), allocatable    :: iuni
        integer(kind=ikind), dimension(:), intent(out), allocatable    :: irec
        integer(kind=ikind)                                            :: i, cnt

        allocate(sprmat(size(matrix(:,1)),size(matrix(1,:))))
        allocate(same(size(matrix(:,1))), iunifull(size(matrix(:,1))))
        allocate(irec(size(matrix(:,1))))
        
        irec = 0
        cnt = 0
            
        do i = 1, size(matrix(:,1))   
            if (irec(i) == 0) then
                cnt = cnt + 1
                iunifull(cnt) = i
                sprmat = spread(matrix(i,:), 1, size(matrix(:,i)))
                same = ceiling(sum(abs(matrix - sprmat), dim=2))
                same = (-1)*same + 1
                same = (abs(same) + same)/2
                irec = irec + same*cnt
            endif
        enddo

        allocate(iuni(cnt))
        allocate(rmat(cnt,size(matrix(1,:))))

        iuni = iunifull(1:cnt)
        rmat = matrix(iuni,:)

    END SUBROUTINE


    SUBROUTINE unique_total(matrix,total,rmat,rtot)
    
        use types
        
        integer(kind=ikind), dimension(:,:), intent(in)                  :: matrix
        real(kind=dp), dimension(:), intent(in)                          :: total
        integer(kind=ikind), dimension(:,:), allocatable                 :: sprmat
        integer(kind=ikind), dimension(:,:), intent(out), allocatable    :: rmat
        real(kind=dp), dimension(:), intent(out), allocatable            :: rtot
        integer(kind=ikind), dimension(:), allocatable                   :: same, iunifull
        real(kind=dp), dimension(:), allocatable                         :: rtotfull
        integer(kind=ikind), dimension(:), allocatable                   :: iuni
        integer(kind=ikind), dimension(:), allocatable                   :: irec
        integer(kind=ikind)                                              :: i, cnt, dim1, dim2


        dim1 = size(matrix(:,1))
        dim2 = size(matrix(1,:))
        print*,dim1,dim2
        allocate(sprmat(dim1,dim2))
        print*, 'first allocation'
        allocate(same(dim1), iunifull(dim1))
        print*,'second allocation'
        allocate(irec(dim1))
        print*,'third allocation'
        allocate(rtotfull(dim1))
        print*,'fourth allocation'
        print*,dim1
        irec = 0
        cnt = 0
            
        do i = 1, dim1
            if (irec(i) == 0) then
                cnt = cnt + 1
                iunifull(cnt) = i
                sprmat = spread(matrix(i,:), 1, size(matrix(:,i)))
                same = sum(abs(matrix - sprmat), dim=2)
                same = (-1)*same + 1
                same = (abs(same) + same)/2
                rtotfull(cnt) = sum(same*total)
                irec = irec + same*cnt
                print*, cnt
            endif
        enddo

        allocate(iuni(cnt))
        allocate(rmat(cnt,size(matrix(1,:))))
        allocate(rtot(cnt))
        
        iuni = iunifull(1:cnt)
        rmat = matrix(iuni,:)
        rtot = rtotfull(1:cnt)
        print*, 'reduced finished'
    END SUBROUTINE


END MODULE uniquemodule



MODULE membermodule

    implicit none 

    contains


    SUBROUTINE ismember(col,matrix,memval,colnum)

        use types
    
        integer(kind=ikind), intent(in), dimension(4)        :: col
        integer(kind=ikind), intent(in), dimension(:,:)      :: matrix
        integer(kind=ikind), dimension(size(matrix(:,1)))    :: subcol
        logical, intent(out)                                 :: memval
        integer(kind=ikind), intent(out)                     :: colnum

        subcol = sum(abs(matrix - spread(col, 1, size(matrix(:,1)))), dim=2)
        colnum = minloc(subcol, dim=1)

        if (subcol(colnum) == 0) then
           memval = .true.
        else
           colnum = -1
           memval = .false.
        endif

    END SUBROUTINE


END MODULE membermodule



MODULE twordmreader

    implicit none

    contains


    SUBROUTINE reduce_density(matst,totalst,m1,m2,m3,m4,total2)

        use types
        use membermodule
        use uniquemodule

        real(kind=dp), intent(in), dimension(:)                    :: totalst
        integer(kind=ikind), intent(in), dimension(:,:) :: matst
        integer(kind=ikind), dimension(:), intent(out), allocatable      :: m1, m2, m3, m4
        real(kind=dp), dimension(:), intent(out),allocatable             :: total2
        integer(kind=ikind)                                              :: i, j, k, l, cc
        integer(kind=ikind)                                              :: n, sdr
        integer(kind=ikind)                                              :: cnt, num
        logical                                                          :: memb

        integer(kind=ikind), dimension(:,:), allocatable                 :: newmat
        real(kind=dp), dimension(:), allocatable      :: stotal,newtotal
        integer(kind=ikind), dimension(:,:), allocatable    :: smat
        integer(kind=ikind)                                              :: mo1, mo2, mo3, mo4
        integer(kind=ikind), dimension(4)                                :: b

        integer(kind=ikind), dimension(:,:), allocatable                 :: mat2

        integer(kind=ikind), dimension(:), allocatable :: num1


!        allocate(newtotal(size(matst(:,1))),newmat(size(matst(:,1)), size(matst(1,:))))
!        newmat=0
!        newtotal=0.0_dp
!        n=maxval(matst)
!        print*, n**4, size(matst(:,1))
!        cnt = 0
!        do i = 1, n
!            do j = 1, n
!                do k = 1, n
!                    do l = 1, n
!                        if (l==1  .and. k==9 .and. j==9 .and. i==8) then
!                            print*, 'Thats the point we want', cnt, newmat(1,:)
!                        end if
!                        call ismember((/ i, j, k, l /), matst, memb, num)
!                        if (memb) then
!                            print*, 'memb var',cnt, memb
!                            cnt = cnt + 1
!                            newtotal(cnt) = totalst(num)
!                            newmat(cnt,:) = (/i, j, k, l /)
!                            print*,cnt,newmat(1,:), newtotal(1)
!                        else
!                            call ismember((/ j, i, l, k /), matst, memb, num)
!                            if (memb) then
!                                cnt = cnt + 1
!                                newtotal(cnt) = totalst(num)
!                                newmat(cnt,:) = (/ j, i, l, k /)
!                                print*,'holaaaa',cnt
!                            else
!                                call ismember((/ l, k, j, i /), matst, memb, num)
!                                if (memb) then
!                                    cnt = cnt + 1
!                                    newtotal(cnt) = totalst(num)
!                                    newmat(cnt,:) = (/ l, k, j, i /)
!                                    print*,'HOLAAAAA',cnt
!                                endif
!                            endif
!
!                        endif
!                        if (l==1  .and. k==9 .and. j==9 .and. i==8) then
!                            print*, 'Thats the point we wanted', cnt, newmat(1,:)
!                        end if
!                    enddo
!                enddo
!            enddo
!        enddo
!        print*,cnt
!        allocate(stotal(cnt), smat(cnt,4))
!        stotal = newtotal(1:cnt)
!        smat = newmat(1:cnt,:)
!
!!       in the following a few things have to be done; for example, the unique function
!!       is called (but we have that only in Python so far?)
!
!        print*,'size',size(smat(:,1))
!        print*,newmat(1,:),smat(1,:)

        call unique_total(matst, totalst, mat2, total2 )
        
        sdr = size(mat2(:,1))

        allocate(m1(sdr), m2(sdr), m3(sdr), m4(sdr))

        !uniquetotal(smat,mat2,stotal,total2)


        m1 = mat2(:,1)
        m2 = mat2(:,2)
        m3 = mat2(:,3)
        m4 = mat2(:,4)
        
        !newmat = (/ transpose(m1), transpose(m2), transpose(m3), transpose(m4) /)

   
    END SUBROUTINE reduce_density

    subroutine maxcoincidence(confs, ep2,ndiff)
        use types
        implicit none
        integer(kind=ikind), intent(in), dimension(:,:) :: confs
        integer(kind=ikind), intent(out), dimension(:,:), allocatable :: ep2, ndiff
        integer(kind=ikind), dimension(size(confs(:,1)), size(confs(1,:))):: matdum
        integer(kind=ikind), dimension(:,:), allocatable :: mat1
        integer(kind=ikind) :: i,j, c1,c2, count


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
            print*, 'mat1 constructed'
            do c1=1, size(confs(:,1))
                do c2=c1+1,size(confs(:,1))
                    do i=1,size(mat1(1,:))
                        if  (mat1(c1,i) /= mat1(c2,i)) then


                            do j=1,size(mat1(1,:))
                                if (mat1(c1,i) /= mat1(c2,j)) then
                                    ep2(c1,c2)=-ep2(c1,c2)
                                    ep2(c2,c1)=ep2(c1,c2)

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

    print*,ndiff(1,2)
    end subroutine maxcoincidence

    subroutine createtwordm(confs,civs,ndiff,ep2,mat,total)

    implicit none

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

    integer(kind=ikind), intent(in),dimension(:,:) :: confs
    integer(kind=ikind), intent(in), dimension(:,:) :: ep2, ndiff
    real(kind=dp), intent(in), dimension(:,:) :: civs
    real(kind=dp), intent(out), dimension(:), allocatable :: total
    integer(kind=ikind), intent(out), dimension(:,:), allocatable :: mat

    real(kind=dp), dimension(:,:,:,:), allocatable :: twordm
    integer(kind=ikind) :: ep,nc1,lconfs,norbs,nc2,sorb,rorb,qorb,porb,p,q,r,s,c1,c2,count1,count

    integer(kind=ikind) :: i,i1,i2,n,count2,eg, ndiff1

    integer(kind=ikind), dimension(:), allocatable :: mat1,mat2
    logical(4) :: sdef, rdef, pdef, qdef
    logical(4), dimension(:), allocatable :: logicaltwordms
    real(kind=dp) :: cutoff
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
    ep=1

    do c1=1,size(confs(:,1))
        do c2=1,size(confs(:,1))

            ndiff1=ndiff(c1,c2)
            ep=ep2(c1,c2)

            if (ndiff(c1,c2)/=0 .and. ndiff1 <=4) then

                sdef = .False.
                rdef = .False.
                pdef = .False.
                qdef = .False.
                if (allocated(diffs1)) deallocate(diffs1)
                if (allocated(diffs2)) deallocate(diffs2)
                if (allocated(spin1)) deallocate(spin1)
                if (allocated(spin2)) deallocate(spin2)
                allocate(diffs1(ndiff(c1,c2)/2), diffs2(ndiff(c1,c2)/2), spin1(ndiff(c1,c2)/2))
                allocate(spin2(ndiff(c1,c2)/2))

                count1=1
                count2=1

                do n=1,size(confs(c1,:))
                    if (confs(c1,n) /= confs(c2,n)) then
                        if (confs(c1,n) /= 0) THEN
                            diffs1(count1)=nint((n) / 2.0_dp + 0.1)
                            spin1(count1)=confs(c1,n)
                            count1=count1+1


                        elseif (confs(c2,n) /= 0) THEN
                            diffs2(count2)=nint((n) / 2.0_dp + 0.1)
                            spin2(count2)=confs(c2,n)
                            count2=count2+1
                        end if
                    end if
                enddo

                if (ndiff(c1,c2) == 4) then

                    sorb = diffs2(1)
                    qorb = diffs2(2)
                    porb = diffs1(2)
                    rorb = diffs1(1)

                    eg = 1.00

                    if (spin2(1) == spin1(1) .and. spin2(2) == spin1(2)) THEN
                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
                                civs(c1,1) * civs(c2,1) * ep * eg

                    end if

                    sorb = diffs2(1)
                    qorb = diffs2(2)
                    rorb = diffs1(2)
                    porb = diffs1(1)
                    eg = -1.00

                    if (spin2(1) == spin1(2) .and. spin2(2) == spin1(1)) then
                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ civs(c1,1) * civs(c2,1) * ep * eg
                    endif
                    qorb = diffs2(1)
                    sorb = diffs2(2)
                    porb = diffs1(2)
                    rorb = diffs1(1)
                    eg = -1.00

                    if (spin1(1) == spin2(2) .and. spin2(1) == spin1(2) ) THEN
                        twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb )+ &
                                civs(c1,1) * civs(c2,1) * ep * eg


                    end if

                    qorb = diffs2(1)
                    sorb = diffs2(2)
                    rorb = diffs1(2)
                    porb = diffs1(1)
                    eg = 1.00

                    if (spin1(2) == spin2(2) .and. spin2(1) == spin1(1)) then
                        twordm(porb , rorb , sorb , qorb )= twordm(porb , rorb , sorb , qorb ) &
                                + civs(c1,1) * civs(c2,1) * ep * eg

                    end if


                elseif (ndiff(c1,c2) == 2) THEN

                    qorb = diffs2(1)
                    porb = diffs1(1)
                    eg = 1.00

                    do i=1,size(confs(c2,:))

                        if (confs(c2,i) /=0) then
                            sdef = .True.

                            sorb =nint(i / 2.0_dp + 0.1)

                            spins = confs(c2,i)
                        end if

                        if (confs(c1,i) /= 0) then

                            rdef = .True.
                            rorb = nint((i) / 2.0_dp + 0.1)
                            spinr = confs(c1,i)
                        end if

                        if (sdef .and. rdef .and. spins == spinr .and. spin2(1) == spin1(1)) then
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.

                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) &
                                    + civs(c1,1) * civs(c2,1) * ep * eg



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
                            spins = confs(c2,i)
                        endif
                        if (confs(c1,i) /= 0) then
                            pdef = .True.
                            porb = nint((i) / 2.0_dp + 0.1)
                            spinp = confs(c1,i)
                        end if

                        if (sdef .and. pdef .and. spin1(1) == spins .and. spin2(1) == spinp) then
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.
                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) + &
                                    civs(c1,1) * civs(c2,1) * ep * eg


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
                            spinq = confs(c2,i)
                        endif
                        if (confs(c1,i) /= 0) then
                            rdef = .True.
                            rorb = nint((i) / 2.0_dp + 0.1)
                            spinr = confs(c1,i)
                        endif
                        if (rdef .and. qdef .and. spins == spinr .and. spin1(1) == spin2(1)) then
                            sdef = .False.
                            rdef = .False.
                            pdef = .False.
                            qdef = .False.
                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) +&
                                    civs(c1,1) * civs(c2,1) * ep * eg


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
                            spinq = confs(c2,i)

                        end if

                        if (confs(c1,i) /= 0) then
                            pdef = .True.
                            porb = nint((i) / 2.0_dp + 0.1)
                            spinp = confs(c1,i)

                        end if

                        if (qdef .and. pdef .and. spinq == spinp .and. spin2(1) == spin1(1)) then
                            qdef = .False.
                            pdef = .False.

                            twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) &
                                    + civs(c1,1) * civs(c2,1) * ep * eg

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

                            if (confs(c1,i1) /= 0 .and. confs(c1,i2) /= 0) then

                                sorb = nint(i1 / 2.0_dp + 0.1)
                                qorb = nint(i2 / 2.0_dp + 0.1)

                                if (confs(c1,i1) == confs(c1,i2)) then
                                    porb = sorb
                                    rorb = qorb

                                    eg = -1.00

                                    twordm(porb , rorb , sorb , qorb ) = twordm(porb , rorb , sorb , qorb ) &
                                            + civs(c1,1) * civs(c2,1) * ep * eg

                                end if
                                porb = qorb
                                rorb = sorb
                                eg = 1.00

                                twordm(porb , rorb , sorb , qorb )  = twordm(porb , rorb , sorb , qorb ) &
                                        + civs(c1,1) * civs(c2,1) * ep * eg
                            endif
                        end if

                    enddo
                end do
            end if


        end do

    end do



    cutoff = 1E-09
    count=0
    count2=1

    allocate(logicaltwordms(norbs**4))
    allocate(totaldum(norbs**4), matdum(norbs**4,4))
    logicaltwordms(:)=.False.
    do p=1,norbs
        do q=1,norbs
            do r=1,norbs
                do s=1,norbs
                    totaldum(count2)=twordm(p,q,r,s)
                    matdum(count2,:)=(/p,s,q,r/)
                    if (abs(twordm(p,q,r,s))>=cutoff) then
                        count=count+1
                        logicaltwordms(count2)=.True.

                    end if
                    count2=count2+1
                end do
            end do
        end do
    end do


    allocate(mat(count, 4), total(count))
    count=1
    do i=1,count2
        if (logicaltwordms(i)) then
            mat(count,:)=matdum(i,:)
            total(count)=totaldum(i)
            count=count+1
        end if
    end do

        print*, 'twordm calculated'



    end subroutine createtwordm

END MODULE twordmreader

