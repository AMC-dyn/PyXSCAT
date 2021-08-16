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



    SUBROUTINE unique(matrix,rmat,iuni,irec)

        integer, dimension(:,:), intent(in)                       :: matrix
        integer, dimension(:,:), allocatable                      :: sprmat
        integer, dimension(:), allocatable                        :: same, iunifull
        integer, dimension(:,:), intent(out), allocatable         :: rmat
        integer, dimension(:), intent(out), allocatable           :: iuni
        integer, dimension(:), intent(out), allocatable           :: irec
        integer                                                   :: i, cnt

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
                same = (-2)*same + 1
                same = (abs(same) + same)/2
                irec = irec + same*cnt
            endif
        enddo

        allocate(iuni(cnt))
        allocate(rmat(cnt,size(matrix(:,1))))

        iuni = iunifull(1:cnt)
        rmat = matrix(iuni,:)

    END SUBROUTINE



END MODULE

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

END MODULE



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
        integer(kind=ikind), dimension(:), allocatable      :: stotal,newtotal,ipos,irep
        integer(kind=ikind), dimension(:,:), allocatable    :: smat
        integer(kind=ikind)                                              :: mo1, mo2, mo3, mo4
        integer(kind=ikind), dimension(4)                                :: b

        integer(kind=ikind), dimension(:,:), allocatable                 :: mat2

        integer(kind=ikind), dimension(:), allocatable :: num1


        allocate(newtotal(size(matst(:,1))),newmat(size(matst(:,1)), size(matst(1,:))))
        cnt = 0
        do i = 1, n

            do j = 1, n
                do k = 1, n 
                    do l = 1, n
                        call ismember((/ i, j, k, l /), matst, memb, num)
                        if (memb) then
                            cnt = cnt + 1 
                            newtotal(cnt) = totalst(num)
                            newmat(cnt,:) = (/i, j, k, l /)
                        else
                            call ismember((/ j, i, l, k /), matst, memb, num)
                            if (memb) then
                                cnt = cnt + 1
                                newtotal(cnt) = totalst(num)
                                newmat(cnt,:) = (/ j, i, l, k /)
                            else
                                call ismember((/ l, k, j, i /), matst, memb, num)
                                if (memb) then
                                    cnt = cnt + 1
                                    newtotal(cnt) = totalst(num)
                                    newmat(cnt,:) = (/ l, k, j, i /)
                                endif
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo

        allocate(stotal(cnt), smat(cnt,4))
        stotal = newtotal(1:cnt)
        smat = newmat(1:cnt,:)

!       in the following a few things have to be done; for example, the unique function
!       is called (but we have that only in Python so far?)


       call unique(smat, mat2, ipos, irep)
        
        sdr = size(mat2, dim=1)

        allocate(m1(sdr), m2(sdr), m3(sdr), m4(sdr), total2(sdr))

        !uniquetotal(smat,mat2,stotal,total2)


        m1 = mat2(:,1)
        m2 = mat2(:,2)
        m3 = mat2(:,3)
        m4 = mat2(:,4)
        
        !newmat = (/ transpose(m1), transpose(m2), transpose(m3), transpose(m4) /)

   
    END SUBROUTINE reduce_density

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



    cutoff = 1E-10
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





    end subroutine createtwordm

END MODULE

