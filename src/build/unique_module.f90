MODULE uniquemodule

    implicit none

    contains


    SUBROUTINE unique_integer(matrix,rmat,iuni,irec)


          INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
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


          INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
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


          INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
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
