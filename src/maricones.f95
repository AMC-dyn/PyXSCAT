!-----------------------------------------------------------------------
! 2RDM Fortran Modules for PyXSCAT code
! Andres Moreno Carrascosa and Mats Simmermacher, 2021
!-----------------------------------------------------------------------

MODULE types
    implicit none
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
END MODULE types



MODULE twordmshit

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



    SUBROUTINE whateverthename(whateverthevariables)

        use types

        newtotal = []
        newmat = []
        
        integer(kind=ikind)                                 :: i, j, k, l, cc
        integer(kind=ikind)                                 :: n, sdr
        integer(kind=ikind)                                 :: cnt, num
        logical                                             :: memb
        integer(kind=ikind), dimension(n**4)                :: newtotal
        integer(kind=ikind), dimension(N**4,4)              :: newmat
        integer(kind=ikind), dimension(:), allocatable      :: stotal
        integer(kind=ikind), dimension(:,:), allocatable    :: smat   
        integer(kind=ikind)                                 :: mo1, mo2, mo3, mo4
        integer(kind=ikind), dimension(4)                   :: b
        integer(kind=ikind), dimension(sdr)                 :: m1, m2, m3, m4, total2
        integer(kind=ikind), dimension(:,:), allocatable    :: mat2

        cnt = 0
        do i = 1, n
            do j = 1, n
                do k = 1, n 
                    do l = 1, n
                        call ismember((/ i, j, k, l /), matst, memb, num)
                        if (memb) then
                            cnt = cnt + 1 
                            newtotal(cnt) = totalst(num)
                            newmat(cnt) = (/i, j, k, l /)
                        else
                            call ismember((/ j, i, l, k /), matst, memb, num)
                            if (memb) then
                                cnt = cnt + 1
                                newtotal(cnt) = totalst(num)
                                newmat(cnt) = (/ j, i, l, k /)
                            else
                                call ismember((/ l, k, j, i /), matst, memb, num)
                                if (memb) then
                                    cnt = cnt + 1
                                    newtotal(cnt) = totalst(num)
                                    newmat(cnt) = (/ l, k, j, i /)
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

        smat = int8(smat)
        mat2 = unique(smat,'rows','stable')
        
        sdr = size(mat2, dim=1)
        
        do cc = 1, sdr
            mo1 = mat2(cc,1)
            mo2 = mat2(cc,2)
            mo3 = mat2(cc,3)
            mo4 = mat2(cc,4)
            b = (/ mo, mo2, mo3, mo4 /)
            call ismember(smat, mat2(cc,:), memb, num)
            m1(cc) = mo1
            m2(cc) = mo2
            m3(cc) = mo3
            m4(cc) = mo4
            total2(cc) = sum(stotal(num))
        enddo
        
        newmat = (/ transpose(m1), transpose(m2), transpose(m3), transpose(m4) /)  ! ???
        newtotal = transpose(total2)  ! ???
   
    END SUBROUTINE



END MODULE

