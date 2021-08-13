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
        
        do i = 1, n
            do j = 1, n
                do k = 1, n 
                    do l = 1, n
                        call ismember((/ i, j, k, l /), matst, memb, num)
                        if (memb) then 
                            newtotal = [newtotal; totalst(num)]
                            newmat = [newmat; [i j k l]]
                        else
                            call ismember((/ j, i, l, k /), matst, memb, num)
                            if (memb) then
                                newtotal = [newtotal; totalst(num)]
                                newmat = [newmat; [j i l k]]
                            else
                                call ismember((/ l, k, j, i /), matst, memb, num)
                                if (memb) then
                                    newtotal = [newtotal; totalst(num)]
                                    newmat = [newmat; [l k j i]]
                                endif
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo

        newmat = int8(newmat)
        mat2 = unique(newmat,'rows','stable')
        [sdR,~] = size(mat2)
        
        do cc = 1, sdr
            mo = mat2(cc,1)
            mo2 = mat2(cc,2)
            mo3 = mat2(cc,3)
            mo4 = mat2(cc,4)
            B = [mo mo2 mo3 mo4]
            call ismember(newmat, mat2(cc,:), tf, num)
            m1(cc) = mo
            m2(cc) = mo2
            m3(cc) = mo3
            m4(cc) = mo4
            total2(cc) = sum(newtotal(tf))
        enddo
        
        newmat = [m1' m2' m3' m4']
        newtotal = total2'
   
    END SUBROUTINE



END MODULE

