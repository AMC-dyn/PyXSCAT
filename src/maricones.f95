!-----------------------------------------------------------------------
! 2RDM Fortran Modules for PyXSCAT code
! Andres Moreno Carrascosa and Mats Simmermacher, 2021
!-----------------------------------------------------------------------

module types
    implicit none
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
end module types

module twordmshit

    implicit none 

    contains



    SUBROUTINE ismember(col,matrix)
        use types

        integer(kind=ikind), intent(in), dimension(4)      :: col
        integer(kind=ikind), intent(in), dimension(:,4)    :: matrix
        integer(kind=ikind), dimension(size(matrix(:,1)))  :: subcol
        logical, intent(out)                               :: memval
        integer(kind=ikind), intent(out)                   :: colnum

        subcol = sum(abs(matrix - spread(col, 1, size(matrix(:,1)))), dim=1)
        colnum = findloc(subcol, 0)
        if (size(colnum) > 0) then
            memval = .true.
        else
            memval = .false.
        endif            

    END SUBROUTINE



    SUBROUTINE whateverthename(whateverthevariables)
        use types

        real(kind=dp), dimension(:), allocatable :: h_saved, pmu, h_sum, h_0, h_1, h_r, muoh,zij, zij2, zkr, zkr2
       
        allocate(zij(size(Z(:,1,1))), zij2(size(Z(:,1,1))), zkr(size(Z(:,1,1))), zkr2(size(Z(:,1,1))))
        
        newtotal = []
        newmat = []
        
        do i = 1, n
            do j = 1, n
                do k = 1, n 
                    do l = 1, n
                        if ~ismember([i,j,k,l],matst,'rows') then
                            if ismember([j i l k],matst,'rows') then
                                [~,num] = ismember([j i l k],matst,'rows')
                                newtotal = [newtotal;totalst(num)]
                                newmat = [newmat;[j i l k]]
                            elseif ismember([l,k,j,i],matst,'rows') then
                                [~,num] = ismember([l,k,j,i],matst,'rows')
                                newtotal = [newtotal;totalst(num)]
                                newmat = [newmat;[l k j i]]
                            endif
                        else
                            [~,num] = ismember([i,j,k,l],matst,'rows')
                            newtotal = [newtotal;totalst(num)]
                            newmat = [newmat;[i j k l]]
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
            [tf, ~]=ismember(newmat,mat2(cc,:),'rows')
            m1(cc) = mo
            m2(cc) = mo2
            m3(cc) = mo3
            m4(cc) = mo4
            total2(cc) = sum(newtotal(tf))
        enddo
        
        newmat = [m1' m2' m3' m4']
        newtotal = total2'
