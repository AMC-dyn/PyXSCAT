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
