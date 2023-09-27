module linspace
    implicit none
    contains
    subroutine linspace_1(from, to, n, array)
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
    integer(kind=ikind), intent(in) :: n
    real(kind = dp), intent(in) :: from, to
    real(kind = dp), intent(out), dimension(:), allocatable :: array
    real(kind = dp) :: range
    integer :: i

    allocate(array(n))


    range = to - from

    do i = 1, n
        array(i) = from + range * (i - 1) / (n - 1)
    end do

    return
end subroutine
end module linspace