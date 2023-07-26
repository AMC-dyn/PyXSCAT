module linspace
              use mod_types
    implicit none
    contains
    subroutine linspace_1(from, to, n, array)
    integer(kind=ikind), intent(in) :: n
    real(kind = dp), intent(in) :: from, to
    real(kind = dp), intent(out), dimension(:), allocatable :: array
    real(kind = dp) :: range
    integer(kind=ikind) :: i

    allocate(array(n))


    range = to - from

    do i = 1, n
        array(i) = from + range * (i - 1) / (n - 1)
    end do
    print*,array
    return
end subroutine
end module linspace
