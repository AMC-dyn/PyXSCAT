Module calculate_form_factors
    implicit none
    contains


    subroutine linspace(from, to, n,array)
     INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
     INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
     integer,intent(in) :: n
    real(kind=dp), intent(in) :: from, to
    real(kind=dp), intent(out),dimension(n) :: array
    real(kind=dp) :: range

    integer :: i

    range = to - from

    do i=1, n
        array(i) = from + range * (i - 1) / (n - 1)
    end do
end subroutine

    subroutine obtain_form_factors(aff1,atom,q1,atoms,a,b,c)
        implicit none
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        real(kind=dp), parameter :: pi=dacos(-1.00_dp)
        character(len=10), intent(in) :: atom
        real(kind=dp), intent(in), dimension(:) :: q1
        real(kind=dp),intent(out),dimension(:),allocatable:: aff1
        character(len=10),intent(in),dimension(151):: atoms
        real(kind=dp), intent(in),dimension(4,151) :: a,b
        real(kind=dp), intent(in),dimension(151):: c
        integer(kind=ikind) :: index,i





        allocate(aff1(size(q1)))

        do i=1,151
           if (atom==atoms(i)) then
               index=i
            end if
        end do

        do i=1,size(q1)
            aff1(i)=c(index)+sum(a(:,index)*exp(-b(:,index)*(q1(i)/(4*pi))**2.0))
        enddo

    end subroutine obtain_form_factors

    subroutine table_of_ff(filename,atoms,a,b,c)
        implicit none
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer :: i
        character(len=16), intent(in) :: filename

        character(len=10),intent(out), dimension(151):: atoms
        real(kind=dp), intent(out),dimension(4,151) :: a,b
        real(kind=dp), intent(out),dimension(151):: c

        open(unit=15,file= filename)

        do i=1,151
            read(15,*)atoms(i),a(1,i) ,b(1,i),a(2,i),b(2,i),a(3,i),b(3,i),a(4,i),b(4,i),c(i)
            print*,i
        end do
        close(unit=15)
    end subroutine table_of_ff




End Module calculate_form_factors