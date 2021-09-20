Module calculate_form_factors
    use interp1D
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

    subroutine obtain_i_form_factors(atoms,atom,iff_or, q1,q_or,iff)
        use interp1D
        implicit none
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        real(kind=dp), parameter :: pi=dacos(-1.00_dp)
        character(len=10), intent(in) :: atom
        real(kind=dp), intent(in), dimension(12)::q_or
        real(kind=dp), intent(in), dimension(:) :: q1
        real(kind=dp), intent(in), dimension(:,:) :: iff_or
        real(kind=dp),intent(out),dimension(:):: iff
        character(len=10),intent(in),dimension(55):: atoms

        integer(kind=ikind) :: index,i







        do i=1,55
           if (atom==atoms(i)) then
               index=i
            end if

        end do
        print*,index,atom
        do i=1,12
            print*,q_or(i), iff_or(index,i)
        end do
        iff=interp_linear_vec(q_or,iff_or(index,:), q1)

    end subroutine obtain_i_form_factors

      function interp_linear_vec(x,y,xout) result(yout)
        ! Interpolate y from ordered x to ordered xout positions

        implicit none

        double precision, dimension(:), intent(IN) :: x, y
        double precision, dimension(:), intent(IN) :: xout
        double precision, dimension(size(xout)) :: yout
        integer :: i, j, n, nout

        n    = size(x)
        nout = size(xout)

!         write(*,*) minval(x), maxval(x), n, nout

        do i = 1, nout
            if (xout(i) < x(1)) then
                yout(i) = y(1)
!                 write(*,*) 1, xout(i)
            else if (xout(i) > x(n)) then
                yout(i) = y(n)
!                 write(*,*) 2, xout(i)
            else
                do j = 1, n
                    if (x(j) >= xout(i)) exit
                end do

                if (j == 1) then
                    yout(i) = y(1)
!                     write(*,*) 3, xout(i)
                else if (j == n+1) then
                    yout(i) = y(n)
!                     write(*,*) 4, xout(i)
                else
                    yout(i) = interp_linear_internal(x(j-1:j),y(j-1:j),xout(i))
!                     write(*,*) 5, xout(i)
                end if
            end if
        end do

        return

      end function interp_linear_vec

       function interp_linear_internal(x,y,xout) result(yout)

        implicit none
         INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        real(dp), intent(IN)  :: x(2), y(2), xout
        real(dp) :: yout
        real(dp) :: alph

        if ( xout .lt. x(1) .or. xout .gt. x(2) ) then
            write(*,*) "interp1: xout < x0 or xout > x1 !"
            write(*,*) "xout = ",xout
            write(*,*) "x0   = ",x(1)
            write(*,*) "x1   = ",x(2)
            stop
        end if

        alph = (xout - x(1)) / (x(2) - x(1))
        yout = y(1) + alph*(y(2) - y(1))

        return

    end function interp_linear_internal
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

        end do
        close(unit=15)
    end subroutine table_of_ff

    subroutine table_of_iff(filename,atoms,ff,q)
        implicit none
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer :: i
        character(len=40), intent(in) :: filename
        character(len=10) :: dumm
        character(len=10),intent(out), dimension(55):: atoms
        real(kind=dp), intent(out),dimension(55,12) :: ff
        real(kind=dp), intent(out),dimension(12):: q

        open(unit=15,file= filename)
        read(15,*)dumm,q(1:12)
        do i=1,55
            read(15,*)atoms(i),ff(i,:)
        end do
        close(unit=15)
    end subroutine table_of_iff


subroutine interp1( xData, yData, xVal, yVal )
! Inputs: xData = a vector of the x-values of the data to be interpolated
!         yData = a vector of the y-values of the data to be interpolated
!         xVal  = a vector of the x-values where interpolation should be performed
! Output: yVal  = a vector of the resulting interpolated values

  implicit none
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
  INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
  real(kind=dp), intent(in) :: xData(:), yData(:), xVal(:)
  real(kind=dp), intent(out) :: yVal(:)
  integer(kind=ikind) :: inputIndex, dataIndex
  real(kind=dp) :: minXdata, minYdata, xRange, weight,maxXDATA

  ! Possible checks on inputs could go here
  ! Things you may want to check:
  !   monotonically increasing xData
  !   size(xData) == size(yData)
  !   size(xVal) == size(yVal)

  minXData = xData(1)
  maxXData = xData(size(xData))
  xRange = maxXData - minXData

  do inputIndex = 1, size(xVal)
      ! possible checks for out of range xVal could go here

      ! this will work if x is uniformly spaced, otherwise increment
      ! dataIndex until xData(dataIndex+1)>xVal(inputIndex)
      dataIndex = floor((xVal(inputIndex)-minXData)/xRange);

      weight = (xVal(inputIndex) - xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex));
      yVal(inputIndex) = (1.0-weight)*yData(dataIndex) +weight*yData(dataIndex+1);
  end do
end subroutine



End Module calculate_form_factors