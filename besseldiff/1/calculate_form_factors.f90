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


    function sinc (a)
             INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
            INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
             real(kind=dp):: sinc, a
             if (abs(a) < 1.0d-10) then
                sinc = 1
             else
                sinc = sin(a) / (a)
             end if
    end function
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
        character(len=22), intent(in) :: filename
        character(len=10) :: dumm
        character(len=10),intent(out), dimension(55):: atoms
        real(kind=dp), intent(out),dimension(55,13) :: ff
        real(kind=dp), intent(out),dimension(13):: q
        q(1)=0.0_dp
        ff(:,1)=0.0_dp
        open(unit=15,file= filename)
        read(15,*)dumm,q(2:13)
        do i=1,55
            read(15,*)atoms(i),ff(i,2:13)
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

subroutine sphbes(n,z,j)
    implicit none
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
    integer(kind=ikind), intent(in) :: n
    real(kind=dp), intent(in) :: z
    real(kind=dp), intent(out) :: j

    if (z==0) then
        if (n==0) then
            j=1.00_dp
        else
            j=0.0_dp
        end if

    else

    j=dsqrt(dacos(-1.0_dp)/(2.0_dp*z))*BESSEL_JN(n, z)
    endif


    end subroutine sphbes

function ttra(b1,b2,q1,q2,anglebond,aq) result(res)
    implicit none
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
    real(dp), intent(IN)  :: q1,q2,anglebond,b1,b2,aq
    real(dp) :: res
    real(dp) :: cutoff,dres
    real(dp) :: sp1,sp2
    integer(ikind) :: l, lmax, lmin

    cutoff=0.00000001
    lmin=12
    lmax=100
    l=2
    sp1=spherical_bessel_jn(0, q1*b1)
    sp2=spherical_bessel_jn(0, q2*b2)

    dres=sp1*sp2


    res=dres
    do while (cutoff<=dres/res .and. l<lmax )
        if (l<lmin) then
            cutoff=0.0_dp
        else
           cutoff=0.00000001
        end if
        sp1=spherical_bessel_jn(l, q1*b1)
        sp2=spherical_bessel_jn(l, q2*b2)

        dres=(2*l+1)*sp1*sp2*pl(anglebond,l)*pl(cos(aq),l)

        res=res+dres

        l=l+2
    end do
    if (isnan(res)) then
            PRINT*,'res',sp1,sp2,b1,b2,q1,q2,l
    end if

    end function ttra

    function pl(x,n)
!======================================
! calculates Legendre polynomials Pn(x)
! using the recurrence relation
! if n > 100 the function retuns 0.0
!======================================
double precision pl
double precision x
double precision pln(0:n)
integer n, k


pln(0) = 1.0
pln(1) = x

if (n <= 1) then
  pl = pln(n)
  else
  do k=1,n-1
    pln(k+1) = ((2.0*k+1.0)*x*pln(k) - float(k)*pln(k-1))/(float(k+1))
  end do
  pl = pln(n)
end if
return
end



real(kind=SELECTED_REAL_KIND(15)) function spherical_bessel_jn(n, x) result(r)
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
integer, intent(in) :: n
real(dp), intent(in) :: x
integer :: nm
real(dp) :: sj(0:n), dj(0:n)
call sphj(n, x, nm, sj, dj)
if (nm /= n) then
    print*,'ERROR'
end if
r = sj(n)
end function
 INTEGER FUNCTION MSTA1(X,MP)
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that the magnitude of
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point
!       ===================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

        IMPLICIT integer (I-N)
        A0=DABS(X)
        N0=INT(1.1D0*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20
           NN=int(N1-(N1-N0)/(1.0D0-F0/F1))
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
           F1=F
 10     END DO
 20     MSTA1=NN
        RETURN
        END

        INTEGER FUNCTION MSTA2(X,N,MP)
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT integer (I-N)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1D0*A0)+1
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=int(N1-(N1-N0)/(1.0D0-F0/F1))
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
           F1=F
10      END DO
20      MSTA2=NN+10
        RETURN
        END

        real(kind=SELECTED_REAL_KIND(15)) function envj(n, x) result(r)
            INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        r = log10(6.28_dp*n)/2 - n*log10(1.36_dp*x/n)
        end function



 SUBROUTINE SPHJ(N,X,NM,SJ,DJ)
!       =======================================================
!       Purpose: Compute spherical Bessel functions jn(x) and
!                their derivatives
!       Input :  x --- Argument of jn(x)
!                n --- Order of jn(x)  ( n = 0,1,… )
!       Output:  SJ(n) --- jn(x)
!                DJ(n) --- jn'(x)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       =======================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

        IMPLICIT integer (I-N)
        DIMENSION SJ(0:N),DJ(0:N)
        NM=N
        IF (DABS(X).LT.1.0D-100) THEN
           DO 10 K=0,N
              SJ(K)=0.0D0
              DJ(K)=0.0D0
10         END DO
           SJ(0)=1.0D0
           IF (N.GT.0) THEN
              DJ(1)=.3333333333333333D0
           ENDIF
           RETURN
        ENDIF
        SJ(0)=DSIN(X)/X
        DJ(0)=(DCOS(X)-DSIN(X)/X)/X
        IF (N.LT.1) THEN
           RETURN
        ENDIF
        SJ(1)=(SJ(0)-DCOS(X))/X
        IF (N.GE.2) THEN
           SA=SJ(0)
           SB=SJ(1)
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F=0.0D0
           F0=0.0D0
           F1=1.0D0-100
           DO 15 K=M,0,-1
              F=(2.0D0*K+3.0D0)*F1/X-F0
              IF (K.LE.NM) SJ(K)=F
              F0=F1
              F1=F
15         END DO
           CS=0.0D0
           IF (DABS(SA).GT.DABS(SB)) CS=SA/F
           IF (DABS(SA).LE.DABS(SB)) CS=SB/F0
           DO 20 K=0,NM
              SJ(K)=CS*SJ(K)
20         END DO
        ENDIF
        DO 25 K=1,NM
           DJ(K)=SJ(K-1)-(K+1.0D0)*SJ(K)/X
25      END DO
        RETURN
        END




SUBROUTINE polint(xa,ya,n,x,y,dy)
      ! Polynomial interpolation. Copied form "Numerical Recipes in
      ! FORTRAN 77"
      !
      ! Bo Terp Paulsen, botp@mek.dtu.dk

      INTEGER n,NMAX
      real(kind=SELECTED_REAL_KIND(15)) :: dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10) ! Largest anticipated value of n.
      !Given arrays xa and ya, each of length n, and given a value x, this routine
      !returns a value y, and an error estimate dy. If P (x) is the polynomial of
      !degree N − 1 such that P (xai) = yai , i = 1, . . . , n, then the returned
      !value y = P (x ).
      INTEGER i,m,ns
      real(kind=SELECTED_REAL_KIND(15)):: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

      ns=1
      dif = abs(x-xa(1))
      DO i=1,n  !Here we find the index ns of the closest table entry,
          dift=abs(x-xa(i))

          IF (dift<dif) THEN
              ns=i
              dif=dift
          ENDIF

          c(i)=ya(i) !and initialize the tableau of c’s and d’s.
          d(i)=ya(i)
      ENDDO

      y=ya(ns) !This is the initial approximation to y.
      ns=ns-1

      DO m=1,n-1 !for each column of the tableau, we loop over the current c’s and d’s and update them.
          DO i=1,n-m
              ho=xa(i)-x
              hp=xa(i+m)-x
              w=c(i+1)-d(i)
              den=ho-hp
              IF(den==0.)then
                 print*, 'failure in polint'
                 stop
              end if
              !This error can occur only if two input xa’s are (to within roundoff) identical.
              den=w/den
              d(i)=hp*den ! Here the c’s and d’s are updated.
              c(i)=ho*den
          ENDDO
          IF (2*ns<n-m)THEN
          !After each column in the tableau is completed, we decide
          !which correction, c or d, we want to add to our accu-
          !mulating value of y, i.e., which path to take through
          !the tableau—forking up or down. We do this in such a
          !way as to take the most “straight line” route through the
          !tableau to its apex, updating ns accordingly to keep track
          !of where we are. This route keeps the partial approxima-
          !tions centered (insofar as possible) on the target x. The
          !last dy added is thus the error indication.
              dy=c(ns+1)
          ELSE
              dy=d(ns)
              ns=ns-1
          ENDIF
      y=y+dy
      ENDDO
      RETURN
      END


    !
!  Function f(x)
!


   subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
implicit none
integer n
double precision x(n), y(n), b(n), c(n), d(n)
integer i, j, gap
double precision h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 3: forward elimination
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine spline

  function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
double precision ispline
integer n
double precision  u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline
End Module calculate_form_factors
