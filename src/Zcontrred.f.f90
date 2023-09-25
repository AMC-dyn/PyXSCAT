program Zcontrred

    use linspace
    implicit none
     INTEGER, PARAMETER :: dpp = SELECTED_REAL_KIND(kind(1.d0))
     INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
     INTEGER, PARAMETER :: ikind = SELECTED_int_KIND(8)
     integer,parameter :: n=100
     real(kind=SELECTED_REAL_KIND(15)):: fact
     real(kind=dp),allocatable,dimension(:)::q
    real(kind=dp), dimension(n):: qe
    real(kind=dp),dimension(500,n)::result
    real(kind=dp),dimension(500,n)::result2
     real(kind=dp),dimension(n)::x2
     real(kind=dp),dimension(500,n)::factvec
      real(kind=dp):: time1, time2, time3
     integer(kind=ikind)::i,j,k,r,nq



    call linspace_1(0.00001_dp,10.0_dp,n,q)
    print*,'after linspace'

    x2(:)=100.0

    qe=q
    i=2
call cpu_time(time1)

    print*,'before dphrec'
  do j=1,100000
    call dphrec(qe*4.d0,result,500,500,n)
enddo
call cpu_time(time2)

    print*,result(6,3)

    call dphrec(qe,result2,500,500,n)
    x2=0.0_dp
     do i=0,200
         factvec(i,:)=(-1.0_dp)**float(i)/fact(float(i))*(4.0**2.0_dp-1.0_dp)**float(i)*(qe/2.d0)**float(i)
    end do
call cpu_time(time2)
    print*,factvec(:,1)
do j=1,100000
     x2=0.0_dp
    do i=0,200

     x2=x2+(4.0**2.0_dp-1.0_dp)*factvec(i,:)*result2(i+6,:)

    end do
    x2=x2*4.0_dp**5.0_dp
enddo
    print*,abs(x2(3)-result(6,3))
call cpu_time(time3)
    print*,time2-time1
    print*,time3-time2

end program Zcontrred


subroutine dphrec(x,rbes,leps1,lmax1,nq)
 INTEGER, PARAMETER :: dpp = SELECTED_REAL_KIND(kind(1.d0))
 INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
 INTEGER, PARAMETER :: ikind = SELECTED_int_KIND(8)
integer(kind=ikind):: K,LU,LEPS1,nq,lmax1
real(kind=dp),DIMENSION(LEPS1,nq) :: RNEU
 real(kind=dp),DIMENSION(LEPS1,nq), intent(out)::RBES
real(kind=dp),DIMENSION(nq):: X, XX, CX, DX, CU, A,SX


XX=X*X
!LEPS=0.5D0*dsqrt(XX/EPSL**(1.D0/3.00)+9.00) + 0.5D0
!IF(LEPS>=LEPS1) GOTO 101

rbes=0.0d0
rneu=0.0d0
RBES(LEPS1,:) =1.0D0
RBES(LEPS1-1,:)=1.0D0
CX=DCOS(X)
SX=DSIN(X)
RNEU(1,:)=CX
RNEU(2,:)=CX+X*SX

DO K=3,LEPS1
LU=LEPS1-K+2
RBES(LU-1,:)=RBES(LU,:)-XX/(4.D0*LU*LU-1.D0)*RBES(LU+1,:)
end do
A=RBES(1,:)*RNEU(2,:)-XX/3.0d0*RBES(2,:)*CX
DO K=1, LEPS1
RBES(K,:)=RBES(K,:)/A
end do
CU= 1.0D0/X

DO K=1, LMAX1
W=2.D0*(K-1)
CU=X/(W+1.D0)*CU
RBES(K,:)= CU*RBES(K,:)
end do


message=0
    return
END

real(kind=SELECTED_REAL_KIND(15))  recursive function fact(n) result (ans)
implicit none
real(kind=SELECTED_REAL_KIND(15)), intent(in) :: n
    if(n == 0) then
        ans=1
        return
    end if

    if (n == 1) then
        ans = 1

    else

        ans = n * fact(n-1)
    end if

    return
end function fact