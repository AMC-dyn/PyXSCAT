
program tryingbess
    implicit none
    integer, parameter:: nq=10
    real(kind=8) ::x(nq)
    integer:: i
    real(kind=8)::rbes2(1001,nq),rbes(116,nq),rneu(1001,nq), time1,time2
    integer:: mssgae



    x=10.d0


 call  dphrec(x,rbes2,1001,2.D-8, 0, 1001,nq)
call cpu_time(time1)



do i=1,1001
call  dphrec(x,rbes,i,2.D-8, 0, i,nq)
    print*,i, sum(abs(rbes(1:16,1)-rbes2(1:16,1)))
enddo

call cpu_time(time2)


print*,time2-time1



end



subroutine dphrec(x,rbes,leps1,epsl,normal,lmax1,nq)

IMPLICIT REAL(kind=8)(A-H,o-Z)
DIMENSION RBES(LEPS1,nq), RNEU(LEPS1,nq)
DIMENSION X(nq), XX(nq), CX(nq), DX(nq), CU(nq), A(nq),SX(nq)

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
RBES(LU-1,:)=RBES(LU,:)-XX/(4.D0*LU**2-1.D0)*RBES(LU+1,:)
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
END

  subroutine van (jbwd,n, x)
        implicit none
        double precision, intent(in):: x
        integer, intent(in) :: n
        integer i, nmax
        double precision, intent(out),dimension(0:N):: jbwd
        double precision ::  ix, jmx1, jmx, j0, j1


        nmax= max(n, ceiling(1.142476370122814*n-4.027048776987268))


        ix = 1 / x
        j0 = sin(x) * ix
        j1 = (j0 - cos(x)) * ix
        !
        ! Backward recursion.
        !
        jmx1 = 0
        jmx = 1
        do i = n + nmax, n + 2, -1
            jbwd(0) = (2 * i + 1) * (ix * jmx - jmx1 / (2 * i + 1))
            jmx1 = jmx
            jmx = jbwd(0)
        end do
        do i = n + 1, 1, -1
            jbwd(i - 1) = (2 * i + 1) * ix * jmx - jmx1
            jmx1 = jmx
            jmx = jbwd(i - 1)
        end do
        if (abs(jmx) >= abs(jmx1)) then
            j0 = j0 / jbwd(0)
        else
            j0 = j1 / jbwd(1)
        end if
        jbwd = j0 * jbwd
    end

