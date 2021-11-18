MODULE besseldiff1

    implicit none 

    contains



    SUBROUTINE rrdj1xy(lmax,x,a)
    ! This subroutine uses a recurrence relation to calculate the
    ! first set of expansion coefficients, a, at given coordinate
    ! x (or y) up to an angular momentum of L 

        integer, parameter    :: dp = selected_real_kind(15)
        integer, intent(in)                                        :: lmax
        integer                                                    :: l,p
        real(kind=dp), intent(in)                                  :: x
        real(kind=dp), dimension(:,:), intent(out), allocatable    :: a

        allocate( a(lmax+1,lmax+2) )
        
        a(1,2) = 1.0_dp
        a(2,3) = -x

        do l = 2, lmax
            a(l+1,2) = -(l-1) * a(l-1,2)
        enddo

        do l = 2, lmax
            do p = 1, l
                a(l+1,p+2) = -(l-1) * a(l-1,p+2) - x * a(l,p+1)
            enddo
        enddo

    END SUBROUTINE

 

    SUBROUTINE rrdj1z(nmax,z,a)
    ! This subroutine uses a recurrence relation to calculate the
    ! first set of expansion coefficients, a, at given coordinate
    ! z up to an angular momentum of L 

        integer, parameter    :: dp = selected_real_kind(15)
        integer, intent(in)                                        :: nmax
        integer                                                    :: n,s
        real(kind=dp), intent(in)                                  :: z
        real(kind=dp), dimension(:,:), allocatable                 :: a0
        real(kind=dp), dimension(:,:), intent(out), allocatable    :: a

        allocate( a0(nmax+2,nmax+3) )
        allocate( a(nmax+1,nmax+2) )
        
        a0(1,2) = -1.0_dp
        a0(2,3) = z

        do n = 2, nmax+1
            a0(n+1,2) = -(n-1) * a0(n-1,2)
        enddo

        do n = 2, nmax+1
            do s = 1, n
                a0(n+1,s+2) = -(n-1) * a0(n-1,s+2) - z * a0(n,s+1)
            enddo
        enddo

        do n = 1, nmax+1
            do s = 1, n+1
                a(n,s) = a0(n+1,s+1) 
            enddo
        enddo

    END SUBROUTINE

    
    
    SUBROUTINE dj1expand(l,m,n,q,x,y,z,ax,ay,az,sig)
    ! This subroutine expands the derivatives for given angular
    ! momenta, L, M, N, coordinates, Hx, Hy, Hz, and momentum transfer,
    ! q. It requires the sets of expansion coefficients, ax, ay, az,
    ! calculated at coordinates Hx, Hy, Hz, respectively,
    ! with the two routines above. 
    
        use calculate_form_factors
       
        integer, parameter    :: dp = selected_real_kind(15)
        integer, intent(in)                          :: l,m,n
        integer                                      :: p,r,s
        integer                                      :: beta
        integer                                      :: plen,rlen,slen
        real(kind=dp), intent(in)                    :: q,x,y,z
        real(kind=dp)                                :: rad
        real(kind=dp)                                :: qrad,hqr
        real(kind=dp)                                :: carg
        real(kind=dp)                                :: hfunc1
        real(kind=dp), dimension(:,:), intent(in)    :: ax,ay,az
        real(kind=dp)                                :: coef
        real(kind=dp), intent(out)                   :: sig
        real(kind=dp)                                :: pi

        pi = dacos(-1.0_dp)

        rad = (x**2 + y**2 + z**2)**(0.5_dp)
        qrad = q * rad
        sig = 0.0_dp

        plen = l+1
        rlen = m+1
        slen = n+1
        
        do p = 0,plen
            do r = 0,rlen
                do s = 0,slen
                    carg = (l+m+n+p+r+s)/2.0_dp - 1.0_dp
                    beta = ceiling(carg)
                    hqr = q**(beta-1) / rad**beta
                    hfunc1 = hqr * spherical_bessel_jn(beta,qrad)
                    coef = ax(l+1,p+1) * ay(m+1,r+1) * az(n+1,s+1)
                    sig = sig + coef * hfunc1
                enddo
            enddo
        enddo

        sig = 0.5_dp * (3.0_dp/pi)**(0.5_dp) * sig
        
    END SUBROUTINE



END MODULE

