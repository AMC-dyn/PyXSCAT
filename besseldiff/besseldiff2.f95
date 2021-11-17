MODULE besseldiff2

    implicit none 

    contains



    SUBROUTINE rrdj2a(lmax,x,a)
    ! This subroutine uses a recurrence relation to calculate the
    ! first set of expansion coefficients, a, at given coordinate
    ! x up to an angular momentum of L 

        integer, parameter    :: dp = selected_real_kind(15)
        integer, intent(in)                                        :: lmax
        integer                                                    :: l, p
        real(kind=dp), intent(in)                                  :: x
        real(kind=dp), dimension(:,:), intent(out), allocatable    :: a

        allocate( a(lmax+1,lmax+3) )
        
        a(1,3) = 1
        a(2,4) = -x

        do l = 2, lmax
            a(l+1,3) = -(l-1) * a(l-1,3)
        enddo

        do l = 2, lmax
            do p = 1, l
                a(l+1,p+3) = -(l-1) * a(l-1,p+3) - x * a(l,p+2)
            enddo
        enddo

    END SUBROUTINE



    SUBROUTINE rrdj2b(a,b)
    ! This subroutine uses the first set of of expansion coefficients,
    ! a, to calculate the second set, b. 
    
        integer, parameter    :: dp = selected_real_kind(15)
        integer                                                    :: lmax
        integer                                                    :: l, p
        real(kind=dp), dimension(:,:), intent(in)                  :: a
        real(kind=dp), dimension(:,:), intent(out), allocatable    :: b
        real(kind=dp)                                              :: carg

        lmax = size(a,1)-1

        allocate( b(lmax+1,lmax+3) )
        
        do l = 1, lmax
            do p = 0, l
                carg = (l+p)/2
                b(l+1,p+1) = 2 * ceiling(carg) * a(l+1,p+3)
            enddo
        enddo

    END SUBROUTINE
    
    
    
    SUBROUTINE dj2expand(l,m,n,q,x,y,z,ax,ay,az,bx,by,bz,sig)
    ! This subroutine expands the derivatives for given angular
    ! momenta, L, M, N, coordinates, Hx, Hy, Hz, and momentum transfer,
    ! q. It requires the sets of expansion coefficients, ax, ay, az,
    ! bx, by, bz, calculated at coordinates Hx, Hy, Hz, respectively,
    ! with the two routines above. 
    
        use calculate_form_factors
       
        integer, parameter    :: dp = selected_real_kind(15)
        integer, intent(in)                          :: l,m,n
        integer                                      :: p,r,s
        integer                                      :: beta
        integer                                      :: plen,rlen,slen
        real(kind=dp), intent(in)                    :: q,x,y,z
        real(kind=dp)                                :: rad,eta
        real(kind=dp)                                :: qrad,hqr
        real(kind=dp)                                :: carg
        real(kind=dp)                                :: hfunc2
        real(kind=dp), dimension(:,:), intent(in)    :: ax,ay,az
        real(kind=dp), dimension(:,:), intent(in)    :: bx,by,bz
        real(kind=dp)                                :: coef
        real(kind=dp), intent(out)                   :: sig
        real(kind=dp)                                :: pi

        pi = dacos(-1.0_dp)

        rad = sqrt(x**2 + y**2 + z**2)
        eta = 3.0_dp * z**2 - rad**2
        qrad = q * rad
        sig = 0

        plen = l+2
        rlen = m+2
        slen = n+2
        
        do p = 0,plen
            do r = 0,rlen
                do s = 0,slen
                    carg = (l+m+n+p+r+s)/2-1
                    beta = ceiling(carg)
                    hqr = q**(beta-2) / rad**beta
                    hfunc2 = hqr * spherical_bessel_jn(beta,qrad)
                    coef = ax(l+1,p+1) * ay(m+1,r+1) * az(n+1,s+1)
                    sig = sig + coef * eta * hfunc2
                    coef = bx(l+1,p+1) * ay(m+1,r+1) * az(n+1,s+1)
                    sig = sig + coef * hfunc2
                    coef = ax(l+1,p+1) * by(m+1,r+1) * az(n+1,s+1)
                    sig = sig + coef * hfunc2
                    coef = ax(l+1,p+1) * ay(m+1,r+1) * bz(n+1,s+1)
                    sig = sig + coef * hfunc2
                enddo
            enddo
        enddo

        sig = (1.0_dp/4.0_dp) * sqrt(5.0_dp/pi) * sig
        
    END SUBROUTINE



END MODULE

