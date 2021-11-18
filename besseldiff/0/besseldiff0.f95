MODULE besseldiff0

    implicit none 

    contains



    SUBROUTINE rrdj0(lmax,x,a)
    ! This subroutine uses a recurrence relation to calculate the
    ! first set of expansion coefficients, a, at given coordinate
    ! x (or y or z) up to an angular momentum of L 

        integer, parameter    :: dp = selected_real_kind(15)
        integer, intent(in)                                        :: lmax
        integer                                                    :: l,p
        real(kind=dp), intent(in)                                  :: x
        real(kind=dp), dimension(:,:), intent(out), allocatable    :: a

        allocate( a(lmax+1,lmax+1) )
        
        a(1,1) = 1.0_dp
        a(2,2) = -x

        do l = 2, lmax
            a(l+1,1) = -(l-1) * a(l-1,1)
        enddo

        do l = 2, lmax
            do p = 1, l
                a(l+1,p+1) = -(l-1) * a(l-1,p+1) - x * a(l,p)
            enddo
        enddo

    END SUBROUTINE

 
    
    SUBROUTINE dj0expand(l,m,n,q,x,y,z,ax,ay,az,sig)
    ! This subroutine expands the derivatives for given angular
    ! momenta, L, M, N, coordinates, Hx, Hy, Hz, and momentum transfer,
    ! q. It requires the sets of expansion coefficients, ax, ay, az,
    ! calculated at coordinates Hx, Hy, Hz, respectively,
    ! with the routine above. 
    
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
        real(kind=dp)                                :: hfunc0
        real(kind=dp), dimension(:,:), intent(in)    :: ax,ay,az
        real(kind=dp)                                :: coef
        real(kind=dp), intent(out)                   :: sig
        real(kind=dp)                                :: pi

        pi = dacos(-1.0_dp)

        rad = (x**2 + y**2 + z**2)**(0.5_dp)
        qrad = q * rad
        sig = 0.0_dp

        plen = l
        rlen = m
        slen = n
        
        do p = 0,plen
            do r = 0,rlen
                do s = 0,slen
                    carg = (l+m+n+p+r+s)/2.0_dp
                    beta = ceiling(carg)
                    hqr = (q / rad)**beta
                    hfunc0 = hqr * spherical_bessel_jn(beta,qrad)
                    coef = ax(l+1,p+1) * ay(m+1,r+1) * az(n+1,s+1)
                    sig = sig + coef * hfunc0
                enddo
            enddo
        enddo

        sig = 0.5_dp * (1.0_dp/pi)**(0.5_dp) * sig
        
    END SUBROUTINE



END MODULE

