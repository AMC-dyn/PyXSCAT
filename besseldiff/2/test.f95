     program test

         use besseldiff2
         use calculate_form_factors

         integer, parameter    :: dp = selected_real_kind(15)
         integer                                       :: lmax
         integer                                       :: l,m,n
         integer                                       :: i
         real                                          :: ul,um,un
         real(kind=dp)                                 :: q,x,y,z
         real(kind=dp), dimension(:,:), allocatable    :: ax,ay,az
         real(kind=dp), dimension(:,:), allocatable    :: bx,by,bz
         real(kind=dp)                                 :: sig

         lmax = 6

         call random_number(ul)
         l = floor(lmax * ul)
         write(*,*) "value of L:"
         write(*,*) l
         call random_number(um)
         m = floor(lmax * um)
         write(*,*) "value of M:"
         write(*,*) m
         call random_number(un)
         n = floor(lmax * un)
         write(*,*) "value of N:"
         write(*,*) n

         call random_number(q)
         q = 0.35_dp + 4.65_dp * q
         write(*,*) "value of q:"
         write(*,*) q

         call random_number(x)
         x = 5.0_dp * (x - 0.5_dp)
         write(*,*) "value of x:"
         write(*,*) x
         call random_number(y)
         y = 5.0_dp * (y - 0.5_dp)
         write(*,*) "value of y:"
         write(*,*) y
         call random_number(z)
         z = 5.0_dp * (z - 0.5_dp)
         write(*,*) "value of z:"
         write(*,*) z

         call rrdj2a(lmax,x,ax)
         call rrdj2a(lmax,y,ay)
         call rrdj2a(lmax,z,az)
         
         call rrdj2b(ax,bx)
         call rrdj2b(ay,by)
         call rrdj2b(az,bz)
         bz = -2 * bz
         
         ! write(*,*) "coefficients ax:"
         ! do i = 1, lmax+1
         !     write(*,*) ax(i,1:i+2)
         ! enddo
         ! write(*,*) "coefficients ay:"
         ! do i = 1, lmax+1
         !     write(*,*) ay(i,1:i+2)
         ! enddo
         ! write(*,*) "coefficients az:"
         ! do i = 1, lmax+1
         !     write(*,*) az(i,1:i+2)
         ! enddo
         
         ! write(*,*) "coefficients bx:"
         ! do i = 1, lmax+1
         !     write(*,*) bx(i,1:i+2)
         ! enddo
         ! write(*,*) "coefficients by:"
         ! do i = 1, lmax+1
         !     write(*,*) by(i,1:i+2)
         ! enddo
         write(*,*) "coefficients bz:"
         do i = 1, lmax+1
             write(*,*) bz(i,1:i+2)
         enddo

         call dj2expand(l,m,n,q,x,y,z,ax,ay,az,bx,by,bz,sig)
         
         write(*,*) "value of derivative:"
         write(*,*) sig

     end program

