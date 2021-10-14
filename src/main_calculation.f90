!-----------------------------------------------------------------------
! Fortran Modules for PyXSCAT code
! Andres Moreno Carrascosa and Mats Simmermacher, 2021
!-----------------------------------------------------------------------

module main_calculation_mod


    implicit none 

    contains

subroutine total_scattering_calculation(type,state1,state2,maxl,ngto,ng,ga,l,m,n,xx,yy,zz, &
        mmod,q,nq, group,&
        cutoffz,cutoffmd,cutoffcentre,confs,civecs,q_abs,result)




    use onerdm

    use twordms

    use variables

    use integrals

    implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

        integer(kind=ikind), intent(in):: ngto, ng,  nq, maxl,type, state1,state2
        integer(kind=ikind), intent(in),dimension(:) :: l, m, n,group
        integer(kind=ikind), dimension(:,:), intent(in):: confs

        real(kind=dp), intent(in),dimension(:) :: ga, xx, yy, zz, q
        real(kind=dp), intent(in),dimension(:,:) :: mmod, civecs
        real(kind=dp),dimension(:,:,:),allocatable :: z1,z2
        real(kind=dp),dimension(ngto,ngto):: z

        real(kind=dp), intent(in) :: cutoffmd, cutoffz,cutoffcentre

        real(kind=dp), intent(out), dimension(nq):: result,q_abs
        REAL(kind=dp), DIMENSION(size(q),4*maxval(l)+1,4*maxval(l)+1,4*maxval(l)+1) :: P0matrix
         real(kind=dp),  dimension(maxl*2+1,maxl+1,maxl+1,ng,ng) :: ddx,ddy,ddz
        real(kind=dp), dimension(ng,ng) :: px,py,pz
        real(kind=dp),  dimension(:,:,:,:), allocatable :: zcontr
        real(kind=dp),  dimension(:,:), allocatable :: onerdm_matrix
        real(kind=dp),  dimension(:), allocatable :: total,newtotal
        real(kind=dp),  dimension(nq,ngto,ngto) :: e12
        INTEGER(kind=ikind), DIMENSION(maxval(group))   :: group_start, group_count
        integer(kind=ikind), dimension(:), allocatable :: m1, m2, m3, m4
        integer(kind=ikind), dimension(:,:), allocatable :: mat,ep3,ndiff2
        integer(kind=ikind):: nmat,i,j,nmomax,LL,MM,NN
        real(kind=dp) :: start,time1,time2,time3,time4,co,wl,rr
        complex(kind=dp), dimension(:,:,:,:), allocatable :: exponent1, exponent2
        complex(kind=dp), dimension(size(q)):: resultaligned
        real(kind=dp),dimension(3,size(q)) :: q_al
        real(KIND=dp), dimension(size(q)) :: ss

        print*, 'here we are '

       ! allocate(ep3(size(confs(:,1)),size(confs(:,1))),ndiff2(size(confs(:,1)),size(confs(:,1))) )

        do i = 1, Ng
            do j = 1, Ngto
                if (group(j) == i) then
                    group_start(i) = j
                    group_count(i) = count(group==i)
                    exit
                end if
            end do
        end do
        q_abs=0
       ! call onerdm_creat(confs,civecs,onerdm_matrix,nmomax)
        if (type==1) then
            call cpu_time(time2)
            P0matrix = 0
            CALL set_P0(P0matrix, 4*maxval(l), q)
            call maxcoincidence(confs,ep3,ndiff2)
            call createtwordm(confs,civecs,ndiff2,ep3,mat,total)
            allocate(m1(size(mat(:,1))), m2(size(mat(:,1))), m3(size(mat(:,1))), m4(size(mat(:,1))))
             m1 = mat(:,1)
             m2 = mat(:,2)
             m3 = mat(:,3)
             m4 = mat(:,4)
             nmat=size(m1)


            call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q,nq)

            call cpu_time(time3)
            print*,'Time variables', time3-time2

            call tot_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z1,z2,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,result)
            print*,'here'
        else if (type==2) then
               P0matrix = 0
               CALL set_P0(P0matrix, 4*maxval(l), q)
            call maxcoincidence(confs,ep3,ndiff2)
            call onerdm_creat(confs,civecs,onerdm_matrix,nmomax)
               print*,'maxnmo',nmomax
               co=0
               do i=1,nmomax
                   co=co+onerdm_matrix(i,i)
               end do
               print*,'trace',co
            call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,onerdm_matrix,nmomax,q,nq)
            call elastic_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,result)


        else if (type==3) then
            call cpu_time(time2)


            call maxcoincidence(confs,ep3,ndiff2)
            call createtwordm(confs,civecs,ndiff2,ep3,mat,total)
            allocate(m1(size(mat(:,1))), m2(size(mat(:,1))), m3(size(mat(:,1))), m4(size(mat(:,1))))
             m1 = mat(:,1)
             m2 = mat(:,2)
             m3 = mat(:,3)
             m4 = mat(:,4)
             nmat=size(m1)

              allocate(exponent1(nq,maxval(l)*4+1,maxval(l)*4+1,maxval(l)*4+1), &
                    exponent2(nq,maxval(l)*4+1,maxval(l)*4+1,maxval(l)*4+1))

             wl=4.0_dp*dacos(-1.0_dp)/maxval(q)

             call linspace(0.0_dp,dacos(-1.0_dp)-dacos(-1.0_dp)/nq,nq,ss);
            print*,'linspace calculated', wl
            do mm=1,nq
                rr=((2.d0*dacos(-1.00_dp)/(wl)))
                q_al(1,mm)=-rr*sin(ss(mm))*cos(0.0_dp)
                q_al(2,mm)=-rr*sin(ss(mm))*sin(0.0_dp)
                q_al(3,mm)=-rr*cos(ss(mm))+rr

                q_abs(mm)=sqrt(q_al(1,mm)**2+q_al(2,mm)**2+q_al(3,mm)**2)
            end do

            call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q_abs,nq)



            print*,'q_al calculated', q_al(1,:)


            do LL=0,maxval(l)*4
                do MM=0,maxval(m)*4
                    do NN=0,maxval(n)*4

                    exponent2(:,LL+1,MM+1,NN+1)=(complex(0,-1.00)*q_al(1,:))**(LL)* &
                            (complex(0,-1.00)*q_al(2,:))**(MM)*(complex(0,-1.00)*q_al(3,:))**(NN)

                    exponent1(:,LL+1,MM+1,NN+1)=(complex(0,1.00)*q_al(1,:))**(LL) &
                            *(complex(0,1.00)*q_al(2,:))**(MM)*(complex(0,1.00)*q_al(3,:))**(NN)
                     enddo
                enddo
            enddo
            print*,'exponents calculated', (exponent1(1,1,1,1))

            call tot_integration_aligned(ng,px,py,pz,l,m,n,ddx,ddy,ddz,z1,z2,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q_al,e12,exponent1,exponent2,resultaligned)

            result=real(resultaligned)



         else if (type==4) then


             call maxcoincidence(confs,ep3,ndiff2)
            call onerdm_creat(confs,civecs,onerdm_matrix,nmomax)
               print*,'maxnmo',nmomax
             wl=4.0_dp*dacos(-1.0_dp)/maxval(q)

             call linspace(0.0_dp,dacos(-1.0_dp)-dacos(-1.0_dp)/nq,nq,ss);
            print*,'linspace calculated', wl
            do mm=1,nq
                rr=((2.d0*dacos(-1.00_dp)/(wl)))
                q_al(1,mm)=-rr*sin(ss(mm))*cos(0.0_dp)
                q_al(2,mm)=-rr*sin(ss(mm))*sin(0.0_dp)
                q_al(3,mm)=-rr*cos(ss(mm))+rr

                q_abs(mm)=sqrt(q_al(1,mm)**2+q_al(2,mm)**2+q_al(3,mm)**2)
            end do
             co=0
               do i=1,nmomax
                   co=co+onerdm_matrix(i,i)
               end do
               print*,'trace',co
             call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,onerdm_matrix,nmomax,q_abs,nq)
            print*, 'variables calculated'
            allocate(exponent1(nq,maxval(l)*4+1,maxval(l)*4+1,maxval(l)*4+1))
             do LL=0,maxval(l)*4
                do MM=0,maxval(m)*4
                    do NN=0,maxval(n)*4


                    exponent1(:,LL+1,MM+1,NN+1)=(complex(0,1.00)*q_al(1,:))**(LL) &
                            *(complex(0,1.00)*q_al(2,:))**(MM)*(complex(0,1.00)*q_al(3,:))**(NN)
                     enddo
                enddo
            enddo
            print*, 'calling elastic'
             call elastic_integration_alig(ng,px,py,pz,l,m,n,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q_al,e12,exponent1,resultaligned)
                 result=abs(resultaligned)
        end if
        end subroutine total_scattering_calculation

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

        end module main_calculation_mod
