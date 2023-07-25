!-----------------------------------------------------------------------
! Fortran Modules for PyXSCAT code
! Andres Moreno Carrascosa and Mats Simmermacher, 2021
!-----------------------------------------------------------------------

module main_calculation


    implicit none 

    contains

subroutine total_scattering_calculation(type,state1,state2,maxl,ngto,ng,ga,l,m,n,xx,yy,zz, &
        mmod,q,nq, group,&
        cutoffz,cutoffmd,cutoffcentre,confs,civecs,result)




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

        real(kind=dp), intent(out), dimension(nq):: result
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
        integer(kind=ikind):: nmat,i,j,nmomax
        real(kind=dp) :: start,time1,time2,time3,time4,co



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


        end if
        end subroutine total_scattering_calculation


        end module main_calculation
